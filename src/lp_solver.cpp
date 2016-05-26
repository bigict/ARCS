#include "lp_solver.h"

#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <glpk.h>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.LPSolver"));

typedef std::map< size_t, int > Children;
typedef std::set< size_t > Parents;
struct Node {
    Children children;
    Parents parents;
    size_t indegree() const {
        return parents.size();
    }
    size_t outdegree() const {
        return children.size();
    }
};
typedef std::map< size_t, Node > Graph;
typedef std::vector< Graph > GraphList;
typedef std::vector< int > PositionList;

static void Graph_addEdge(Graph& graph, size_t xi, size_t xj, long val) {
    // xi
    {
        Graph::iterator i = graph.find(xi);
        if (i == graph.end()) {
            graph[xi] = Node();
        }
    }
    // xj
    {
        Graph::iterator j = graph.find(xj);
        if (j == graph.end()) {
            graph[xj] = Node();
        }
    }
    // children
    {
        Graph::iterator i = graph.find(xi);
        BOOST_ASSERT(i != graph.end());
        i->second.children[xj] = val;
    }
    // parents
    {
        Graph::iterator j = graph.find(xj);
        BOOST_ASSERT(j != graph.end());
        j->second.parents.insert(xi);
    }
}

static bool Graph_read(Graph& graph, std::istream& stream) {
    if (!stream) {
        return false;
    }

    // s.t. con3122 : x_2250 - x_2866 + e_2866_2250 = 190;
    static boost::regex reg("s\\.t\\. con(\\d+)\\s*:\\s*x_(\\d+)\\s*-\\s*x_(\\d+)\\s*\\+\\s*e_(\\d+)_(\\d+)\\s*=\\s*((\\+|-)?\\d+)\\s*;");

    std::string line;
	while (std::getline(stream, line)) {
        boost::algorithm::trim(line);
        boost::smatch what;
        if (boost::regex_match(line, what, reg)) {
            size_t xj = boost::lexical_cast< size_t >(what[2]);
            size_t xi = boost::lexical_cast< size_t >(what[3]);
            int err = boost::lexical_cast< int >(what[6]);

            LOG4CXX_TRACE(logger, boost::format("regex: %d %d %d %s") % xi % xj % err % line);

            Graph_addEdge(graph, xi, xj, err);
		}	
	}

    return true;
}

static bool Graph_read(Graph& graph, const std::string& file) {
    std::ifstream stream(file.c_str());
    return Graph_read(graph, stream);
}

static void Graph_compress(const Graph& graph, Graph& new_graph) {
    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        if (i->second.indegree() == 1 && i->second.outdegree() == 1) {
            continue;
        }

        for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            int delta = j->second;
            Graph::const_iterator k = graph.find(j->first);
            while (k != graph.end() && k->second.indegree() == 1 && k->second.outdegree() == 1) {
                Children::const_iterator c = k->second.children.begin();
                if (c->first == i->first) { // ring
                    break;
                }
                delta += c->second;
                k = graph.find(c->first);
            }
            BOOST_ASSERT(k != graph.end());

            //std::cout << boost::format("compressed\t%d\t%d\t%d") % i->first % k->first % delta << std::endl;
            Graph_addEdge(new_graph, i->first, k->first, delta);
        }
    }
}

static bool PositionList_write(const PositionList& position_tbl, std::ostream& stream) {
    if (!stream) {
        return false;
    }
    BOOST_FOREACH(int x, position_tbl) {
        stream << x << "\n";
    }
    return true;
}

static bool Graph_solve(Graph& graph, size_t loops, PositionList* position_tbl) {
    glp_prob* lp = glp_create_prob();
    glp_set_prob_name(lp, "scaffold");
    glp_set_obj_dir(lp, GLP_MIN);

    size_t rows = 0, cols = 0, vals = 0;
    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        ++cols; // var x_i
        for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            rows += 3;
            cols += 2; // var e_i_j; var E_i_j;
            vals += 7;
        }
    }

    glp_add_rows(lp, rows);
    glp_add_cols(lp, cols);

    std::map< std::pair< size_t, size_t >, size_t > mapping;
    {
        size_t row = 1, col = 1;
        for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
            mapping[std::make_pair(i->first, -1)] = col;

            glp_set_col_bnds(lp, col++, GLP_LO, 0.0, 0.0); // x_i >= 0
            for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                glp_set_row_bnds(lp, row++, GLP_FX, j->second, j->second);
                glp_set_row_bnds(lp, row++, GLP_LO, 0.0, 0.0);
                glp_set_row_bnds(lp, row++, GLP_LO, 0.0, 0.0);

                mapping[std::make_pair(i->first, j->first)] = col;

                ++col; // var e_i
                glp_set_obj_coef(lp, col++, 1.0);
            }
        }
    }

    LOG4CXX_TRACE(logger, boost::format(" rows = %d, cols = %d, vals = %d") % rows % cols % vals);

    int* ia = new int[vals + 1];
    int* ja = new int[vals + 1];
    double* ra = new double[vals + 1];

    {
        size_t l = 1, row = 1, col = 1;
        for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
            for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                // x_j - x_i + e_i_j = d_i_j
                ia[l] = row;
                ja[l] = mapping[std::make_pair(j->first, -1)];
                ra[l] = 1.0;
                ++l;

                ia[l] = row;
                ja[l] = mapping[std::make_pair(i->first, -1)];
                ra[l] = -1.0;
                ++l;

                ia[l] = row;
                ja[l] = mapping[std::make_pair(i->first, j->first)] + 0;
                ra[l] = 1.0;
                ++l;

                ++row;

                // E_i_j + e_i_j >= 0
                ia[l] = row;
                ja[l] = mapping[std::make_pair(i->first, j->first)] + 1;
                ra[l] = 1.0;
                ++l;

                ia[l] = row;
                ja[l] = mapping[std::make_pair(i->first, j->first)] + 0;
                ra[l] = 1.0;
                ++l;

                ++row;

                // E_i_j - e_i_j >= 0
                ia[l] = row;
                ja[l] = mapping[std::make_pair(i->first, j->first)] + 1;
                ra[l] = 1.0;
                ++l;

                ia[l] = row;
                ja[l] = mapping[std::make_pair(i->first, j->first)] + 0;
                ra[l] = -1.0;
                ++l;

                ++row;
            }
        }
    }

    glp_load_matrix(lp, vals, ia, ja, ra);

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.it_lim = loops;
    //parm.pricing = GLP_PT_PSE;
    //parm.presolve = GLP_ON;
    parm.msg_lev = GLP_MSG_ERR;
    glp_simplex(lp, &parm);

    double z = glp_get_obj_val(lp);
    LOG4CXX_TRACE(logger, boost::format("z = %f") % z);

    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        size_t val = glp_get_col_prim(lp, mapping[std::make_pair(i->first, -1)]);
        if (position_tbl != NULL) {
            (*position_tbl)[i->first] = val;
        }
        LOG4CXX_TRACE(logger, boost::format("x\t%d\t%d") % i->first % val);
    }

    delete[] ra;
    delete[] ja;
    delete[] ia;

    glp_delete_prob(lp);
    return true;
}

static bool Graph_divide(Graph& graph, size_t loops, PositionList* position_tbl) {
    typedef std::set< size_t > NodeList;
    NodeList nodes;

    // nodes
    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        nodes.insert(i->first);
    }

    while (!nodes.empty()) {
        // BFS
        Graph component;
        std::deque< size_t > Q = boost::assign::list_of(*nodes.begin());
        while (!Q.empty()) {
            size_t xi = Q.front();
            Q.pop_front();

            if (nodes.find(xi) == nodes.end()) {
                continue;
            }

            nodes.erase(xi);

            Graph::const_iterator i = graph.find(xi);
            if (i != graph.end()) {
                for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                    Graph_addEdge(component, xi, j->first, j->second);
                    Q.push_back(j->first);
                }
                for (Parents::const_iterator j = i->second.parents.begin(); j != i->second.parents.end(); ++j) {
                    Q.push_back(*j);
                }
            }
        }

        LOG4CXX_TRACE(logger, boost::format("component:%d/%d") % component.size() % graph.size());

        if (!Graph_solve(component, loops, position_tbl)) {
            LOG4CXX_ERROR(logger, "solve component failed");
            return false;
        }
    }

    return true;
}
 

int LPSolver::run(const Properties& options, const Arguments& arguments) {
    int r = 0;

    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    LOG4CXX_INFO(logger, "solverLP begin");

    size_t loops = options.get< size_t >("l", INT_MAX);
    // input & output opened
    std::istream* is = &std::cin;
    if (options.find("i") != options.not_found()) {
        std::string file = options.get< std::string >("i");
        is = new std::ifstream(file.c_str());

        LOG4CXX_DEBUG(logger, boost::format("input: %s") % file);
    }
    if (is == &std::cin) {
        std::cin.sync_with_stdio(false);
    }
    std::ostream* os = &std::cout;
    if (options.find("o") != options.not_found()) {
        std::string file = options.get< std::string >("o");
        os = new std::ofstream(file.c_str());

        LOG4CXX_DEBUG(logger, boost::format("output: %s") % file);
    }

    Graph graph;

	// read graph data
    if (r == 0 && !Graph_read(graph, *is)) {
        LOG4CXX_ERROR(logger, boost::format("failed to read graph from stream!!!"));
        r = 1;
    }

    PositionList position_tbl(options.get< size_t >("EDGE_CLUSTER_NUM"), 0);
    // divide into connected components and solve them one by one.
    if (r == 0 && !Graph_divide(graph, loops, &position_tbl)) {
        LOG4CXX_ERROR(logger, boost::format("solve LP failed!!!"));
        r = 2;
    }
    // write graph data
    if (r == 0 && !PositionList_write(position_tbl, *os)) {
        LOG4CXX_ERROR(logger, boost::format("failed to write position to stream!!!"));
        r = 3;
    }

    // input & output closed
    if (os != &std::cout) {
        delete os;
    }
    if (is != &std::cin) {
        delete is;
    }

    LOG4CXX_INFO(logger, "solveLP end");
    return r;
}

LPSolver LPSolver::_runner;

LPSolver::LPSolver() : Runner("c:s:i:o:l:h") {
    RUNNER_INSTALL("solveLP", this, "solveLP");
}

int LPSolver::printHelps() const {
    std::cout << "usage: arcs solveLP [arguments]" << std::endl;
    std::cout << std::endl;
    std::cout << "\t-c[=<file>]    log config file, default ./log4cxx.properties" << std::endl;
    std::cout << "\t-s[=<file>]    <scaffold_parameter_\%d>" << std::endl;
    std::cout << "\t-i[=<file>]    lp problem file <position_lp_\%d.math>" << std::endl;
    std::cout << "\t-o[=<file>]    output file <edge_cluster_pos_\%d>" << std::endl;
    std::cout << "\t-l[=<number>]  lp prarameter, default INT_MAX(recommend)" << std::endl;
    std::cout << "\t-h             help" << std::endl;
    std::cout << std::endl;
    return 256;
}

int LPSolver::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
