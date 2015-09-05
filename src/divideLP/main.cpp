//
// depart mixture LP to many one small LP
//
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("dividelP.main"));

typedef std::map< size_t, int > Children;
typedef std::set< size_t > Parents;
struct Node {
    Children children;
    Parents parents;
};
typedef std::map< size_t, Node > Graph;
typedef std::vector< Graph > GraphList;

void Graph_addEdge(Graph& graph, size_t xi, size_t xj, long val) {
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

bool Graph_read(Graph& graph, std::istream& stream) {
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

bool Graph_read(Graph& graph, const std::string& file) {
    std::ifstream stream(file.c_str());
    return Graph_read(graph, stream);
}

bool Graph_write(const Graph& graph, std::ostream& stream) {
    if (!stream) {
        return false;
    }

    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        stream << boost::format("var x_%d;") % i->first << std::endl;
    }
    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            stream << boost::format("var e_%d_%d;") % i->first % j->first << std::endl;
            stream << boost::format("var E_%d_%d;") % i->first % j->first << std::endl;
        }
    }

    stream << std::endl;
    stream << "minimize z:   ";
    size_t count = 0;
    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            stream << boost::format(" E_%d_%d + ") % i->first % j->first;
            ++count;
		    if (count % 10 == 0) {
                stream << std::endl;
            }
        }
    }
    stream << "0;" << std::endl << std::endl;

    size_t index = 1;
    for (Graph::const_iterator i = graph.begin(); i != graph.end(); ++i) {
        for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            stream << boost::format("s.t. con%d : x_%d - x_%d + e_%d_%d = %d;") % index++ % j->first % i->first % i->first % j->first % j->second << std::endl;
            stream << boost::format("s.t. con%d : E_%d_%d + e_%d_%d >= 0;") % index++ % i->first % j->first % i->first % j->first << std::endl;
            stream << boost::format("s.t. con%d : E_%d_%d - e_%d_%d >= 0;") % index++ % i->first % j->first % i->first % j->first << std::endl;
        }
    }

    stream << std::endl;
    stream << "end;";

    return true;
}

bool Graph_write(Graph& graph, const std::string& file) {
    std::ofstream stream(file.c_str());
    return Graph_write(graph, stream);
}

void Graph_divide(Graph& graph, GraphList& components) {
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

        components.push_back(component);
    }
}
 
int main(int argc, char **argv) {
	if (argc != 3) {
        std::cerr << boost::format("Usage: %s [big_lp_file] [small_lp_dir]") % argv[0] << std::endl;
        return 1;
	}
    log4cxx::BasicConfigurator::configure();

    boost::filesystem::path workdir = std::string(argv[2]);
    if (!boost::filesystem::exists(workdir) && !boost::filesystem::create_directory(workdir)) {
        LOG4CXX_ERROR(logger, boost::format("failed to create directory: %s") % workdir);
        return 1;
    }

    Graph graph;

	// read graph data
    if (!Graph_read(graph, argv[1])) {
        LOG4CXX_ERROR(logger, boost::format("failed to read graph from file: %s") % argv[1]);
        return 1;
    }
    // divide into connected components
    GraphList components;
    Graph_divide(graph, components);
    // write graph data
    for (size_t i = 0; i < components.size(); ++i) {
        const Graph& component = components[i];
        boost::filesystem::ofstream stream(workdir / boost::str(boost::format("lp%d.math") % i));
        if (!Graph_write(component, stream)) {
            LOG4CXX_ERROR(logger, boost::format("failed to write graph to file: lp%d.math") % i);
            return 1;
        }
    }

	return 0;
}
