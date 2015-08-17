#include "uniq_edge_graph.h"

#include <algorithm>
#include <deque>
#include <fstream>
#include <numeric>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("remove_repeats.uniq_edge_graph"));

class ConflictResolver {
public:
    struct Node {
        Node(size_t id, int position, int length) : id(id), position(position), length(length) {
        }
        bool operator < (const Node& o) const {
            return position < o.position;
        }
        size_t id;
        int position;
        int length;
    };
    typedef std::vector< Node > Scaffold;
    typedef std::vector< Scaffold > ScaffoldList;

    ConflictResolver(UniqEdgeGraph* graph) : _graph(graph) {
    }

    void resolve();
    void scaffolds(ScaffoldList& scaffoldlist) const;
private:
    UniqEdgeGraph* _graph;

    friend std::ostream& operator << (std::ostream& os, const ConflictResolver& resoler);
};

void ConflictResolver::scaffolds(ScaffoldList& scaffolds) const {
    LOG4CXX_DEBUG(logger, boost::format("initialize scaffolds"));
	
	scaffolds.clear();

    // BFS
    std::map< size_t, size_t > flag;
    //for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
    for (size_t i = 0; i < _graph->_position_tbl.size(); ++i) {
        if (flag.find(i) == flag.end()) {
            Scaffold component;

            std::deque< size_t > Q = boost::assign::list_of(i);
            while (!Q.empty()) {
                size_t node = Q.front();
                Q.pop_front();
                flag[node] = 2;

                Node ele(node, _graph->_position_tbl[node], _graph->_length_tbl[node]);
                component.push_back(ele);

                UniqEdgeGraph::NodeList::const_iterator k = _graph->_nodelist.find(node);
                if (k != _graph->_nodelist.end()) {
                    for (UniqEdgeGraph::Children::const_iterator j = k->second.children.begin(); j != k->second.children.end(); ++j) {
                        if (flag.find(j->first) == flag.end()) {
                            Q.push_back(j->first);
                            flag[j->first] = 1;
                        }
                    }
                    for (UniqEdgeGraph::Parents::const_iterator j = k->second.parents.begin(); j != k->second.parents.end(); ++j) {
                        if (flag.find(*j) == flag.end()) {
                            Q.push_back(*j);
                            flag[*j] = 1;
                        }
                    }
                }
            }

            scaffolds.push_back(component);
        }
    }
    
    LOG4CXX_DEBUG(logger, boost::format("#### postion=[%d] scaffolds=[%d] flag[%d],_nodelist[%d]") % _graph->_position_tbl.size() % scaffolds.size() % flag.size() % _graph->_nodelist.size());
	
    BOOST_FOREACH(Scaffold& component, scaffolds) {
        std::sort(component.begin(), component.end());
    }
}

void ConflictResolver::resolve() {
    typedef boost::tuple< size_t, size_t, size_t > ConflictItem;
    typedef std::vector< ConflictItem > ConflictList;
    ConflictList overlaplist;

    // Find overlaps
    //
    {
        ScaffoldList scaffolds;
        this->scaffolds(scaffolds);

        int conflict_scaf = 0;

        for (size_t k = 0; k < scaffolds.size(); ++k) {
            if (!scaffolds.empty()) {
                bool find_cft = false;

                for (size_t i = 0; i < scaffolds[k].size() - 1; ++i) {
                    for(size_t j = i + 1; j < scaffolds[k].size(); ++j) {
                        if (scaffolds[k][i].position + scaffolds[k][i].length > scaffolds[k][j].position && 
                                !_graph->hasEdge(scaffolds[k][i].id, scaffolds[k][j].id) &&
                                !_graph->hasEdge(scaffolds[k][j].id, scaffolds[k][i].id)) {
                            if (scaffolds[k][i].position + scaffolds[k][i].length < scaffolds[k][j].position + scaffolds[k][j].length) { // overlap
                                if (_graph->_max_overlap <= scaffolds[k][i].position + scaffolds[k][i].length  - scaffolds[k][j].position) {
                                    overlaplist.push_back(boost::make_tuple(k, scaffolds[k][i].id, scaffolds[k][j].id));
                                    find_cft = true;
                                }

                            } else if (_graph->_max_overlap <= scaffolds[k][j].length) {
                                overlaplist.push_back(boost::make_tuple(k, scaffolds[k][i].id, scaffolds[k][j].id));
                                find_cft = true;
                            }
                        } else {
                            break;
                        }
                    }
                }

                if (find_cft) {
                    conflict_scaf ++;
                }
            }
        }

        LOG4CXX_DEBUG(logger, boost::format("scaffolds containing conflicts = %d") % conflict_scaf);
        LOG4CXX_DEBUG(logger, boost::format("conflict overlap contigs pair number = %d") % overlaplist.size());
    }

    // Resolve
    BOOST_FOREACH(const ConflictItem& overlap, overlaplist) {
		if (_graph->hasEdge(overlap.get< 1 >(), overlap.get< 2 >()) || _graph->hasEdge(overlap.get< 2 >(), overlap.get< 1 >())) {
			continue;
		}

		int temp1 = _graph->getAncestor(overlap.get< 1 >(), overlap.get< 2 >());	
		while (temp1 >= 0) {
            UniqEdgeGraph::EdgeInfo edgeinfo = _graph->getMinEdge(temp1);
			if (edgeinfo.from >= 0) {
				_graph->removeEdge(edgeinfo.from, edgeinfo.to);
                LOG4CXX_DEBUG(logger, boost::format("backward chimeric link %d,%d") % edgeinfo.from % edgeinfo.to); 
			} else {
				break;
			}
			temp1 = _graph->getAncestor(overlap.get< 1 >(), overlap.get< 2 >());
		}
		
		while (temp1 >= 0) {
			LOG4CXX_DEBUG(logger, boost::format("%d\tancestor of\t(%d,%d)") % temp1 % overlap.get< 1 >() % overlap.get< 2 >()); 
			_graph->removeNode(temp1);

			temp1 = _graph->getAncestor(overlap.get< 1 >(), overlap.get< 2 >());
		}
		
		int temp2 = _graph->getDescendant(overlap.get< 1 >(), overlap.get< 2 >());
		while (temp2 >= 0) {
            UniqEdgeGraph::EdgeInfo edgeinfo = _graph->getMinEdge(temp2);
			if (edgeinfo.from >= 0) {
				_graph->removeEdge(edgeinfo.from, edgeinfo.to);
				LOG4CXX_DEBUG(logger, boost::format("forward chimeric link %d,%d") % edgeinfo.from % edgeinfo.to); 
			} else {
				break;
			}
			temp2 = _graph->getDescendant(overlap.get< 1 >(), overlap.get< 2 >());
		}

		while (temp2 >= 0) {
			LOG4CXX_DEBUG(logger, boost::format("%d\tdescendant of\t(%d,%d)") % temp2 % overlap.get< 1 >() % overlap.get< 2 >()); 
			_graph->removeNode(temp2);

			temp2 = _graph->getDescendant(overlap.get< 1 >(), overlap.get< 2 >());
		}
	}
}

bool UniqEdgeGraph::input_inner_component() {
    
    std::string file = boost::str(boost::format("component_%ld") % _iteration);
    std::ifstream stream(file.c_str());
	
	if (!stream) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
        return false;
	}

    ComponentReader reader(stream);
    Component component;
    while (reader.read(component)) {
        _component_tbl.push_back(component);
    }

    return true;
}

bool UniqEdgeGraph::input_edge_position() {

    std::string file = boost::str(boost::format("edge_cluster_pos_%ld") % _iteration);
    std::ifstream stream(file.c_str());
	if (!stream) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
        return false;
	}
	
    std::string line;
	while (std::getline(stream, line)) {
        _position_tbl.push_back(boost::lexical_cast< size_t >(line));
	}
    return true;
}

bool UniqEdgeGraph::input_edge_length() {

    std::string file = boost::str(boost::format("edge_cluster_len_%ld") % _iteration);
    std::ifstream stream(file.c_str());
	if (!stream) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
        return false;
	}

    std::string line;
	while (std::getline(stream, line)) {	
        _length_tbl.push_back(boost::lexical_cast< size_t >(line));
	}
    return true;
}

bool UniqEdgeGraph::input_edge_link(const std::string& name) {
    std::string file = boost::str(boost::format("%s_%ld") % name % _iteration);
    std::ifstream stream(file.c_str());
	
	if (!stream) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
        return false;
	}

 	int from, to, dis, c;
	double score;
    std::string line;
	while (std::getline(stream, line)) {
		sscanf(line.c_str(), "%d\t%d\t%d\t%d\t%lf", &from, &to, &dis, &c, &score);
		addEdge(from, to, dis, c, score);
	}
    return true;
}

void UniqEdgeGraph::linearize() {
    ConflictResolver resovler(this);

    resovler.resolve();

    ConflictResolver::ScaffoldList scaffolds;
    resovler.scaffolds(scaffolds);

    std::string file = boost::str(boost::format("component_%ld") % (_iteration + 1));
    std::ofstream out(file.c_str());

	for (size_t i = 0; i < scaffolds.size(); ++i) {
		out << ">component " << i << std::endl;

		for (size_t j = 0; j < scaffolds[i].size(); ++j) {
			//cout << "line component : " <<  _component[i][j].id << endl;
			for (size_t k = 0; k < _component_tbl[scaffolds[i][j].id].items.size(); ++k) {
				out << _component_tbl[scaffolds[i][j].id].items[k].contig << " ";
			}		
		}
		out << std::endl;
		for (int k = 0; k < _component_tbl[scaffolds[i][0].id].items.size(); k ++)
		{
			out << _component_tbl[scaffolds[i][0].id].items[k].gap << " ";
		}
		for (int j = 1; j < scaffolds[i].size(); j ++)
		{	
			//out << _position_tbl[_component[i][j].id] - _position_tbl[_component[i][j-1].id] - _length_tbl[_component[i][j-1].id] + K << " " ;
			int tmp_dis = getDistance(scaffolds[i][j-1].id, scaffolds[i][j].id);
			if (tmp_dis >= 0)
			{
				out << tmp_dis - scaffolds[i][j-1].length + _K  << " " ;
			}else
			{
				out << _position_tbl[scaffolds[i][j].id] - _position_tbl[scaffolds[i][j-1].id] - _length_tbl[scaffolds[i][j-1].id] + _K << " " ;
			}

			for (int k = 0; k < _component_tbl[scaffolds[i][j].id].items.size(); k ++)
			{
				out << _component_tbl[scaffolds[i][j].id].items[k].gap << " ";
			}
		}
		out << std::endl;
	}

	//initialize_component("new_component");
	//remove_arc_con_edge_from_overlap_pair();	
	
    //resovler.init();
	//initialize_component("new_component");

	//output_graph("contig_arc_graph_after_repeats_removing");	
	//tran_to_line();
}

void UniqEdgeGraph::tran_to_line() {
	//cout << "begin transform to line" << endl;

    std::string file = boost::str(boost::format("component_%ld") % (_iteration + 1));
    std::ofstream out(file.c_str());

	for(int i = 0; i < _component.size(); i ++)
	{
		out << ">component " << i << std::endl;
		 	
		for (int j = 0; j < _component[i].size(); j ++)
		{
			//cout << "line component : " <<  _component[i][j].id << endl;
			for (int k = 0; k < _component_tbl[_component[i][j].id].items.size(); k ++)
			{
				out << _component_tbl[_component[i][j].id].items[k].contig << " ";
			}		
		}
		out << std::endl;
		for (int k = 0; k < _component_tbl[_component[i][0].id].items.size(); k ++)
		{
			out << _component_tbl[_component[i][0].id].items[k].gap << " ";
		}
		for (int j = 1; j < _component[i].size(); j ++)
		{	
			//out << _position_tbl[_component[i][j].id] - _position_tbl[_component[i][j-1].id] - _length_tbl[_component[i][j-1].id] + K << " " ;
			int tmp_dis = getDistance(_component[i][j-1].id, _component[i][j].id);
			if (tmp_dis >= 0)
			{
				out << tmp_dis - _component[i][j-1].len + _K  << " " ;
			}else
			{
				out << _position_tbl[_component[i][j].id] - _position_tbl[_component[i][j-1].id] - _length_tbl[_component[i][j-1].id] + _K << " " ;
			}

			for (int k = 0; k < _component_tbl[_component[i][j].id].items.size(); k ++)
			{
				out << _component_tbl[_component[i][j].id].items[k].gap << " ";
			}
		}
		out << std::endl;
	}
	out.close();
	
    file = boost::str(boost::format("component_cluster_%ld") % _iteration);
	out.open(file.c_str());
	
	for(int i = 0; i < _component.size(); i ++)
	{
		out << ">component " << i << std::endl;
		
		for (int j = 0; j < _component[i].size(); j ++)
		{
			out << _component[i][j].id << "|" << _length_tbl[_component[i][j].id] << "|" << _position_tbl[_component[i][j].id] << "\t";
		}
		out << std::endl;
		for (int j = 1; j < _component[i].size(); j ++)
		{	
			int tmp_dis = getDistance(_component[i][j-1].id, _component[i][j].id);
			if (tmp_dis >= 0)
			{
				out << tmp_dis - _component[i][j-1].len + _K  << " " ;
			}else
			{
				out << _position_tbl[_component[i][j].id] - _position_tbl[_component[i][j-1].id] - _length_tbl[_component[i][j-1].id] + _K << " " ;
			}
		}
		out << std::endl;
	}
	out.close();
}

int cmp( Edge_Seq_Element es1, Edge_Seq_Element es2)
{
	return es1.pos < es2.pos;
}

void UniqEdgeGraph::initialize_component(const std::string& cmp_name) {
    std::cout << "initialize scaffolds" << std::endl;
	
	_component.clear();

    // BFS
    std::map< size_t, size_t > flag;
    //for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
    for (size_t i = 0; i < _position_tbl.size(); ++i) {
        if (flag.find(i) == flag.end()) {
            std::vector< Edge_Seq_Element > elements;

            std::deque< size_t > Q = boost::assign::list_of(i);
            while (!Q.empty()) {
                size_t node = Q.front();
                Q.pop_front();
                flag[node] = 2;

                Edge_Seq_Element ele(node, _position_tbl[node], _length_tbl[node]);
                elements.push_back(ele);

                NodeList::const_iterator k = _nodelist.find(node);
                if (k != _nodelist.end()) {
                    for (Children::const_iterator j = k->second.children.begin(); j != k->second.children.end(); ++j) {
                        if (flag.find(j->first) == flag.end()) {
                            Q.push_back(j->first);
                            flag[j->first] = 1;
                        }
                    }
                    for (Parents::const_iterator j = k->second.parents.begin(); j != k->second.parents.end(); ++j) {
                        if (flag.find(*j) == flag.end()) {
                            Q.push_back(*j);
                            flag[*j] = 1;
                        }
                    }
                }
            }

            _component.push_back(elements);
        }
    }
    
    LOG4CXX_DEBUG(logger, boost::format("#### postion=[%d] componets=[%d] flag[%d],_nodelist[%d]") % _position_tbl.size() % _component.size() % flag.size() % _nodelist.size());
	
	for (int i = 0; i < _component.size(); i++) {
		sort(_component[i].begin(), _component[i].end(), cmp);
	}

	LOG4CXX_DEBUG(logger, boost::format("\tscaffolds number = %d") % _component.size());

}

int UniqEdgeGraph::EdgeScorePlus(int l, const EdgeInfo& r) const {
    return l + r.edge.score;
}

UniqEdgeGraph::EdgeInfo UniqEdgeGraph::getMinEdge(size_t index) const {
    std::list< EdgeInfo > edgelist;
	EdgeInfo edgeinfo;
	 
    // Children
    {
        NodeList::const_iterator i = _nodelist.find(index);
        if (i != _nodelist.end()) {
            for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                edgeinfo.from = index;
                edgeinfo.to = j->first;
                edgeinfo.edge = j->second;
                edgelist.push_back(edgeinfo);
            }
        }
    }
    // Parents
    {
        NodeList::const_iterator i = _nodelist.find(index);
        if (i != _nodelist.end()) {
            for (Parents::const_iterator j = i->second.parents.begin(); j != i->second.parents.end(); ++j) {
                NodeList::const_iterator k = _nodelist.find(*j);
                BOOST_ASSERT(k != _nodelist.end());
                BOOST_ASSERT(k->second.children.find(i->first) != k->second.children.end());
                edgeinfo.from = *j;
                edgeinfo.to = index;
                Children::const_iterator l = k->second.children.find(i->first);
                edgeinfo.edge = l->second;
                edgelist.push_back(edgeinfo);

            }
        }
    }

    //std::sort(edgelist.begin(), edgelist.end());
    edgelist.sort();
	
	if (edgelist.size() < 2) {
		edgeinfo.from = -1;
		edgeinfo.to = -1;
		return edgeinfo;
	}
	edgeinfo = edgelist.front();
    edgelist.pop_front();

    int avg = 0;
    if (!edgelist.empty()) {
        avg = std::accumulate(edgelist.begin(), edgelist.end(), 0, boost::bind(&UniqEdgeGraph::EdgeScorePlus, this, _1, _2)) / edgelist.size();
    }

	if (edgeinfo.edge.score <  avg) {
		return edgeinfo;
	}

    edgeinfo.from = -1;
    edgeinfo.to = -1;
    return edgeinfo;
}

void UniqEdgeGraph::remove_arc_con_edge_from_overlap_pair() {
    typedef boost::tuple< size_t, size_t > OverlapItem;
    typedef std::vector< OverlapItem > OverlapList;
    OverlapList overlaplist;

    size_t conflict_scaf = 0;

	for (size_t index = 0; index < _component.size(); ++index) {
        bool find_cft = false;

		for (size_t i = 0; i < _component[index].size() - 1; ++i) {
			for(size_t j = i + 1; j < _component[index].size(); ++j) {
				if (_component[index][i].pos + _component[index][i].len > _component[index][j].pos && 
                        !hasEdge(_component[index][i].id, _component[index][j].id) && 
                        !hasEdge(_component[index][j].id, _component[index][i].id)) {
					if (_component[index][i].pos + _component[index][i].len < _component[index][j].pos + _component[index][j].len ) {
						if (_max_overlap <= _component[index][i].pos + _component[index][i].len	- _component[index][j].pos) {
							overlaplist.push_back(boost::make_tuple(_component[index][i].id, _component[index][j].id));
                            find_cft = true;
						}
					} else if (_max_overlap <= _component[index][j].len) {
                        overlaplist.push_back(boost::make_tuple(_component[index][i].id, _component[index][j].id));
                        find_cft = true;
					}
				} else {
					break;
				}
			}
		}

        if (find_cft) {
            conflict_scaf ++;
        }
	}

    LOG4CXX_DEBUG(logger, boost::format("\tscaffolds containing conflicts = %d") % conflict_scaf);
	LOG4CXX_DEBUG(logger, boost::format("\tconflict overlap contigs pair number = %d") % overlaplist.size());

    BOOST_FOREACH(const OverlapItem& overlap, overlaplist) {
		if (hasEdge(overlap.get< 0 >(), overlap.get< 1 >()) || hasEdge(overlap.get< 1 >(), overlap.get< 0 >())) {
			continue;
		}

		int temp1 = getAncestor(overlap.get< 0 >(), overlap.get< 1 >());	
		while (temp1 >= 0) {
			EdgeInfo edgeinfo = getMinEdge(temp1);
			if (edgeinfo.from >= 0) {
				removeEdge(edgeinfo.from, edgeinfo.to);
				LOG4CXX_DEBUG(logger, boost::format("#%d backward chimeric link %d,%d") % _iteration % edgeinfo.from % edgeinfo.to); 
			} else {
				break;
			}
			temp1 = getAncestor(overlap.get< 0 >(), overlap.get< 1 >());
		}
		
		while (temp1 >= 0) {
			LOG4CXX_DEBUG(logger,  boost::format("#%d %d\tancestor of\t(%d,%d)") % _iteration % temp1 % overlap.get< 0 >() % overlap.get< 1 >()); 
			removeNode(temp1);
			temp1 = getAncestor(overlap.get< 0 >(), overlap.get< 1 >());	
		}
		
		int temp2 = getDescendant(overlap.get< 0 >(), overlap.get< 1 >());
		while (temp2 >= 0) {
			EdgeInfo edgeinfo = getMinEdge(temp2);
			if (edgeinfo.from >= 0) {
				//cout << "del " << edgeinfo.from << "\t" << edgeinfo.to << endl;
				removeEdge(edgeinfo.from, edgeinfo.to);
				LOG4CXX_DEBUG(logger, boost::format("#%d forward chimeric link %d,%d") % _iteration % edgeinfo.from % edgeinfo.to); 
			} else {
				break;
			}
			temp2 = getDescendant(overlap.get< 0 >(), overlap.get< 1 >());
		}

		while (temp2 >= 0) {
			LOG4CXX_DEBUG(logger, boost::format("#%d %d\tdescendant of\t(%d,%d)") % _iteration % temp2 % overlap.get< 0 >() % overlap.get< 1 >()); 
			removeNode(temp2);
			temp2 = getDescendant(overlap.get< 0 >(), overlap.get< 1 >());
		}
	}
}

int UniqEdgeGraph::getAncestor(size_t x, size_t y) const {
    std::map< size_t, int > flag;

    // BFS
    // node=>x
    {
        int depth = 2, width = 0, count = 1;
        std::deque< size_t > Q = boost::assign::list_of(x);
        while (!Q.empty()) {
            size_t node = Q.front();
            Q.pop_front();

            flag[node] = depth;
            NodeList::const_iterator k = _nodelist.find(node);
            if (k != _nodelist.end()) {
                for (Parents::const_iterator i = k->second.parents.begin(); i != k->second.parents.end(); ++i) {
                    if (flag.find(*i) == flag.end()) {
                        Q.push_back(*i);
                        flag[*i] = 1;
                        ++width;
                    }
                }
            }
            if (--count == 0) {
                count = width;
                ++depth;
            }
        }
    }
    // node=>y
    {
        int depth = -1, width = 0, count = 1;
        size_t id = -1;
        int min_depth = 0;
        std::deque< size_t > Q = boost::assign::list_of(y);
        while (!Q.empty()) {
            size_t node = Q.front();
            Q.pop_front();

            if (flag.find(node) != flag.end() && flag[node] >= 2) {
                return node;
            }

            flag[node] = depth;
            NodeList::const_iterator k = _nodelist.find(node);
            if (k != _nodelist.end()) {
                for (Parents::const_iterator i = k->second.parents.begin(); i != k->second.parents.end(); ++i) {
                    if (flag.find(*i) == flag.end()) {
                        Q.push_back(*i);
                        flag[*i] = 1;
                        ++width;
                    } else if (flag[*i] >= 2) {
                        if (id == -1 || flag[*i] < min_depth) {
                            id = *i;
                            min_depth = flag[*i];
                        }
                    }
                }
            }

            if (--count == 0) {
                if (id != -1) {
                    return id;
                }
                count = width;
                --depth;
            }
        }
    }

	return -1;
}

int UniqEdgeGraph::getDescendant(size_t x, size_t y) const {
    std::map< size_t, int > flag;

    // BFS
    // node=>x
    {
        int depth = 2, width = 0, count = 1;
        std::deque< size_t > Q = boost::assign::list_of(x);
        while (!Q.empty()) {
            size_t node = Q.front();
            Q.pop_front();
            flag[node] = 2;

            NodeList::const_iterator k = _nodelist.find(node);
            if (k != _nodelist.end()) {
                for (Children::const_iterator i = k->second.children.begin(); i != k->second.children.end(); ++i) {
                    if(flag[i->first] == 0) {
                        flag[i->first] = 1;
                        Q.push_back(i->first);
                        width++;
                    }
                }
            }

            if (--count == 0) {
                count = width;
                depth++;
            }
        }
    }
    // node=>y
    {
        int depth = -1, width = 0, count = 1;
        std::deque< size_t > Q = boost::assign::list_of(y);
        size_t id = 0;
        int min = 0;
        while (!Q.empty()) {
            size_t node = Q.front();
            Q.pop_front();
            if (flag[node] >= 2) {
                return node;
            }
            flag[node] = 3;

            NodeList::const_iterator k = _nodelist.find(node);
            if (k != _nodelist.end()) {
                for (Children::const_reverse_iterator i = k->second.children.rbegin(); i != k->second.children.rend(); ++i) {
                    if (flag[i->first] == 0) {
                        flag[i->first] = 1;
                        Q.push_back(i->first);
                        width++;
                    } else if (flag[i->first] >= 2) {
                        if (id == 0 || flag[i->first] > min) {
                            id = i->first;
                            min = flag[id];
                        }
                    }
                }
            }

            if (--count == 0) {
                if (id != 0)
                    return id;
                count = width;
                depth--;
            }
        }
    }
	return -1;
}

void UniqEdgeGraph::output_graph(const std::string& filename) const {
	//cout << "begin output edge graph" << endl;
    std::string file = boost::str(boost::format("%s_%ld") % filename % _iteration);
    std::ofstream out(file.c_str());

    for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
        for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
				out << i->first << "|" << _length_tbl[i->first] << "|" << _position_tbl[i->first] << "\t" << j->first << "|" << _length_tbl[j->first] << "|" << _position_tbl[j->first] << "\t" << j->second.distance << "|" << j->second.count << "|" << j->second.score << std::endl;
        }
    }
	out.close();
	//cout << "end out put graph" << endl;
}

void UniqEdgeGraph::addEdge(size_t from, size_t to, int distance, int count, int score) {
    // Add nodes
    {
        if (_nodelist.find(from) == _nodelist.end()) {
           _nodelist[from] = Node();
        }
        if (_nodelist.find(to) == _nodelist.end()) {
           _nodelist[to] = Node();
        }
    }
    // Set Children
    {
        Node& node = _nodelist[from];
        if (node.children.find(to) == node.children.end()) {
            node.children[to] = Edge(distance, count, score);
        }
    }
    // Set Parents
    {
        Node& node = _nodelist[to];
        node.parents.insert(from);
    }

    LOG4CXX_TRACE(logger, boost::format("addEdge: %d,%d,%d,%d,%d") % from % to % distance % _nodelist[from].children.size() % _nodelist[to].parents.size());
}

void UniqEdgeGraph::removeEdge(size_t from, size_t to) {
    LOG4CXX_TRACE(logger, boost::format("removeEdge: from=%d,to=%d") % from % to);
    // Remove child
    {
        NodeList::iterator i = _nodelist.find(from);
        if (i != _nodelist.end()) {
            Children::iterator j = i->second.children.find(to);
            if (j != i->second.children.end()) {
                i->second.children.erase(j);
            }
        }
    }
    // Remove parent
    {
        NodeList::iterator i = _nodelist.find(to);
        if (i != _nodelist.end()) {
            Parents::iterator j = i->second.parents.find(from);
            if (j != i->second.parents.end()) {
                i->second.parents.erase(j);
            }
        }
    }
    // Clear
    {
        NodeList::iterator i = _nodelist.find(from);
        if (i != _nodelist.end() && !i->second) {
            _nodelist.erase(i);
        }
        NodeList::iterator j = _nodelist.find(to);
        if (j != _nodelist.end() && !j->second) {
            _nodelist.erase(j);
        }
    }
}

void UniqEdgeGraph::removeNode(size_t node) {
    LOG4CXX_TRACE(logger, boost::format("removeNode: %d") % node);

    NodeList::iterator i = _nodelist.find(node);
    if (i != _nodelist.end()) {
        std::vector< size_t > checklist;
        for (Children::iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
            NodeList::iterator k = _nodelist.find(j->first);
            if (k != _nodelist.end()) {
                k->second.parents.erase(i->first);
                checklist.push_back(k->first);
            }
        }
        for (Parents::iterator j = i->second.parents.begin(); j != i->second.parents.end(); ++j) {
            NodeList::iterator k = _nodelist.find(*j);
            if (k != _nodelist.end()) {
                k->second.children.erase(i->first);
                checklist.push_back(k->first);
            }

        }

        _nodelist.erase(i);

        BOOST_FOREACH(size_t id, checklist) {
            NodeList::iterator k = _nodelist.find(id);
            if (k != _nodelist.end() && !k->second) {
                _nodelist.erase(k);
            }
        }
    }
}

bool UniqEdgeGraph::hasEdge(size_t from, size_t to) const {
    NodeList::const_iterator i = _nodelist.find(from);
    if (i != _nodelist.end()) {
        return i->second.children.find(to) != i->second.children.end();
    }
    return false;
}

int UniqEdgeGraph::getDistance(size_t from, size_t to) const {
    NodeList::const_iterator i = _nodelist.find(from);
    if (i != _nodelist.end()) {
        Children::const_iterator j = i->second.children.find(to);
        if (j != i->second.children.end()) {
            return j->second.distance;
        }
    }
    return -1;
}
