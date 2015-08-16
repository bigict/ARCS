#include "uniq_edge_graph.h"

#include <deque>
#include <fstream>
#include <algorithm>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("remove_repeats.uniq_edge_graph"));

bool UniqEdgeGraph::input_inner_component() {
    
    std::string file = boost::str(boost::format("component_%ld") % _iteration);
	ifstream stream(file.c_str());
	
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
	ifstream stream(file.c_str());
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
	ifstream stream(file.c_str());
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
	ifstream stream(file.c_str());
	
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
	initialize_component("new_component");
	remove_arc_con_edge_from_overlap_pair();	
	
	initialize_component("new_component");

	output_graph("contig_arc_graph_after_repeats_removing");	
	tran_to_line();
}

void UniqEdgeGraph::tran_to_line()
{
	//cout << "begin transform to line" << endl;

	line_component = _component;
	
    std::vector< int > flag(_position_tbl.size(), 0);

	for (int i = 0; i < line_component.size(); i ++)
	{
		for (int j = 0; j < line_component[i].size(); j ++)
		{
			flag[line_component[i][j].id] = 2;
		}
	}
	int count = 0;
	for (int i = 0; i < flag.size(); i++)
	{
		if(flag[i] == 0)
		{
			count ++;
			removeNode(i);
		}
	}
	//cout << "\ttran to line delete " << count << " contigs " << endl;
	//fout.close();
	
	//cout << "output line component" << endl;
	//cout << "inner component size : " << inner_component.size() << endl;
	//cout << "line component size : " << line_component.size() << endl;

    std::string file = boost::str(boost::format("component_%ld") % (_iteration + 1));
    ofstream out(file.c_str());

	for(int i = 0; i < line_component.size(); i ++)
	{
		out << ">component " << i << endl;
		 	
		for (int j = 0; j < line_component[i].size(); j ++)
		{
			//cout << "line component : " <<  line_component[i][j].id << endl;
			for (int k = 0; k < _component_tbl[line_component[i][j].id].items.size(); k ++)
			{
				out << _component_tbl[line_component[i][j].id].items[k].contig << " ";
			}		
		}
		out << endl;
		for (int k = 0; k < _component_tbl[line_component[i][0].id].items.size(); k ++)
		{
			out << _component_tbl[line_component[i][0].id].items[k].gap << " ";
		}
		for (int j = 1; j < line_component[i].size(); j ++)
		{	
			//out << _position_tbl[line_component[i][j].id] - _position_tbl[line_component[i][j-1].id] - _length_tbl[line_component[i][j-1].id] + K << " " ;
			int tmp_dis = getDistance(line_component[i][j-1].id, line_component[i][j].id);
			if (tmp_dis >= 0)
			{
				out << tmp_dis - line_component[i][j-1].len + _K  << " " ;
			}else
			{
				out << _position_tbl[line_component[i][j].id] - _position_tbl[line_component[i][j-1].id] - _length_tbl[line_component[i][j-1].id] + _K << " " ;
			}

			for (int k = 0; k < _component_tbl[line_component[i][j].id].items.size(); k ++)
			{
				out << _component_tbl[line_component[i][j].id].items[k].gap << " ";
			}
		}
		out << endl;
	}
	out.close();
	
    file = boost::str(boost::format("component_cluster_%ld") % _iteration);
	out.open(file.c_str());
	
	for(int i = 0; i < line_component.size(); i ++)
	{
		out << ">component " << i << endl;
		
		for (int j = 0; j < line_component[i].size(); j ++)
		{
			out << line_component[i][j].id << "|" << _length_tbl[line_component[i][j].id] << "|" << _position_tbl[line_component[i][j].id] << "\t";
		}
		out << endl;
		for (int j = 1; j < line_component[i].size(); j ++)
		{	
			int tmp_dis = getDistance(line_component[i][j-1].id, line_component[i][j].id);
			if (tmp_dis >= 0)
			{
				out << tmp_dis - line_component[i][j-1].len + _K  << " " ;
			}else
			{
				out << _position_tbl[line_component[i][j].id] - _position_tbl[line_component[i][j-1].id] - _length_tbl[line_component[i][j-1].id] + _K << " " ;
			}
		}
		out << endl;
	}
	out.close();
}

int cmp( Edge_Seq_Element es1, Edge_Seq_Element es2)
{
	return es1.pos < es2.pos;
}

void UniqEdgeGraph::initialize_component(string cmp_name) {
	cout << "initialize scaffolds" << endl;
	
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
    
    std::cout << boost::format("#### postion=[%d] componets=[%d] flag[%d],_nodelist[%d]") % _position_tbl.size() % _component.size() % flag.size() % _nodelist.size() << std::endl;
	
	for (int i = 0; i < _component.size(); i++)
	{
		sort(_component[i].begin(), _component[i].end(), cmp);
	}

	cout << "\tscaffolds number = "<< _component.size() << endl;

	ofstream out(cmp_name.c_str());
	
	for(int i = 0; i < _component.size(); i ++)
	{
		out << ">component " << i << endl;
		
		for (int j = 0; j < _component[i].size(); j ++)
		{
			out << _component[i][j].id << " ";
		}	
		out << endl;
		for (int j = 1; j < _component[i].size(); j ++)
		{
			out << (_component[i][j].pos - _component[i][j - 1].pos - _length_tbl[_component[i][j - 1].id] + _K ) << " ";
		}
		out << endl;
	}
	out.close();

    int conflict_scaf = 0;
    bool find_cft = false;

	for (int index = 0; index < _component.size(); index ++)
	{

        find_cft = false;

		for (int i = 0; i < _component[index].size() - 1; i ++)
		{
			for(int j = i + 1; j < _component[index].size(); j ++)
			{
				if (_component[index][i].pos + _component[index][i].len > _component[index][j].pos&& !hasEdge(_component[index][i].id, _component[index][j].id) &&
				!hasEdge(_component[index][j].id, _component[index][i].id))
				{
					if (_component[index][i].pos + _component[index][i].len < _component[index][j].pos + _component[index][j].len )
					{
						if (_max_overlap <= _component[index][i].pos + _component[index][i].len 
								- _component[index][j].pos)
						{
							overlap_pair.push_back(make_pair(_component[index][i].id, _component[index][j].id));
							overlap_com_id.push_back(index);
                            find_cft = true;
						}

					}else
					{
						if (_max_overlap <= _component[index][j].len)
						{
							overlap_pair.push_back(make_pair(_component[index][i].id, _component[index][j].id));
							overlap_com_id.push_back(index);
                            find_cft = true;

						}
					}
				}else
				{
					break;
				}
			}
		}

        if (find_cft == true)
        {
            conflict_scaf ++;
        }

	}
    cout << "\tscaffolds containing conflicts = " << conflict_scaf << endl;
	cout << "\tconflict overlap contigs pair number = " << overlap_pair.size() << endl;

}

int cmp_link( Ed e1, Ed e2)
{
	return e1.score < e2.score;
}

int avg_score(const vector<Ed> &tmp)
{
	if (tmp.size() == 0)
		return 0;
	int sum = 0;
	for (int i = 0; i < tmp.size(); i++ )
	{
		sum += tmp[i].score;
	}
	return sum / tmp.size();
}

Ed UniqEdgeGraph::get_min_ed(int index) {
	vector<Ed> tmp_ed_set;
	Ed tmp_ed;
	 
    // Children
    {
        NodeList::const_iterator i = _nodelist.find(index);
        if (i != _nodelist.end()) {
            for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
                tmp_ed.from = index;
                tmp_ed.to = j->first;
                tmp_ed.score = j->second.score;
                tmp_ed_set.push_back(tmp_ed);
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
                tmp_ed.from = *j;
                tmp_ed.to = index;
                Children::const_iterator l = k->second.children.find(i->first);
                tmp_ed.score = l->second.score;
                tmp_ed_set.push_back(tmp_ed);

            }
        }
    }

	sort(tmp_ed_set.begin(), tmp_ed_set.end(), cmp_link);
	
	if( tmp_ed_set.size() < 2)
	{
		tmp_ed.from = -1;
		tmp_ed.to = -1;
		return tmp_ed;
	}
	tmp_ed = tmp_ed_set[0];
	for (int i =1; i < tmp_ed_set.size(); i++)
	{
		tmp_ed_set[i-1] = tmp_ed_set[i];
	}
	tmp_ed_set.resize(tmp_ed_set.size() - 1);

	int avg = avg_score(tmp_ed_set);
	if (tmp_ed.score <  avg)
	{
		return tmp_ed;
	}else
	{
		tmp_ed.from = -1;
		tmp_ed.to = -1;
		return tmp_ed;
	}
}

void UniqEdgeGraph::remove_arc_con_edge_from_overlap_pair()
{
	int temp1, temp2;
	Ed tmp_ed;

	//ofstream fout("contig_removed");
	
    std::string file = boost::str(boost::format("chimeric_link_repeat_removed_log_%ld") % _iteration);
	ofstream fout(file.c_str());

	for (int i = 0; i < overlap_pair.size(); i++)
	{
		//while (temp1 >= 0 && temp1 != overlap_pair[i].first && temp1 != overlap_pair[i].second && !hasEdge(temp1, overlap_pair[i].first) && !hasEdge(temp1, overlap_pair[i].second))
		if ((hasEdge(overlap_pair[i].first, overlap_pair[i].second)) 
				|| 
				(hasEdge(overlap_pair[i].second, overlap_pair[i].first)) )
		{
			continue;
		}

		temp1 = getAncestor(overlap_pair[i].first, overlap_pair[i].second);	
		while (temp1 >= 0)
		{
			tmp_ed = get_min_ed(temp1);
			if (tmp_ed.from >= 0)
			{
				//fout << "del " << tmp_ed.from << "\t" << tmp_ed.to << endl;
				removeEdge(tmp_ed.from, tmp_ed.to);
				fout <<  "backward chimeric link " << tmp_ed.from << "," << tmp_ed.to << endl; 
			}else
			{
				break;
			}
			temp1 = getAncestor(overlap_pair[i].first, overlap_pair[i].second);
		}
		
		while (temp1 >= 0) 
		{
			fout << temp1 << "\tancestor of\t(" << overlap_pair[i].first << "," << overlap_pair[i].second << ")" << endl; 
			removeNode(temp1);

			temp1 = getAncestor(overlap_pair[i].first, overlap_pair[i].second);	
		}
		//while (temp2 >= 0 && temp2 != overlap_pair[i].first && temp2 != overlap_pair[i].second && !hasEdge( overlap_pair[i].first, temp2) && !hasEdge( overlap_pair[i].second, temp2)) 
		
		temp2 = getDescendant(overlap_pair[i].first, overlap_pair[i].second);
		while (temp2 >= 0)
		{
			tmp_ed = get_min_ed(temp2);
			if (tmp_ed.from >= 0)
			{
				//cout << "del " << tmp_ed.from << "\t" << tmp_ed.to << endl;
				removeEdge(tmp_ed.from, tmp_ed.to);
				fout <<  "forward chimeric link " << tmp_ed.from << "," << tmp_ed.to << endl; 
			}else
			{
				break;
			}
			temp2 = getDescendant(overlap_pair[i].first, overlap_pair[i].second);
		}

		while (temp2 >= 0) 
		{
			fout << temp2 << "\tdescendant of\t(" << overlap_pair[i].first << "," << overlap_pair[i].second << ")" << endl; 
			removeNode(temp2);

			temp2 = getDescendant(overlap_pair[i].first, overlap_pair[i].second);
		}
		
	}
	fout.close();
//	cout << "overlap pair num = " << overlap_pair.size() << endl;
	overlap_pair.clear();
	overlap_com_id.clear();
//	cout << "overlap pair num = " << overlap_pair.size() << endl;

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

	int depth = 2,count = 0, width = 0;
    std::deque<int> Q;
	Q.push_back(x);
	count = 1;

    // BFS
    // node=>x
    {
        int depth = 2, width = 0, count = 1;
        while (!Q.empty()) {
            size_t node = Q.front();
            Q.pop_front();
            flag[node] = 2;
        NodeList::const_iterator k = _nodelist.find(node);
        if (k != _nodelist.end()) {
            for (Children::const_iterator i = k->second.children.begin(); i != k->second.children.end(); ++i) {
                if(flag[i->first] == 0)
                {
                    flag[i->first] = 1;
                    Q.push_back(i->first);
                    width++;
                }
            }
        }
            if((--count) == 0){
                count = width;
                depth++;
            }
        }
    }

	Q.push_back(y);
	unsigned int id = 0;
	depth = -1;
	count = 1;
	width = 0;
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
			if (flag[i->first] == 0)
			{
				flag[i->first] = 1;
				Q.push_back(i->first);
				width++;
			}else if (flag[i->first] >= 2)
			{
				if(id == 0 || flag[i->first] > min){
					id = i->first;
					min = flag[id];
				}
			}
		}
    }
		if((--count) == 0){
			if(id != 0)
				return id;
			count = width;
			depth--;
		}
	}

	return -1;
}

void UniqEdgeGraph::output_graph(const std::string& filename) const {
	//cout << "begin output edge graph" << endl;
    std::string file = boost::str(boost::format("%s_%ld") % filename % _iteration);
	ofstream out(file.c_str());

    for (NodeList::const_iterator i = _nodelist.begin(); i != _nodelist.end(); ++i) {
        for (Children::const_iterator j = i->second.children.begin(); j != i->second.children.end(); ++j) {
				out << i->first << "|" << _length_tbl[i->first] << "|" << _position_tbl[i->first] << "\t" << j->first << "|" << _length_tbl[j->first] << "|" << _position_tbl[j->first] << "\t" << j->second.distance << "|" << j->second.count << "|" << j->second.score << endl;
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
