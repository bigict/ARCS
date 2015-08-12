#include "uniq_edge_graph.h"

#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <assert.h>
#include <stack>
#include <queue>
#include <string.h>
//#include <stdio.h>

#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("remove_repeats.uniq_edge_graph"));

void UniqEdgeGraph::add_a_len(int a_len)
{
	len.push_back(a_len);
}

void UniqEdgeGraph::add_a_pos(int pos_tmp)
{
	pos.push_back(pos_tmp);
}

void UniqEdgeGraph::add_a_dis(int i, int j, int c_dis, int cc, int sc)
{
	if(i >= con.size() | j >= con.size())
	{
		cout << "i or j out of range" << endl;
		exit(0);
	}

	list<Dis_Node>::iterator it;
	bool find = false;
	for(it = con[i].begin(); it != con[i].end(); it++)
	{
		if(it->id == j)
		{
			find = true;
			it->dis = c_dis;
			it->c = cc;
			it->score = sc;
			break;
		}	
	}
	if(find == false)
	{
		Dis_Node nd;
		nd.id = j;
		nd.dis = c_dis;
		nd.c = cc;
		nd.score = sc;
		con[i].push_back(nd);
	}
}

bool UniqEdgeGraph::is_a_dis(int i, int j)
{
	if(i >= con.size() | j >= con.size())
	{
		return false;
	}

	list<Dis_Node>::iterator it;

	for(it = con[i].begin(); it != con[i].end(); it++)
	{
		if(it->id == j)
		{
			return true;
		}	
	}
	return false;
}

void UniqEdgeGraph::add_a_arc(int i, int j, int c_dis, int cc, int sc)
{
	if(i >= arc.size() | j >= arc.size())
	{
		cout << "i or j out of range" << endl;
		exit(0);
	}
	list<Dis_Node>::iterator it;
	bool find = false;
	for(it = arc[i].begin(); it != arc[i].end(); it++)
	{
		if(it->id == j)
		{
			find = true;
			it->dis = c_dis;
			it->c = cc;
			it->score = sc;
			break;
		}	
	}

	if(find == false)
	{
		Dis_Node nd;
		nd.id = j;
		nd.dis = c_dis;
		nd.c = cc;
		nd.score = sc;
		arc[i].push_back(nd);
	}
}

void UniqEdgeGraph::add_a_rev_arc(int i, int j, int c_dis, int cc, int sc)
{
	if(i >= rev_arc.size() | j >= rev_arc.size())
	{
		return;
	}
	list<Dis_Node>::iterator it;
	bool find = false;
	for(it = rev_arc[i].begin(); it != rev_arc[i].end(); it++)
	{
		if(it->id == j)
		{
			find = true;
			it->dis = c_dis;
			it->c = cc;
			it->score = sc;
			break;
		}	
	}
	if(find == false)
	{
		Dis_Node nd;
		nd.id = j;
		nd.dis = c_dis;
		nd.c = cc;
		nd.score = sc;
		rev_arc[i].push_back(nd);
	}
}

int UniqEdgeGraph::get_dis(int i, int j)
{
	if(i >= arc.size() | j >= arc.size())
	{
		return 0;
	}
	list<Dis_Node>::iterator it;
	for(it = arc[i].begin(); it != arc[i].end(); it++)
	{
		if(it->id == j)
		{
			return it->dis;
		}
	}
	return -1;
}


void UniqEdgeGraph::del_a_dis(int i, int j)
{
	if(i >= con.size() | j >= con.size())
	{
		cout << "out of range" << endl;
	}
	list<Dis_Node>::iterator it = con[i].begin();
	while( it != con[i].end())
	{
		if(it->id == j)
		{
			it = con[i].erase(it);
		}else
		{
			it ++;
		}
	}
}

void UniqEdgeGraph::del_a_arc(int i, int j)
{
	if(i >= arc.size() | j >= arc.size())
	{
		cout << "out of range" << endl;
	}
	list<Dis_Node>::iterator it = arc[i].begin();
	while( it != arc[i].end())
	{
		if(it->id == j)
		{
			it = arc[i].erase(it);
		}else
		{
			it ++;
		}
	}
}

void UniqEdgeGraph::del_a_rev_arc(int i, int j)
{
	if(i >= rev_arc.size() | j >= rev_arc.size())
	{
		cout << "out of range" << endl;
	}
	list<Dis_Node>::iterator it = rev_arc[i].begin();
	while( it != rev_arc[i].end())
	{
		if(it->id == j)
		{
			it = rev_arc[i].erase(it);
		}else
		{
			it ++;
		}
	}
}

void UniqEdgeGraph::input_inner_component() {
    
    std::string file = boost::str(boost::format("component_%ld") % _iteration);
	ifstream in(file.c_str());
	
	if (!in) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
		exit(0);
	}
	
	string line;
	unsigned int count = 0;
	
	while(!in.eof())
	{
		getline(in, line);
		count ++;
	}
	inner_component.resize(count / 3);
	gap.resize(count / 3);
	in.close();
	in.open(file.c_str());
	
	char *word;
	count = 0;
	
	while(!in.eof())
	{
		getline(in, line);
		getline(in, line);
		word = strtok(const_cast<char*>(line.c_str()), " \t,");
		while(word)
		{
			//cout << word << endl;
			inner_component[count].push_back(atoi(word));
			word = strtok(NULL, " \t,");
		}
		getline(in, line);
		word = strtok(const_cast<char*>(line.c_str()), " \t,");
		while(word)
		{
			//cout << word << endl;
			gap[count].push_back(atoi(word));
			word = strtok(NULL, " \t,");
		}
		count ++;
	}
/*
	ofstream out("inner_component");
	for (int i = 0; i < inner_component.size(); i ++)
	{
		out << ">component " << i << endl;
		for (int j = 0; j < inner_component[i].size(); j ++)
		{
			out << inner_component[i][j] << " ";
		}
		out << endl;
		
		for(int j = 0; j < gap[i].size(); j ++)
		{
			out << gap[i][j] << " ";
		}
		out << endl;
	}
	out.close();
*/
}

void UniqEdgeGraph::input_edge_pos() {

    std::string file = boost::str(boost::format("edge_cluster_pos_%ld") % _iteration);
	ifstream fin(file.c_str());
	if(!fin) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
		exit(0);
	}
	
	int pos_tmp = 0;
	string line;
	int eql_pos = 0;

	while(getline(fin, line))
	{
		pos_tmp = atoi(line.c_str());
		add_a_pos(pos_tmp);
	}
}

void UniqEdgeGraph::input_edge_len() {

    std::string file = boost::str(boost::format("edge_cluster_len_%ld") % _iteration);
	ifstream fin(file.c_str());
	if(!fin) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
		exit(0);
	}
	int len_tmp = 0;
	string line;
	while(getline(fin, line))
	{	
		len_tmp = atoi(line.c_str());
		add_a_len(len_tmp);
	}
}

void UniqEdgeGraph::input_edge_link(string name) {
    std::string file = boost::str(boost::format("%s_%ld") % name % _iteration);
	ifstream fin(file.c_str());
	
	if (!fin) {
        LOG4CXX_ERROR(logger, boost::format("%s open failed") % file);
		exit(0);
	}
 	int from, to, dis, c;
	double score;
	string line;
	while (getline(fin, line))
	{
		sscanf(line.c_str(), "%d\t%d\t%d\t%d\t%lf", &from, &to, &dis, &c, &score);
		add_a_dis(from, to, dis, c, score);
		add_a_arc(from, to, dis, c, score);
		add_a_rev_arc(to, from, -dis, c, score);
	}
	
	fin.close();
}

void UniqEdgeGraph::linearize()
{
	
	//output_graph("tmp_graph_before_repeats_removing.data");
	transform_bidirection();
	initialize_component("new_component");
	transform_single_direction();
	remove_arc_con_edge_from_overlap_pair();	
	//cout << "this time remove overlap end" << endl;
	
	transform_bidirection();
	initialize_component("new_component");
	transform_single_direction();

	output_graph("contig_arc_graph_after_repeats_removing");	
	tran_to_line();
}

void UniqEdgeGraph::tran_to_line()
{
	//cout << "begin transform to line" << endl;

    vector<int> flag;
	flag.resize(pos.size());
	line_component.resize(component.size());
	
	//cout << "component size = " << component.size() << endl;

	for(int i = 0; i < flag.size(); i++)
	{
		flag[i] = 0;
	}
	
	//ofstream fout("no_end");

	stack<int> st;
	for(int i = 0; i < component.size(); i ++)
	{
		int size = component[i].size();
		int first = component[i][0].id;
		int end = component[i][size - 1].id;
		
		if(first == end)
		{
			line_component[i].resize(1);
			line_component[i][0].id = first;
			line_component[i][0].len = len[first];
			line_component[i][0].pos = pos[first];
			continue;
		}

		while (!st.empty())
		{
			st.pop();
		}

		// for(int j=0; j<size; j++){
		// 	st.push(component[i][j].id);
		// }

		flag[first] = 1;
		int temp;
		bool find_next = false;
		bool find_tar = false;

		while(!st.empty())
		{
			temp = st.top();
			find_next = false;
			for(list<Dis_Node>::iterator it = con[temp].begin(); it != con[temp].end(); it ++)
			{
				if(flag[it->id] == 0 && it->id == end)
				{
					st.push(it->id);
					flag[it->id] = 1;
					find_next = true;
					find_tar = true;
					break;
				}else if (flag[it->id] == 0)
				{
					st.push(it->id);
					flag[it->id] = 1;
					find_next = true;
					break;
				}
			}

			if(find_tar == true)
			{
				break;
			}

			if(find_next == true)
			{
				continue;
			}else
			{
				flag[temp] = 2;
				st.pop();
			}
		}
		
		if (find_tar == false)
		{
	//		fout << first << endl;
			line_component[i] = component[i];
		}else
		{
			line_component[i].resize(st.size());
			for (int j = st.size() - 1; j >= 0; j--)
			{
				Edge_Seq_Element e_s_e;
				e_s_e.id = st.top();
				e_s_e.len = len[e_s_e.id];
				e_s_e.pos = pos[e_s_e.id];
				st.pop();
				line_component[i][j] = e_s_e;
			}
		}
	}
	
	for(int i = 0; i < flag.size(); i++)
	{
		flag[i] = 0;
	}

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
			del_a_edge(i);
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
			for (int k = 0; k < inner_component[line_component[i][j].id].size(); k ++)
			{
				out << inner_component[line_component[i][j].id][k] << " ";
			}		
		}
		out << endl;
		for (int k = 0; k < gap[line_component[i][0].id].size(); k ++)
		{
			out << gap[line_component[i][0].id][k] << " ";
		}
		for (int j = 1; j < line_component[i].size(); j ++)
		{	
			//out << pos[line_component[i][j].id] - pos[line_component[i][j-1].id] - len[line_component[i][j-1].id] + K << " " ;
			int tmp_dis = get_dis(line_component[i][j-1].id, line_component[i][j].id);
			if (tmp_dis >= 0)
			{
				out << tmp_dis - line_component[i][j-1].len + _K  << " " ;
			}else
			{
				out << pos[line_component[i][j].id] - pos[line_component[i][j-1].id] - len[line_component[i][j-1].id] + _K << " " ;
			}

			for (int k = 0; k < gap[line_component[i][j].id].size(); k ++)
			{
				out << gap[line_component[i][j].id][k] << " ";
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
			out << line_component[i][j].id << "|" << len[line_component[i][j].id] << "|" << pos[line_component[i][j].id] << "\t";
		}
		out << endl;
		for (int j = 1; j < line_component[i].size(); j ++)
		{	
			int tmp_dis = get_dis(line_component[i][j-1].id, line_component[i][j].id);
			if (tmp_dis >= 0)
			{
				out << tmp_dis - line_component[i][j-1].len + _K  << " " ;
			}else
			{
				out << pos[line_component[i][j].id] - pos[line_component[i][j-1].id] - len[line_component[i][j-1].id] + _K << " " ;
			}
		}
		out << endl;
	}
	out.close();
}

void UniqEdgeGraph::transform_single_direction()
{
	list<Dis_Node>::iterator it;
	for (int i = 0; i < con.size(); i++)
	{
		it = con[i].begin();
		while(it != con[i].end())
		{
			if(it->dis < 0)
			{
				it = con[i].erase(it);
			}else
			{
				it ++;
			}
		}
	}
}

int cmp( Edge_Seq_Element es1, Edge_Seq_Element es2)
{
	return es1.pos < es2.pos;
}

void UniqEdgeGraph::initialize_component(string cmp_name)
{
	cout << "initialize scaffolds" << endl;
	
	vector<int> flag;
	component.clear();

    // BFS
	flag.resize(con.size());
	for(int i = 0; i < flag.size(); i ++)
	{
		flag[i] = 0;
	}
	vector<Edge_Seq_Element> e_s;
	Edge_Seq_Element e_s_e;
	list<int> qu;
	int tmp;
	list<Dis_Node>::iterator it;
	for(int i = 0; i < flag.size(); i++)
	{
		if(flag[i] == 0)
		{
			e_s.clear();
			qu.push_back(i);
			while(!qu.empty())
			{
				tmp = qu.front();
				qu.pop_front();
				e_s_e.id = tmp;
				e_s_e.pos = pos[tmp];
				e_s_e.len = len[tmp];
				e_s.push_back(e_s_e);	
				flag[tmp] = 2;
				for(it = con[tmp].begin(); it != con[tmp].end(); it++)
				{
					if(flag[it->id] == 0)
					{
						flag[it->id] = 1;
						qu.push_back(it->id);
					}
				}			
			}
			component.push_back(e_s);
		}
	}
	
	for (int i = 0; i < component.size(); i++)
	{
		sort(component[i].begin(), component[i].end(), cmp);
	}

	vector<vector<int> > pre;
	pre.resize(con.size());
	for(int i = 0; i < con.size(); i++)
	{
		for(list<Dis_Node>::iterator it = con[i].begin(); it != con[i].end(); it ++)
		{
			if(it->dis > 0)
			{
				pre[it->id].push_back(i);
			}
		}
	}

	cout << "\tscaffolds number = "<< component.size() << endl;

	ofstream out(cmp_name.c_str());
	
	for(int i = 0; i < component.size(); i ++)
	{
		out << ">component " << i << endl;
		
		for (int j = 0; j < component[i].size(); j ++)
		{
			out << component[i][j].id << " ";
		}	
		out << endl;
		for (int j = 1; j < component[i].size(); j ++)
		{
			out << (component[i][j].pos - component[i][j - 1].pos - len[component[i][j - 1].id] + _K ) << " ";
		}
		out << endl;
	}
	out.close();

    int conflict_scaf = 0;
    bool find_cft = false;

	for (int index = 0; index < component.size(); index ++)
	{

        find_cft = false;

		for (int i = 0; i < component[index].size() - 1; i ++)
		{
			for(int j = i + 1; j < component[index].size(); j ++)
			{
				if (component[index][i].pos + component[index][i].len > component[index][j].pos&& !is_a_dis(component[index][i].id, component[index][j].id) &&
				!is_a_dis(component[index][j].id, component[index][i].id))
				{
					if (component[index][i].pos + component[index][i].len < component[index][j].pos + component[index][j].len )
					{
						if (_max_overlap <= component[index][i].pos + component[index][i].len 
								- component[index][j].pos)
						{
							overlap_pair.push_back(make_pair(component[index][i].id, component[index][j].id));
							overlap_com_id.push_back(index);
                            find_cft = true;
						}

					}else
					{
						if (_max_overlap <= component[index][j].len)
						{
							overlap_pair.push_back(make_pair(component[index][i].id, component[index][j].id));
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

Ed UniqEdgeGraph::get_min_ed(int index)
{
	vector<Ed> tmp_ed_set;
	Ed tmp_ed;
	list<Dis_Node>::iterator it;
	 
	for(it = arc[index].begin(); it != arc[index].end(); it++)
	{
		tmp_ed.from = index;
		tmp_ed.to = it->id;
		tmp_ed.score = it->score;
		tmp_ed_set.push_back(tmp_ed);
	}
	
	for(it = rev_arc[index].begin(); it != rev_arc[index].end(); it++)
	{
		tmp_ed.from = it->id;
		tmp_ed.to = index;
		tmp_ed.score = it->score;
		tmp_ed_set.push_back(tmp_ed);
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
		//while (temp1 >= 0 && temp1 != overlap_pair[i].first && temp1 != overlap_pair[i].second && !is_a_dis(temp1, overlap_pair[i].first) && !is_a_dis(temp1, overlap_pair[i].second))
		if ((is_a_dis(overlap_pair[i].first, overlap_pair[i].second)) 
				|| 
				(is_a_dis(overlap_pair[i].second, overlap_pair[i].first)) )
		{
			continue;
		}

		temp1 = get_ancestor(overlap_pair[i].first, overlap_pair[i].second);	
		while (temp1 >= 0)
		{
			tmp_ed = get_min_ed(temp1);
			if (tmp_ed.from >= 0)
			{
				//fout << "del " << tmp_ed.from << "\t" << tmp_ed.to << endl;
				del_a_dis(tmp_ed.from, tmp_ed.to);
				del_a_arc(tmp_ed.from, tmp_ed.to);
				del_a_rev_arc(tmp_ed.to, tmp_ed.from);
				fout <<  "backward chimeric link " << tmp_ed.from << "," << tmp_ed.to << endl; 
			}else
			{
				break;
			}
			temp1 = get_ancestor(overlap_pair[i].first, overlap_pair[i].second);
		}
		
		while (temp1 >= 0) 
		{
			fout << temp1 << "\tancestor of\t(" << overlap_pair[i].first << "," << overlap_pair[i].second << ")" << endl; 
			del_a_edge(temp1);

			temp1 = get_ancestor(overlap_pair[i].first, overlap_pair[i].second);	
		}
		//while (temp2 >= 0 && temp2 != overlap_pair[i].first && temp2 != overlap_pair[i].second && !is_a_dis( overlap_pair[i].first, temp2) && !is_a_dis( overlap_pair[i].second, temp2)) 
		
		temp2 = get_descendant(overlap_pair[i].first, overlap_pair[i].second);
		while (temp2 >= 0)
		{
			tmp_ed = get_min_ed(temp2);
			if (tmp_ed.from >= 0)
			{
				//cout << "del " << tmp_ed.from << "\t" << tmp_ed.to << endl;
				del_a_dis(tmp_ed.from, tmp_ed.to);
				del_a_arc(tmp_ed.from, tmp_ed.to);
				del_a_rev_arc(tmp_ed.to, tmp_ed.from);
				fout <<  "forward chimeric link " << tmp_ed.from << "," << tmp_ed.to << endl; 
			}else
			{
				break;
			}
			temp2 = get_descendant(overlap_pair[i].first, overlap_pair[i].second);
		}

		while (temp2 >= 0) 
		{
			fout << temp2 << "\tdescendant of\t(" << overlap_pair[i].first << "," << overlap_pair[i].second << ")" << endl; 
			del_a_edge(temp2);

			temp2 = get_descendant(overlap_pair[i].first, overlap_pair[i].second);
		}
		
	}
	fout.close();
//	cout << "overlap pair num = " << overlap_pair.size() << endl;
	overlap_pair.clear();
	overlap_com_id.clear();
//	cout << "overlap pair num = " << overlap_pair.size() << endl;

}

int UniqEdgeGraph::get_ancestor(int i, int j)
{
    std::vector<int> flag(rev_arc.size(), 0);
	int depth = 2,count = 0,k = 0;
    std::queue<int> qu;
	qu.push(i);
	count = 1;

    // BFS
	int tmp;
	list<Dis_Node>::iterator it;
	while (!qu.empty())
	{
		tmp = qu.front();
		qu.pop();
		flag[tmp] = depth;
		for (it = rev_arc[tmp].begin(); it != rev_arc[tmp].end(); it++)
		{
			if(flag[it->id] == 0)
			{
				flag[it->id] = 1;
				qu.push(it->id);
				k++;
			}
		}
		if((--count) == 0){
			count = k;
			depth++;
		}
	}

	qu.push(j);
	unsigned int id = 0;
	depth = -1;
	count = 1;
	k = 0;
	int min = 0;
	while (!qu.empty())
	{
		tmp = qu.front();
		qu.pop();
		if(flag[tmp] >= 2)
		{
			return tmp;
		}
		flag[tmp] = depth;
		for(it = rev_arc[tmp].begin(); it != rev_arc[tmp].end(); it++)
		{
			if (flag[it->id] == 0)
			{
				flag[it->id] = 1;
				qu.push(it->id);
				k++;
			}else if (flag[it->id] >= 2)
			{
				if(id == 0 || flag[it->id] < min){
					id = it->id;
					min = flag[id];
				}
			}
		}
		if((--count) == 0){
			if(id != 0)
				return id;
			count = k;
			depth--;
		}
	}
	return -1;
}

int UniqEdgeGraph::get_descendant(int i, int j)
{
	vector<int> flag;
	int depth = 2,count = 0,k = 0;
	flag.resize(arc.size());
	for(int index = 0; index < flag.size(); index ++)
	{
		flag[index] = 0;
	}
	queue<int> qu;
	qu.push(i);
	count = 1;

	int tmp;
	list<Dis_Node>::iterator it;
	while (!qu.empty())
	{
		tmp = qu.front();
		qu.pop();
		flag[tmp] = 2;
		for(it = arc[tmp].begin(); it != arc[tmp].end(); it++)
		{
			if(flag[it->id] == 0)
			{
				flag[it->id] = 1;
				qu.push(it->id);
				k++;
			}
		}
		if((--count) == 0){
			count = k;
			depth++;
		}
	}

	qu.push(j);
	unsigned int id = 0;
	depth = -1;
	count = 1;
	k = 0;
	int min = 0;
	while (!qu.empty())
	{
		tmp = qu.front();
		qu.pop();
		if (flag[tmp] >= 2)
		{
			return tmp;
		}
		flag[tmp] = 3;
		for(it = arc[tmp].begin(); it != arc[tmp].end(); it++)
		{
			if (flag[it->id] == 0)
			{
				flag[it->id] = 1;
				qu.push(it->id);
				k++;
			}else if (flag[it->id] >= 2)
			{
				if(id == 0 || flag[it->id] > min){
					id = it->id;
					min = flag[id];
				}
			}
		}
		if((--count) == 0){
			if(id != 0)
				return id;
			count = k;
			depth--;
		}
	}
	return -1;
}

void UniqEdgeGraph::del_a_edge(int i)
{
	if(i >= con.size() || i < 0)
	{
		cout << "out of range" << endl;
		return; 
	}
	list<Dis_Node>::iterator it;
	list<Dis_Node>::iterator it_next;
	for(it = rev_arc[i].begin(); it != rev_arc[i].end(); it++)
	{
		for(it_next = con[it->id].begin(); it_next != con[it->id].end(); it_next++)
		{
			if(it_next->id == i)
			{
				con[it->id].erase(it_next);
				break;
			}	
		}
	}
	con[i].clear();
	
	for(it = rev_arc[i].begin(); it != rev_arc[i].end(); it++)
	{
		for(it_next = arc[it->id].begin(); it_next != arc[it->id].end(); it_next++)
		{
			if(it_next->id == i)
			{
				arc[it->id].erase(it_next);
				break;
			}	
		}
	}

	for(it = arc[i].begin(); it != arc[i].end(); it++)
	{
		for(it_next = rev_arc[it->id].begin(); it_next != rev_arc[it->id].end(); it_next++)
		{
			if(it_next->id == i)
			{
				rev_arc[it->id].erase(it_next);
				break;
			}	
		}
	}
	rev_arc[i].clear();
	arc[i].clear();
}

void UniqEdgeGraph::transform_bidirection()
{
//	cout << "uniq edge graph tran bidirect" << endl;
	list<Dis_Node>::iterator it;
	for(int i = 0; i < con.size(); i++)
	{
		for(it = con[i].begin(); it != con[i].end(); it++)
		{
			add_a_dis(it->id, i, -(it->dis), it->c, it->score);
		}
	}
//	cout << "tran end" << endl;
}


void UniqEdgeGraph::resize(int n)
{
	arc.resize(n);
	con.resize(n);
	//pos.resize(n);
	rev_arc.resize(n);
}

void UniqEdgeGraph::output_graph(string s)
{
	//cout << "begin output edge graph" << endl;
    std::string file = boost::str(boost::format("%s_%ld") % s % _iteration);
	ofstream out(file.c_str());

	list<Dis_Node>::iterator it;
	for(int i = 0; i < con.size(); i ++)
	{
		for(it = con[i].begin();it != con[i].end(); it++)
		{
			//if (it->dis > 0)
			{
				out << i << "|" << len[i] << "|" << pos[i] << "\t" << it->id << "|" << len[it->id] << "|" << pos[it->id] << "\t" << it->dis << "|" << it->c << "|" << it->score << endl;
			}
		}
	}
	out.close();
	//cout << "end out put graph" << endl;
}


