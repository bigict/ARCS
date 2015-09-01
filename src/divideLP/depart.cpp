//depart mixture LP to many one small LP
#include<iostream>
#include<sstream>
#include<cstdio>
#include<fstream>
#include<cmath>
#include<string>
#include<cstring>
#include<algorithm>
#include<iterator>
#include<map>
#include<set>
#include<queue>
#include<vector>
using namespace std;
#define M 999999

int edge=0; //count the number of edges
const int STD = 10;

set<int> z_mp;

vector<string> data_file;

vector<vector<pair<int,int> > > z_s;
vector<vector<int > > z_ss;
map<int, vector<pair<int,int> > > dis;
	// set<int> s[M];
set<int>::iterator is;
//get pair<x,d> related to x[i]
// struct node
// {
// 	vector<pair<int,int> >vxd;
// }vx[M];
map<int, vector<int> > z_vx;
int maxI = 0;

int x[M];
int result[M];

void read(string file_path)
{
	ifstream input;
	string data_in_file;
	input.open(file_path.c_str());
	while(!input.eof())
	{
		getline(input,data_in_file);
		data_file.push_back(data_in_file);
	}
	vector<string>::iterator i;
	string sr;
	
	for(i=data_file.begin();i!=data_file.end();i++)
	{
		sr=(*i);
		int p=sr.find(" : x_");
		if(p>0)
		{
			edge++;
			int p1=p+5;
			int p2=sr.find("- x_")+4;
			int p3=sr.find(" = ")+3;
			int tem1=0;
			while(sr[p1]>='0'&&sr[p1]<='9')
			{
				tem1=tem1*10+(sr[p1]-'0');
				p1++;
			}
			int tem2=0;
			while(sr[p2]>='0'&&sr[p2]<='9')
			{
				tem2=tem2*10+(sr[p2]-'0');
				p2++;
			}
			z_mp.insert(tem1);
			z_mp.insert(tem2);
			int tem3=0;
			int flag = 0;
			if(sr[p3] == '-')
				flag = 1;
			else
				tem3 = sr[p3] - '0';
			while(sr[++p3]>='0'&&sr[p3]<='9')
			{
				tem3=tem3*10+(sr[p3]-'0');
				// p3++;
			}
			if(flag)
				tem3 = -tem3;
			if(dis.find(tem2) == dis.end())
				dis[tem2] = vector<pair<int, int> >();
			dis[tem2].push_back(make_pair(tem1,tem3));
   			//mp.insert(pair<int,int>(tem1,1));
			// mp.insert(pair<int,int>(tem2,1));
			
			if( z_vx.find(tem1) == z_vx.end() ) {
				z_vx[tem1] = vector<int>();
			}
			z_vx[tem1].push_back(tem2);
			if( z_vx.find(tem2) == z_vx.end() ) {
				z_vx[tem2] = vector<int>();
			}
			z_vx[tem2].push_back(tem1);

			int tmp = tem1 > tem2 ? tem1 : tem2;
			maxI = maxI > tmp ? maxI : tmp;
			// vx[tem1].vxd.push_back(pair<int,int>(tem2,1));
			// vx[tem2].vxd.push_back(pair<int,int>(tem1,1));
		}	
	}
}

//init x[]
int ct=0; //count the num of LP to be solved
void init()
{
	queue<int> q;
	unsigned int cnt=z_mp.size();
	vector<bool> visited(maxI+1, false);
	// z_s.resize(maxI+1);
	// memset(visited,0,sizeof(visited));
	set<int>::iterator it = z_mp.begin();
	//int num = maxI;
	unsigned int count = 0;
	unsigned int maxNum = 0;
	cout << cnt << endl;
	while(count < cnt)
	{
		// vector<pair<int,int> > tmp;
		vector<int> tmp1;
		for(;it!=z_mp.end();it++)
		{	
			if(visited[*it]==0)
			{
				q.push(*it);
				visited[*it] = 1;  //define the first vertex's position is 1
				// z_s[ct].insert(*it);
				tmp1.push_back(*it);
				// cnt--;
				count ++;
				break;
			}
		}
		while(!q.empty())
		{
			unsigned int t=q.front();
			q.pop();
			for(unsigned int i=0;i<z_vx[t].size();i++)
			{
				int j = z_vx[t][i];
				if(visited[j]==0)
				{
					// tmp.push_back(make_pair(t, j));
					tmp1.push_back(j);
					q.push(j);
					visited[j]=1;
					count ++;
				}	
			}
		}
		// cnt -= tmp.size();
		// z_s.push_back(tmp);
		z_ss.push_back(tmp1);
		if(tmp1.size() > 1)
			ct++;
		maxNum = maxNum > tmp1.size()? maxNum : tmp1.size();
	}
	cout<<"[info] LP number is:"<< ct <<endl;
	cout<<"[info] Max nodes :"<< maxNum <<endl;
}


string intTostr(int k)
{
	string s;
	stringstream ss;
	ss<<k;
	s=ss.str();
	return s;
}
void output_lp(string out_path)
{
	cout << "output linear programing mod" << endl;
	//char buff[100];

	string ss;
	for(int k=0; k<ct; k++) {
		if(z_ss[k].size() <= 1)
			continue;
		ss = out_path + "/lp";
		ss.append(intTostr(k)).append(".math");
		ofstream fout(ss.c_str());

		for (unsigned int i = 0; i < z_ss[k].size(); i++)
		{
			fout << "var x_" << z_ss[k][i] << ";" << endl;
		}
		int x1,x2;
		for(unsigned int i = 0; i < z_ss[k].size(); i++)
		{
			if(dis.find(z_ss[k][i]) == dis.end())
				continue;
			for(unsigned int j=0; j < dis[z_ss[k][i]].size(); j++){
				x1 = z_ss[k][i];
				x2 = dis[z_ss[k][i]][j].first;			
				fout << "var e_" << x1 << "_" << x2  << ";" << endl; 
				fout << "var E_" << x1 << "_" << x2  << ";" << endl;
			}
		}
		fout << endl;
		fout << "minimize z:   ";

		int count = 0;
		for(unsigned int i = 0; i < z_ss[k].size(); i++)
		{
			if(dis.find(z_ss[k][i]) == dis.end())
				continue;
			for(unsigned int j=0; j < dis[z_ss[k][i]].size(); j++){
				x1 = z_ss[k][i];
				x2 = dis[z_ss[k][i]][j].first;			
			
				fout << " E_" << x1 << "_" << x2  << " + ";
			    count ++;
			    if (count % 10 == 0)
				    fout << endl;
			}
		}
		fout << "0;" << endl << endl;


		int index = 1;
		for(unsigned int i = 0; i < z_ss[k].size(); i++)
		{
			if(dis.find(z_ss[k][i]) == dis.end())
				continue;
			for(unsigned int j=0; j < dis[z_ss[k][i]].size(); j++){
				x1 = z_ss[k][i];
				x2 = dis[z_ss[k][i]][j].first;	
				
				fout << "s.t. con" << index ++ << " : x_" << x2 << " - x_" << x1 << " + e_" << x1 << "_" << x2 << " = " << dis[z_ss[k][i]][j].second << ";" << endl;
				fout << "s.t. con" << index ++ << " : E_" << x1 << "_" << x2 << " + e_" << x1 << "_" << x2 << " >= 0;" << endl;
				fout << "s.t. con" << index ++ << " : E_" << x1 << "_" << x2 << " - e_" << x1 << "_" << x2 << " >= 0;" << endl;
			}
		}

		fout << endl << "end;";

		fout.close();
	}
}


 
int main(int argv, char **args) {
	if (argv != 3 {
		cout << "Usage: " << args[0] << " bigLPFile outSmallLPPath" << endl;
		exit(1); 
	}
	//depart data
	read(string( args[1] ));
	//init data
	init();
	// pushToVector( string(args[1]) );
	// newFile( string(args[2]));
	output_lp( string(args[2]) );
	return 0;
}
