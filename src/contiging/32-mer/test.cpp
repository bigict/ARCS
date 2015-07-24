#include "Kmer.h"
#include "Kmer_Hash_Md5.h"
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;

int K;

unsigned long test_get_kmer_num(string file, Kmer_Hash_Md5 &kmer_hash){
	ifstream in(file.c_str());

	string line;
	unsigned long count = 0;
	while(getline(in, line)){
		getline(in, line);
		if (!in){
			break;
		}
		for(int i=0; i<line.length()-K+1; i++){
			string tem = line.substr(i, K);
			int num = kmer_hash.get_kmer_num( Kmer(tem) );
			count ++;
		}
		getline(in, line);
		getline(in, line);
		if (!in){
			break;
		}
	}
	in.close();
	return count;
}

int main (int argc, char ** argv){
	if (argc != 3){
		cout << "Usage: " << argv[0] << " K readfile" << endl;
		exit(1);
	}
	
	time_t start, end;
	start = time(NULL);
	Kmer_Hash_Md5 kmer_hash(2);
	K = atoi(argv[1]);
	ifstream in(argv[2]);
	string line;
	while(getline(in, line)){
		getline(in, line);
		if (!in){
			break;
		}
//		cout << "[" << line << "]" << endl;
		for(int i=0; i<line.length()-K+1; i++){
			string tem = line.substr(i, K);
//			cout << tem << endl;
			kmer_hash.push_a_kmer( Kmer(tem) );
		}
		getline(in, line);
		getline(in, line);
		if (!in){
			break;
		}
	}

	kmer_hash.show_all_kmers();
	kmer_hash.filter();
//	kmer_hash.show_all_kmers();
	end = time(NULL);
	cout << "kmer hash time=" << difftime(end, start) << " seconds" << endl;
	in.close();
	
	start = time(NULL);
	unsigned long num = test_get_kmer_num( argv[2] , kmer_hash);
	end = time(NULL);
	cout << "call get_kmer_num " << num << " times use time=" << difftime(end, start) << " seconds" << endl;

	kmer_hash.clear();
	return 0;
}
