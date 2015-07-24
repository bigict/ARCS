/***************************************************************************
 * Description:
 *          
 * Author  : wangbing
 * Language: C++
 * Date    : 2015-01-29
 ***************************************************************************/
#include "Kmer_Hash_Md5.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <new>

using namespace std;

int Kmer_Hash_Md5::bufferid = 0;
int Kmer_Hash_Md5::nodeid = 0;
KNode* Kmer_Hash_Md5::ret = NULL;


Kmer_Hash_Md5::Kmer_Hash_Md5(unsigned int half_first_alloc):hash_table(NULL), last_alloc_num(half_first_alloc), 
											last_alloc_used(half_first_alloc), prime_pos(16), kmer_count(0),
											load_factor(0.9), hash_table_used(0ul){
	
	hash_table = new Kmer_Node_Ptr[ _prime_list[prime_pos] ];
	fill(hash_table, hash_table+_prime_list[prime_pos], static_cast<Kmer_Node_Ptr>(NULL));

	cout << "[INFO] construct Kmer_Hash_Md5" << endl;
}

unsigned long Kmer_Hash_Md5::size(){
	return kmer_count;
}

unsigned long Kmer_Hash_Md5::hash_table_length(){
	return _prime_list[prime_pos];
}

void Kmer_Hash_Md5::clear(){
	//cout << "[INFO] clear length=" << _prime_list[prime_pos] << " kmer_num=" << kmer_count << endl;
	if ( hash_table == NULL){
		return ;
	}
	
	for(int i=0; i<buffer_node.size(); i++){
		free( buffer_node[i] );
	}
	
	buffer_node.clear();
	last_alloc_num = 512;
	last_alloc_used = 512;
	delete[] hash_table;
	
	prime_pos = kmer_count = hash_table_used= 0;
	hash_table = NULL;
	cout << "[INFO] Kmer Hash cleared." << endl;
	return ;
}

Kmer_Hash_Md5::~Kmer_Hash_Md5(){
	//cout << "[INFO] distruct" << endl;
	if ( hash_table == NULL){
		return ;
	}
	for(int i=0; i<buffer_node.size(); i++){
		free( buffer_node[i] );
	}
	buffer_node.clear();
	last_alloc_num = 512;
	last_alloc_used = 512;
	delete[] hash_table;
	prime_pos = kmer_count = hash_table_used= 0;
	hash_table = NULL;
	return ;
}

Kmer_Hash_Md5::Kmer_Hash_Md5(const Kmer_Hash_Md5 &kh){
	cout << "can not be called" << endl;
}

Kmer_Hash_Md5& Kmer_Hash_Md5::operator=(const Kmer_Hash_Md5 &kh){
	cout << "can not be called" << endl;
	return *this;
}

void Kmer_Hash_Md5::rehash(){
	if (prime_pos == _primes_num-1){
		cout << "[WARNING] reach rehash limit" << endl;
		return ;
	}
	int new_prime_pos = prime_pos + 1;
	cout << "[INFO] rehash length = " << _prime_list[new_prime_pos] << endl;
	Kmer_Node_Ptr * tem_table = new Kmer_Node_Ptr[ _prime_list[new_prime_pos] ];
	memset(tem_table, 0, _prime_list[new_prime_pos]*sizeof( Kmer_Node_Ptr ));
	if (tem_table == NULL){
		cout << "bad_alloc" << endl;
		exit(1);
	}

	hash_table_used = 0;
	for(unsigned long i=0; i<_prime_list[prime_pos]; i++){
		if (hash_table[i] == NULL){
			continue;
		}
		Kmer_Node_Ptr tem = hash_table[i];
		Kmer_Node_Ptr next = tem->next;
		while(tem){
			unsigned long pos = tem->get_code() % (_prime_list[new_prime_pos]);
			if (tem_table[pos] == NULL){
				hash_table_used ++;
			}
			tem->next = tem_table[pos];
			tem_table[pos] = tem;
			tem = next;
			if (tem){
				next = tem->next;
			}
		}
	}

	prime_pos = new_prime_pos;
	delete[] hash_table;
	hash_table = tem_table;
	cout  << "[INFO] rehash Done " << endl;
	return ;
}

Kmer_Node_Ptr Kmer_Hash_Md5::has_kmer(const Kmer &kmer){
//	cout << "[INFO] has_kmer" << endl;
	unsigned long pos = kmer.get_code() % _prime_list[prime_pos];
	Kmer_Node_Ptr tem = hash_table[pos];
	while(tem){
		if (tem->kmer == kmer){
			return tem;
		}
		tem = tem->next;
	}
	return NULL;
}

static unsigned long cou = 1;
static unsigned long kmer_hashed = 0;
/*
void Kmer_Hash_Md5::push_a_kmer(const Kmer &kmer){
//	cout << "[INFO] push_a_kmer length=" << _prime_list[prime_pos] << " kmer_num=" << kmer_count << " last_alloc=" << last_alloc_num << " last_alloc_used=" << last_alloc_used << endl;
	kmer_hashed ++;
	if (kmer_hashed % cou == 0){
		cout << "[INFO] " << kmer_hashed << " kmers have already been hashed..." << endl;
		cou *= 10;
	}
	if ( (double)(hash_table_used + 1) >= _prime_list[ prime_pos ]*load_factor ){
		rehash();
	}
	if (last_alloc_num == last_alloc_used){
		Kmer_Node_Ptr buffer = (Kmer_Node_Ptr)malloc( last_alloc_num*2*sizeof(Kmer_Node));
		if (buffer == NULL){
			cout << "bad alloc" << endl;
			exit(1);
		}
		last_alloc_num *= 2;
		last_alloc_used = 0;
		buffer_node.push_back(buffer);
	}

	unsigned long pos = kmer.get_code() % _prime_list[prime_pos];
	Kmer_Node_Ptr tem = hash_table[pos];

	if(hash_table[pos] == NULL)
		hash_table_used ++;
	else
		do{
			if (tem->kmer.first_code == kmer.first_code){
				tem->num ++;
				return ;
			}
			tem = tem->next;
		}while(tem);

	Kmer_Node_Ptr tem_node = new( buffer_node[ buffer_node.size()-1 ]+last_alloc_used) Kmer_Node(kmer);

	tem_node->next = hash_table[pos];
	hash_table[pos] = tem_node;
	last_alloc_used ++;
	kmer_count ++;
}
*/

void Kmer_Hash_Md5::push_a_kmer(unsigned long first, unsigned long second, unsigned long third){
	kmer_hashed ++;
	if (kmer_hashed % cou == 0){
		cout << "[INFO] " << kmer_hashed << " kmers have already been hashed..." << endl;
		cou *= 10;
	}

	unsigned long pos = (first | second | third)  % _prime_list[prime_pos];
	Kmer_Node_Ptr tem = hash_table[pos];

	if(hash_table[pos] == NULL)
		hash_table_used ++;
	else
		do{
			if (tem->kmer.first_code == first && tem->kmer.second_code == second && tem->kmer.third_code == third){
				tem->num ++;
				return ;
			}
			tem = tem->next;
		}while(tem);

	if ( (double)(hash_table_used + 1) >= _prime_list[ prime_pos ]*load_factor ){
		rehash();
	}
	if (last_alloc_num == last_alloc_used){
		last_alloc_num *= 1.5;
		last_alloc_num = last_alloc_num > ONEG ? ONEG : last_alloc_num;
		Kmer_Node_Ptr buffer = (Kmer_Node_Ptr)malloc( last_alloc_num*sizeof(Kmer_Node));
		// cout << "[TEST] void Kmer_Hash_Md5::push_a_kmer(unsigned long) : buffer malloc. (" << last_alloc_num*2*sizeof(Kmer_Node)/(1024*1024) << "M) \n" ; 
		if (buffer == NULL){
			cout << "bad alloc" << endl;
			exit(1);
		}
		last_alloc_used = 0;
		buffer_node.push_back(buffer);
		buffer_size.push_back(0);
	}


	Kmer_Node_Ptr tem_node = new( buffer_node[ buffer_node.size()-1 ]+last_alloc_used) Kmer_Node(Kmer(first, second, third));
	buffer_size[ buffer_node.size()-1 ]++;

	tem_node->next = hash_table[pos];
	hash_table[pos] = tem_node;
	last_alloc_used ++;
	kmer_count ++;
}

unsigned int Kmer_Hash_Md5::get_kmer_num(const Kmer &kmer){
	Kmer_Node_Ptr tem = has_kmer(kmer);
	if (tem){
		return tem->num;
	}
	return 0;
}


Kmer_Node_Ptr& Kmer_Hash_Md5::operator[](unsigned long i){
	return hash_table[i];
}

void Kmer_Hash_Md5::filter(){
	cout << "[INFO] Filtering out Kmer with coverage < 2 before filet kmer num=" << kmer_count << endl;
	for(unsigned long i=0; i<_prime_list[prime_pos]; i++){
		if (hash_table[i] == NULL){
			continue;
		}
		Kmer_Node_Ptr tem = hash_table[i];
		while(tem){
			if (tem->num < 2){
				hash_table[i] = tem->next;
				tem = tem->next;
				kmer_count --;
			}else{
				break;
			}
		}
		if (hash_table[i] == NULL){
			continue;
		}
		Kmer_Node_Ptr pre = hash_table[i];
		tem = pre->next;
		while(tem){
			if (tem->num < 2){
				pre->next = tem->next;
				tem = tem->next;
				kmer_count --;
			}else{
				pre = tem;
				tem = tem->next;
			}
		}
	}
	cout << "[INFO] after filter kmer num =" << kmer_count << endl;
}

void Kmer_Hash_Md5::show_all_kmers(){
	cout << "[INFO] hash table size=" << _prime_list[prime_pos] << " kmer size=" << kmer_count << endl;
	cout << "[TEST] show all kmers" << endl;
	int list_size_morethan1 = 0;
	int empty_node = 0;
	int count[11];
	memset(count, 0, sizeof(count));
	for(unsigned long i=0; i<_prime_list[prime_pos]; i++){
		if (hash_table[i] ==NULL){
			empty_node ++;
			continue;
		}
		int list_size = 0;
		Kmer_Node_Ptr tem = hash_table[i];
		while(tem){
			list_size ++;
			tem = tem->next;
		}
		if (list_size > 1){
			list_size_morethan1 ++;
		}
		if (list_size >= 10){
			count[10] ++;
		}else{
			count[ list_size ]++;
		}
	}
	cout << "[INFO] " << list_size_morethan1 << " hash node have more than 1 kmer, " << empty_node << " empty node, " << _prime_list[prime_pos] << "total node" << endl;
	for(int i=1; i<11; i++){
		cout << count[i] << " ";
	}
	cout << endl;
}

void Kmer_Hash_Md5::cal_lamada(double &avg,double &var)
{
	//cout << "Calculating k-mer coverage distribution" << endl;
	int count = 0;
	double sum = 0;
	vector<int> cov_vec;
	unsigned long i;
	unsigned long hash_length = hash_table_length();
	for(i=0; i<hash_length/3; i++){
		if (hash_table[i] ==NULL){
			continue;
		}
		Kmer_Node_Ptr tem = hash_table[i];
		while(tem){
			if(tem->num > 3){
				count ++;
				sum += tem->num;
				cov_vec.push_back(tem->num);
			}
			tem = tem->next;
		}
	}
	avg = sum / count;
	cout << "\t[test]Initial k-mer coverage = " << avg << endl;

	if(avg > 100){
		sum = 0;
		count = 0;
		cov_vec.clear();
		for(i=0; i<hash_length; i++){
			if (hash_table[i] ==NULL){
				continue;
			}
			Kmer_Node_Ptr tem = hash_table[i];
			while(tem){
				if(tem->num >= avg/2.0 && tem->num <= avg*2){
					count ++;
					sum += tem->num;
					cov_vec.push_back(tem->num);
				}
				tem = tem->next;
			}
		}
	}
	avg = sum / count;
	cout << "\t[test]Ultimate k-mer coverage = " << avg << endl;
	sum = 0;
	for (int i = 0; i < cov_vec.size(); i++)
	{
		sum += pow(cov_vec[i] - avg, 2);
	}
	sum = sum/(cov_vec.size() - 1);
	var = sqrt(sum);

	cov_vec.clear();
}


void Kmer_Hash_Md5::filter2()	// change buffer size
{
	cout << "[INFO] Filtering out Kmer with coverage < 2 before filet kmer num=" << kmer_count << endl;
	int n = buffer_size.size();
	KNode *cnt = reinterpret_cast<KNode*>(buffer_node[0]);
	int index = 0;
	int last_buffer = 0;
	Kmer_Node_Ptr p = NULL;
	for(int i=0; i<n; i++){
		int m = buffer_size[i];
		int lastm = buffer_size[last_buffer]*1.2;
		for(int j=0; j<m; j++){
			p = buffer_node[i] + j;
			if(p->num > 1){
				if(index >=lastm){	// 40B -> 32B
					buffer_size[last_buffer] = lastm;
					last_buffer++;
					cnt = reinterpret_cast<KNode*>(buffer_node[last_buffer]);
					index = 0;
				}
				new (cnt + index++) KNode(*p);
			}else{
				kmer_count--;
			}
		}
	}

	buffer_size[last_buffer] = index;
	for(int i=n-1; i>=last_buffer+1; i--){
		cout << "[TEST] free buffer " << i << endl;
		Kmer_Node_Ptr p = buffer_node[i];
		buffer_node.pop_back();
		buffer_size.pop_back();
		free(p);
	}
	cout << "[INFO] after filter kmer num =" << kmer_count << endl;
}


KNode* Kmer_Hash_Md5::next_node() // traverse all node
{	
	if(ret == NULL)	ret = reinterpret_cast<KNode*>(buffer_node[0]);
	if(nodeid >= buffer_size[bufferid]){
		bufferid++;
		ret = reinterpret_cast<KNode*>(buffer_node[bufferid]);
		nodeid = 0;
	}
	if(bufferid >= buffer_size.size()){
		bufferid = nodeid = 0;
		ret = NULL;
		return NULL;
	}
	return ret + nodeid++;
}
