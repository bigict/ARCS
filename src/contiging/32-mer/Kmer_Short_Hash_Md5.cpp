/***************************************************************************
 * Description:
 *          
 * Author  : wangbing
 * Language: C++
 * Date    : 2015-01-29
 ***************************************************************************/
#include "Kmer_Short_Hash_Md5.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <new>

using namespace std;



Kmer_Short_Hash_Md5::Kmer_Short_Hash_Md5(unsigned int half_first_alloc):hash_table(NULL), last_alloc_num(half_first_alloc), last_alloc_used(half_first_alloc), prime_pos(0), kmer_count(0),load_factor(0.88), hash_table_used(0){
	
	hash_table = new Kmer_Short_Node_Ptr[ _prime_list[prime_pos] ];
	memset(hash_table,0,sizeof(Kmer_Short_Node_Ptr)*_prime_list[prime_pos]);
	buffers = 0;
	cout << "[INFO] construct Kmer_Short_Hash_Md5" << endl;
}

unsigned long Kmer_Short_Hash_Md5::size(){
	return kmer_count;
}

unsigned long Kmer_Short_Hash_Md5::hash_table_length(){
	return _prime_list[prime_pos];
}

void Kmer_Short_Hash_Md5::clear(){
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
	cout << "[INFO] Short Kmer Hash cleared." << endl;
	return ;
}

Kmer_Short_Hash_Md5::~Kmer_Short_Hash_Md5(){
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

Kmer_Short_Hash_Md5::Kmer_Short_Hash_Md5(const Kmer_Short_Hash_Md5 &kh){
	cout << "can not be called" << endl;
}

Kmer_Short_Hash_Md5& Kmer_Short_Hash_Md5::operator=(const Kmer_Short_Hash_Md5 &kh){
	cout << "can not be called" << endl;
	return *this;
}

void Kmer_Short_Hash_Md5::rehash(){
//	cout << "[INFO] rehash length=" << _prime_list[prime_pos] << endl;
	if (prime_pos == _primes_num-1){
		cout << "[WARNING] reach rehash limit" << endl;
		return ;
	}
	int new_prime_pos = prime_pos + 1;
	Kmer_Short_Node_Ptr * tem_table = new Kmer_Short_Node_Ptr[ _prime_list[new_prime_pos] ];
	memset(tem_table, 0, _prime_list[new_prime_pos]*sizeof( Kmer_Short_Node_Ptr ));
	if (tem_table == NULL){
		cout << "bad_alloc" << endl;
		exit(1);
	}

	hash_table_used = 0;
	for(unsigned long i=0; i<_prime_list[prime_pos]; i++){
		if (hash_table[i] == NULL){
			continue;
		}
		Kmer_Short_Node_Ptr tem = hash_table[i];
		Kmer_Short_Node_Ptr next = tem->next;
		while(tem){
//			unsigned long pos = tem->get_md5() % (_prime_list[new_prime_pos]);
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
//	cout  << "[INFO] rehash end length=" << _prime_list[prime_pos] << endl;
	return ;
}

Kmer_Short_Node_Ptr Kmer_Short_Hash_Md5::has_kmer(const Kmer_Short &kmer){
	unsigned long pos = kmer.get_code() % _prime_list[prime_pos];
	Kmer_Short_Node_Ptr tem = hash_table[pos];
	while(tem){
		if (tem->kmer_short == kmer){
			return tem;
		}
		tem = tem->next;
	}
	return NULL;
}

static unsigned long cou = 1;
static unsigned long kmer_hashed = 0;
void Kmer_Short_Hash_Md5::push_a_kmer_short(const Kmer_Short &kmer){
	kmer_hashed ++;
	if (kmer_hashed % cou == 0){
		cout << "[INFO] " << kmer_hashed << " kmers have already been hashed..." << endl;
		cou *= 10;
	}
	Kmer_Short_Node_Ptr tem = has_kmer( kmer );
	if (tem){
		// tem->num ++;
		return ;
	}
	if ( (double)(hash_table_used + 1) >= _prime_list[ prime_pos ]*load_factor ){
		rehash();
	}
	if (last_alloc_num == last_alloc_used){
		last_alloc_num *= 1.5;
		last_alloc_num = last_alloc_num  > ONEG/2 ? ONEG/2 : last_alloc_num;

		Kmer_Short_Node_Ptr buffer = (Kmer_Short_Node_Ptr)malloc( last_alloc_num*sizeof(Kmer_Short_Node));
		cout << "[TEST] void Kmer_Short_Hash_Md5::push_a_kmer(unsigned long) : buffer malloc. (" << last_alloc_num*sizeof(Kmer_Short_Node)/(1024*1024) << "M) " << endl ;

		if (buffer == NULL){
			cout << "bad alloc" << endl;
			exit(1);
		}
//		last_alloc_num *= 1.5;
		last_alloc_used = 0;
		buffer_node.push_back(buffer);
		buffer_size.push_back(0);
		buffers ++;
	}
	tem = new( buffer_node[ buffers - 1 ]+last_alloc_used) Kmer_Short_Node(kmer);
	tem->pos = kmer_count;
	buffer_size[buffers - 1]++;
	last_alloc_used ++;
	kmer_count ++;
	unsigned long pos = tem->get_code() % _prime_list[prime_pos];
	if (hash_table[pos] == NULL){
		hash_table_used ++;
	}
	tem->next = hash_table[pos];
	hash_table[pos] = tem;
	return ;
}

void Kmer_Short_Hash_Md5::push_a_kmer_short(const Kmer &kmer, unsigned int numb){
	kmer_hashed ++;
	if (kmer_hashed % cou == 0){
		cout << "[INFO] " << kmer_hashed << " kmers have already been hashed..." << endl;
		cou *= 10;
	}
	Kmer_Short tmp = kmer.first_code >> 2;
	// kmer.get_pre_node();
	int next = int(kmer.get_code() & 0x3);
	Kmer_Short_Node_Ptr tem = has_kmer( tmp );
	if (tem){
		tem->set( next, numb );
		return ;
	}
	if ( (double)(hash_table_used + 1) >= _prime_list[ prime_pos ]*load_factor ){
		rehash();
	}
	if (last_alloc_num == last_alloc_used){
		last_alloc_num *= 1.5;
		last_alloc_num = last_alloc_num  > ONEG/2 ? ONEG/2 : last_alloc_num;

		Kmer_Short_Node_Ptr buffer = (Kmer_Short_Node_Ptr)malloc( last_alloc_num*sizeof(Kmer_Short_Node));
		cout << "[TEST] void Kmer_Short_Hash_Md5::push_a_kmer(unsigned long) : buffer malloc. (" << last_alloc_num*sizeof(Kmer_Short_Node)/(1024*1024) << "M) " << endl ;

		if (buffer == NULL){
			cout << "bad alloc" << endl;
			exit(1);
		}
//		last_alloc_num *= 1.5;
		last_alloc_used = 0;
		buffer_node.push_back(buffer);
		buffer_size.push_back(0);
		buffers ++;
	}
	tem = new( buffer_node[ buffers-1 ]+last_alloc_used) Kmer_Short_Node(tmp, next, numb);
	tem->pos = kmer_count;
	buffer_size[buffers - 1]++;
	last_alloc_used ++;
	kmer_count ++;
	unsigned long pos = tem->get_code() % _prime_list[prime_pos];
	if (hash_table[pos] == NULL){
		hash_table_used ++;
	}
	tem->next = hash_table[pos];
	hash_table[pos] = tem;
	return ;
}


unsigned short Kmer_Short_Hash_Md5::get_kmer_num(const Kmer_Short &kmer, int index){
	Kmer_Short_Node_Ptr tem = has_kmer(kmer);
	if (tem){
		return tem->num[index];
	}
	return 0;
}
unsigned short* Kmer_Short_Hash_Md5::get_kmer_nums(const Kmer_Short &kmer){
	Kmer_Short_Node_Ptr tem = has_kmer(kmer);
	if (tem){
		return tem->num;
	}
	return NULL;
}

Kmer_Short_Node_Ptr& Kmer_Short_Hash_Md5::operator[](unsigned long i){
	/*	if (i<0 || i>hash_table_length()){
			return NULL:
		}*/
	return hash_table[i];
}

// void Kmer_Short_Hash_Md5::filter(){
// 	cout << "[INFO] Filtering out Kmer_Short with coverage < 2 before filet kmer num=" << kmer_count << endl;
// 	for(unsigned long i=0; i<_prime_list[prime_pos]; i++){
// 		if (hash_table[i] == NULL){
// 			continue;
// 		}
// 		Kmer_Short_Node_Ptr tem = hash_table[i];
// 		while(tem){
// 			if (tem->num < 2){
// 				hash_table[i] = tem->next;
// 				tem = tem->next;
// 				kmer_count --;
// 			}else{
// 				break;
// 			}
// 		}
// 		if (hash_table[i] == NULL){
// 			continue;
// 		}
// 		Kmer_Short_Node_Ptr pre = hash_table[i];
// 		tem = pre->next;
// 		while(tem){
// 			if (tem->num < 2){
// 				pre->next = tem->next;
// 				tem = tem->next;
// 				kmer_count --;
// 			}else{
// 				pre = tem;
// 				tem = tem->next;
// 			}
// 		}
// 	}
// 	cout << "[INFO] after filter kmer num =" << kmer_count << endl;
// }

void Kmer_Short_Hash_Md5::show_all_kmers(){
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
		Kmer_Short_Node_Ptr tem = hash_table[i];
		while(tem){
			list_size ++;
		//	tem->show();
			tem = tem->next;
		}
	//	cout << "\t " << list_size << " kmers on this hash_node" << endl;
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

void Kmer_Short_Hash_Md5::esti_pos(){
	unsigned long tmp_first_code = 0, next_node = 0;
	unsigned long code = 0;
	
	for(unsigned long i=0; i<_prime_list[prime_pos]; i++){
		if (hash_table[i] ==NULL){
			continue;
		}
		Kmer_Short_Node_Ptr tem = hash_table[i];
		while(tem){
			code = tem->get_code();
			tmp_first_code = (code << (64 - 2*K + 4)) >> (64 - 2*K + 2);
			for(int j = 0; j < 4; j++)
			{
				if(tem->num[j] == 0)
					continue;
				// tmp_first_code = tmp_first_code >> (64 - 2*K + 2);
				next_node = (tmp_first_code | j);

				tem->nex[j] = has_kmer(next_node)->pos;
			}
			tem = tem->next;
		}
	}
}


void  Kmer_Short_Hash_Md5::transform2graph(vector<Node>& ret)
{
	// g = vector<Node>();
	int size = ret.size();
	for(int i=0; i<buffers; i++){
		int n = buffer_size[i];
		Kmer_Short_Node_Ptr p = buffer_node[i];
		for(int j=0; j<n; j++){
			ret.push_back(Node(*(p+j)));
		}
		free(p);		
		cout << "[TEST] short kmer hash : free buffer " << i << endl;

	}
	buffers = 0;
	buffer_node.clear();
	buffer_size.clear();
	// return ret;
}