
#include "Kmer_Hash.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <new>

using namespace std;

template <typename T>
Kmer_Hash<T>::Kmer_Hash(unsigned int half_first_alloc):hash_table(NULL), last_alloc_num(half_first_alloc), last_alloc_used(half_first_alloc), prime_pos(7), kmer_count(0),load_factor(0.75), hash_table_used(0){
	
	hash_table = new T*[ _prime_list[prime_pos] ];
	cout << "[INFO] construct Kmer_Hash" << endl;
}

template <typename T>
unsigned long Kmer_Hash<T>::size(){
	return kmer_count;
}

template <typename T>
unsigned long Kmer_Hash<T>::hash_table_length(){
	return _prime_list[prime_pos];
}

template <typename T>
void Kmer_Hash<T>::clear(){
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
	cout << "[INFO] Hash cleared." << endl;
	return ;
}

template <typename T>
Kmer_Hash<T>::~Kmer_Hash(){
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

template <typename T>
Kmer_Hash<T>::Kmer_Hash(const Kmer_Hash<T> &kh){
	cout << "can not be called" << endl;
}

template <typename T>
Kmer_Hash<T>& Kmer_Hash<T>::operator=(const Kmer_Hash<T> &kh){
	cout << "can not be called" << endl;
	return *this;
}

template <typename T>
void Kmer_Hash<T>::rehash(){
//	cout << "[INFO] rehash length=" << _prime_list[prime_pos] << endl;
	if (prime_pos == _primes_num-1){
		cout << "[WARNING] reach rehash limit" << endl;
		return ;
	}
	int new_prime_pos = prime_pos + 1;
	T* * tem_table = new T*[ _prime_list[new_prime_pos] ];
	memset(tem_table, 0, _prime_list[new_prime_pos]*sizeof( T* ));
	if (tem_table == NULL){
		cout << "bad_alloc" << endl;
		exit(1);
	}

	hash_table_used = 0;
	for(unsigned long i=0; i<_prime_list[prime_pos]; i++){
		if (hash_table[i] == NULL){
			continue;
		}
		T* tem = hash_table[i];
		T* next = tem->next;
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
//	cout  << "[INFO] rehash end length=" << _prime_list[prime_pos] << endl;
	return ;
}

template <typename T>
T* Kmer_Hash<T>::has_kmer(const Kmer &kmer){
	unsigned long pos = kmer.get_code() % _prime_list[prime_pos];
	T* tem = hash_table[pos];
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
template <typename T>
void Kmer_Hash<T>::push_a_kmer(const Kmer &kmer,int i, int j, int k){
//	cout << "[INFO] push_a_kmer length=" << _prime_list[prime_pos] << " kmer_num=" << kmer_count << " last_alloc=" << last_alloc_num << " last_alloc_used=" << last_alloc_used << endl;
	kmer_hashed ++;
//	if (kmer_hashed % cou == 0){
//		cout << "[INFO] " << kmer_hashed << " kmers have already been hashed..." << endl;
//		cou *= 10;
//	}
	T* tem;// = has_kmer( kmer );

	if ( (double)(hash_table_used + 1) >= _prime_list[ prime_pos ]*load_factor ){
		rehash();
	}
	if (last_alloc_num == last_alloc_used){
		T* buffer = (T*)malloc( last_alloc_num*2*sizeof(T));
		if (buffer == NULL){
			cout << "bad alloc" << endl;
			exit(1);
		}
		last_alloc_num *= 2;
		last_alloc_used = 0;
		buffer_node.push_back(buffer);
	}
	tem = new( buffer_node[ buffer_node.size()-1 ]+last_alloc_used) T(kmer,i,j,k);
	last_alloc_used ++;
	kmer_count ++;
	unsigned long pos = tem->get_code() % _prime_list[prime_pos];
	if (hash_table[pos] == NULL){
		hash_table_used ++;
	}
	tem->next = hash_table[pos];
	hash_table[pos] = tem;
//	cout << "[INFO] push_a_kmer length=" << _prime_list[prime_pos] << " kmer_num=" << kmer_count << " last_alloc=" << last_alloc_num << " last_alloc_used=" << last_alloc_used << endl;
	return ;
}

template <typename T>
void Kmer_Hash<T>::push_a_path(int a,int b,int c){
		kmer_hashed ++;
		T* tem;
		if(c==-1){
			tem = has_path2(a,b);
		}else
			tem = has_path3(a,b,c );

		if(tem){
			tem->num ++;
			return;
		}

		if ( (double)(hash_table_used + 1) >= _prime_list[ prime_pos ]*load_factor ){
			rehash();
		}
		if (last_alloc_num == last_alloc_used){
			T* buffer = (T*)malloc( last_alloc_num*2*sizeof(T));
			if (buffer == NULL){
				cout << "bad alloc" << endl;
				exit(1);
			}
			last_alloc_num *= 2;
			last_alloc_used = 0;
			buffer_node.push_back(buffer);
		}
		tem = new( buffer_node[ buffer_node.size()-1 ]+last_alloc_used) T(a,b,c);
		last_alloc_used ++;
		kmer_count ++;
		unsigned long pos = tem->B % _prime_list[prime_pos];
		if (hash_table[pos] == NULL){
			hash_table_used ++;
		}
		tem->next = hash_table[pos];
		hash_table[pos] = tem;
		return ;
}
template <typename T>
T* Kmer_Hash<T>::has_path1(int b,int c){
	unsigned long pos = b % _prime_list[prime_pos];
	T* tem = hash_table[pos];
	while(tem){
		if (tem->C == c){
			return tem;
		}
		tem = tem->next;
	}
	return NULL;
}
template <typename T>
T* Kmer_Hash<T>::has_path2(int a,int b){
	unsigned long pos = b % _prime_list[prime_pos];
	T* tem = hash_table[pos];
	while(tem){
		if (tem->A == a){
			return tem;
		}
		tem = tem->next;
	}
	return NULL;
}
template <typename T>
T* Kmer_Hash<T>::has_path3(int a,int b,int c){
	unsigned long pos = b % _prime_list[prime_pos];
	T* tem = hash_table[pos];
	while(tem){
		if (tem->A == a && tem->C == c){
			return tem;
		}
		tem = tem->next;
	}
	return NULL;
}

template <typename T>
int* Kmer_Hash<T>::get_kmer_pos(const Kmer &kmer){
	T* tem = has_kmer(kmer);
	if (tem){
		return &(tem->pos[0]);
	}
	return NULL;
}

template <typename T>
T* Kmer_Hash<T>::operator[](unsigned long i){
	return hash_table[i];
}

template <typename T>
void Kmer_Hash<T>::show_all_kmers(){
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
		T* tem = hash_table[i];
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

