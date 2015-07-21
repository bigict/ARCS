#ifndef _CONSTANT_H
#define _CONSTANT_H

#include <iostream>
#include <string>

using namespace std;

const int ALPHABET = 5;               //alphabet
const long MAX_INSERT_LEN = 10000;    //max insert len(for next-generation sequencing data)
const string R_NEXT = "=";            //the next pair of paired reads
const double ISNAN = 1e-320;          //minimum probability

//constant index of alphabet
const char BASE[ALPHABET] = {'A','T','C','G','N'}; 
const int A = 0;
const int T = 1;
const int C = 2;
const int G = 3;
const int N = 4;

const int DELTA = 0;                  //read alignment sequence length delta

enum F_OR_R {FORWARD,REVERSE};        //paired reads or its reverse complement

const int MAX_N_RATIO = 0.8;          //the max ratio of 'N' in correct reads 

//const string NO_ERROR_MD("MD:Z:36");

//sam segments
const int SIZE_OF_SAM = 13;
enum SEG_INDEX{
	READ_NAME=0, 
	LABEL,
	SEQ_NAME,
	LEFT_END_POS, 
	MAP_QUALITY, 
	CIGAR, 
	R_NEXT_INDEX,
	RIGHT_END_POS,
	INSERT_LEN, 
	READ_SEQ,
	QUALITY,
	NM_CHAR,
	MD_CHAR
};

const char * const NO_SEQ_LABEL = "*";
const char * const DELIM_OF_READ = "IDMS^\t\n	";
const char * const DELIM_OF_REF = "ACGTN^\t\n	";

#endif /*_CONSTANT_H*/ 
