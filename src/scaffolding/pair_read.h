#ifndef pair_read_h__
#define pair_read_h__

#include <string>
#include <iostream>
#include <vector>

#include "kseq.h"
#include "kmer.h"
#include "kmer_tbl.h"
#include "graph.h"

struct PairRead {
    std::string set1;
    std::string set2;
};

class PairReadSet {
public:
    PairReadSet(std::istream& stream_1, std::istream& stream_2, size_t K, size_t INSERT_SIZE);
    void readFile();
    void buildConnectGraph(Graph& g, KmerTable& tbl, const std::vector<Component>& component);
    void estimateInsertSize(const KmerTable& tbl);
    size_t size() {
        return _pair_reads.size(); 
    }

    size_t INSERT_SIZE;
    double DELTA;
    size_t PAIR_KMER_NUM;

private:
    void estimateOnePR(const std::string& read1, const std::string& read2, std::vector<int>& insert_len_dis, const KmerTable& tbl);
    std::string make_complement(std::string seq);
    void findLink(const std::string& read1, const std::string& read2, Graph& graph, const KmerTable& tbl, const std::vector<Component>& components);    
    friend std::ostream& operator<<(std::ostream& os, const PairReadSet& p_r);
    std::vector< PairRead > _pair_reads;
    std::istream& _stream_1;
    std::istream& _stream_2;
    size_t _k;
};

//pair read file is two FASTQ Format file
// FASTQ Format Specification
// Syntax
//    <fastq>	:=	<block>+
//    <block>	:=	@<seqname>\n<seq>\n+[<seqname>]\n<qual>\n
//    <seqname>	:=	[A-Za-z0-9_.:-]+
//    <seq>	:=	[A-Za-z\n\.~]+
//    <qual>	:=	[!-~\n]+
// Requirements
//    1. The <seqname> following '+' is optional, but if it appears right after '+', it should be 
//    identical to the <seqname> following '@'.
//    2. The length of <seq> is identical the length of <qual>. Each character in <qual> represents 
//    the phred quality of the corresponding nucleotide in <seq>.
//    3. If the Phred quality is $Q, which is a non-negative integer, the corresponding quality 
//    character can be calculated with the following Perl code:
//        $q = chr(($Q<=93? $Q : 93) + 33);
//    where chr() is the Perl function to convert an integer to a character based on the ASCII  
//    table.
//    4. Conversely, given a character $q, the corresponding Phred quality can be calculated with:
//        $Q = ord($q) - 33;
//    where ord() gives the ASCII code of a character.

#endif
