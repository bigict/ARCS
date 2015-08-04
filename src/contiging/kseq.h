#ifndef kseq_h_
#define kseq_h_

#include <string>
#include <iostream>

//
// DNASeq represents a DNA sequence.
//
class DNASeq {
public:
    DNASeq() {}
    DNASeq(const std::string& name, const std::string& seq, const std::string& quality) : name(name), seq(seq), quality(quality) {}
    virtual ~DNASeq() {}

    std::string name;
    std::string seq;
    std::string quality;
private:
    friend class DNASeqReader;
    friend std::ostream& operator << (std::ostream& os, const DNASeq& seq);
};

//
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
// 
class DNASeqReader {
public:
    DNASeqReader(std::istream& stream, double percent=1.0, size_t read_cutoff=-1) : _stream(stream), _percent(percent), _read_cutoff(read_cutoff) {
    }
    
    bool read(DNASeq& sequence);
private:
    void cutoff(DNASeq& sequence) const;

    std::istream& _stream;
    double _percent;
    size_t _read_cutoff;
};

#endif // kseq_h_
