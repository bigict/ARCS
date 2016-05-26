#include "gap_filling.h"
#include "gap_filler.h"
#include "constant.h"

#include <string>
#include <fstream>
#include <tuple>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("arcs.GapFilling"));

int GapFilling::run(const Properties& options, const Arguments& arguments) {
    int r = 0;
    if ((r = checkOptions(options)) != 0) {
        return r;
    }

    boost::filesystem::path workdir(options.get< std::string >("d", "."));
    if (!boost::filesystem::exists(workdir) && !boost::filesystem::create_directory(workdir)) {
        LOG4CXX_ERROR(logger, boost::format("failed to create directory: %s") % workdir);
        return 1;
    }

    LOG4CXX_DEBUG(logger, "gap_filling begin");
		
    std::string initial_contig_file_name;
    std::string scaffold_file_name;
    std::string condensed_contig_file_name;
    std::string work_dir;
    size_t K = options.get< size_t >("K", kKmerSize);
    size_t OVERLAP = 10;
    if (options.find("MAX_OVERLAP") != options.not_found()) {
        OVERLAP = options.get< size_t >("MAX_OVERLAP");
    }
    if (options.find("C") != options.not_found()) {
        condensed_contig_file_name = options.get< std::string >("C");
    }
    if (options.find("l") != options.not_found()) {
        scaffold_file_name = options.get< std::string >("l");
    }
    if (options.find("I") != options.not_found()) {
        initial_contig_file_name = options.get< std::string >("I");
    }
	
    size_t INSERT_SIZE = options.get< size_t >("INSERT_SIZE");
    double DELTA = options.get< double >("DELTA");

	GapFiller gf(K, OVERLAP, INSERT_SIZE, DELTA);
    // load data
    {
        typedef bool(GapFiller::*LoadDataPtr)(const std::string&);
        typedef std::tuple< std::string, LoadDataPtr > File2FuncPtr;
        std::vector< File2FuncPtr > file2func_list = boost::assign::list_of
            (std::make_tuple(condensed_contig_file_name, &GapFiller::input_contigs))
            (std::make_tuple(initial_contig_file_name,   &GapFiller::input_debruijn))
            (std::make_tuple(scaffold_file_name,         &GapFiller::input_scaffold))
            ;
        BOOST_FOREACH(const File2FuncPtr& file2func, file2func_list) {
            std::string file = std::get< 0 >(file2func);
            if (!boost::bind(std::get< 1 >(file2func), &gf, _1)(file)) {
                LOG4CXX_ERROR(logger, boost::format("failed to load %s") % file);
                return 2;
            }
        }
    }

	gf.fill();

    // write data
    {
        std::string file = boost::str(boost::format("%dmer.scaf_seq_with_gaps") % K);
        boost::filesystem::ofstream stream(workdir / file);
        if (!stream) {
            LOG4CXX_ERROR(logger, boost::format("Create %s error!!!") % file);
            return 2;
        }
        stream << gf;
    }

    LOG4CXX_DEBUG(logger, "gap_filling end");

    return 0;
}

GapFilling GapFilling::_runner;

GapFilling::GapFilling() : Runner("c:s:d:K:O:C:I:l:h", boost::assign::map_list_of('O', "MAX_OVERLAP")) {
    RUNNER_INSTALL("gapfill", this, "fill intra-scaffold gaps");
}

int GapFilling::printHelps() const {
    std::cout << "usage: arcs gapfill [arguments]" << std::endl;
    std::cout << std::endl;
    std::cout << "\t-c[=<file>]    log config file, default ./log4cxx.properties" << std::endl;
    std::cout << "\t-s[=<file>]    scaffold_parameter_0" << std::endl;
    std::cout << "\t-d[=<workdir>] word dir, default ." << std::endl;
    std::cout << "\t-K[=<number>]  kmer size, default 31" << std::endl;
    std::cout << "\t-O[=<number>]  minimum same bases number to judge overlap, default 10" << std::endl;
    std::cout << "\t-C[=<file>]    contig file <cdbg_copy_number.fa>" << std::endl;
    std::cout << "\t-l[=<file>]    <component_last>" << std::endl;
    std::cout << "\t-I[=<file>]    <condensed_de_bruijn_graph_before_trimming.data>" << std::endl;
    std::cout << "\t-h             help" << std::endl;
    std::cout << std::endl;
    return 256;
}

int GapFilling::checkOptions(const Properties& options) const {
    if (options.find("h") != options.not_found()) {
        return printHelps();
    }
    return 0;
}
