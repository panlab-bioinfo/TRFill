# include "statistic_test.hpp"
# include <iostream>
# include "clipp.h"
# include <vector>
using namespace clipp; using std::cout; using std::string; using std::vector;


int main(int argc, char *argv[]){
    string kmerfile;        // jellyfish result dump
    string ref_fasta_file;
    vector<string> qry_files;
    string paf_file;
    string output = "./";          //output path
    string chr;
    int t = 8;
    int k_size = 21;        //k-mer length
    uint32_t gap_start,gap_end;
    double my_delta = 0.963534; // 0.963368
    double my_sigma = 336.210817; // 265.221971;
    double significance,z_value;

    auto cli = (
        option("-t", "--threads").doc("Number of threads [8]") & value("threads", t),
        option("-k", "--kmer").doc("K-mer length [21]") & value("k-size", k_size),
        option("-o").doc("Output path [./]") & value("outpath", output),
		option("-d", "--delta").doc("delta [0.963368]") & value("delta", my_delta),
		option("-s", "--sigma").doc("sigma [265.221971]") & value("sigma", my_sigma),
        value("target chr", chr).doc("target chr on reference with gap to patch"),
        value("gap_start", gap_start).doc("gap start coordinate on reference(0 based)"),
        value("gap_end", gap_end).doc("gap end coordinate on reference(0 based)"),
        value("significance", significance).doc("significance level for statistic test"),
        value("reference genome",ref_fasta_file).doc("reference genome fasta file path"),
        value("paf result", paf_file).doc("paf file path of qry aligning to ref"),
        value("uniquekmer file", kmerfile).doc("unique Dump file of jellyfish result"),
        value("reads sequence", qry_files).doc("Query read sequence files path").repeatable(true)
    );

    auto fmt = doc_formatting{}
		.first_column(8)                           //left border column for text body
		.doc_column(30)                            //column where parameter docstring starts
		.last_column(100);

    if(!parse(argc, const_cast<char **>(argv), cli)) {
		cout << "Usage:\n" << usage_lines(cli, "statistic")
		<< "\nOptions:\n" << documentation(cli,fmt) << "\nERROR: Required parameter missing\n";
		// throw "Division by zero condition!";
		exit(0);
    }

    read_kmer(kmerfile.c_str(),t,k_size);
    build_fasta(ref_fasta_file.c_str(), chr.c_str(), gap_start, gap_end, qry_files, paf_file.c_str(), output, k_size, significance, my_delta, my_sigma, 1, t);
    return 0;

    
}
