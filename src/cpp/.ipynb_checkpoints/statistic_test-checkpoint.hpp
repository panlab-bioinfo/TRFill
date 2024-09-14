#include <string>
#include <vector>
using namespace std;


int read_kmer(const char *k_path, int t,int k);
// search the poses for certain k-mers; Using 64 bit binary int to store the k-mer 

int build_fasta(const char *ref_fasta_file, const char *chr, uint32_t gap_start, uint32_t gap_end, vector<string> qry_files, const char *paf_file, string out_path, int k, double significance, double my_delta, double my_sigma, int t1, int t2);