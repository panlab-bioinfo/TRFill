#include <string>
#include <vector>
using namespace std;


int read_kmer(const char *k_path, int t,int k);
// search the poses for certain k-mers; Using 64 bit binary int to store the k-mer 

int build_pos(vector<string> fasta_file, string out_path, int k, bool strand, int t1, int t2);