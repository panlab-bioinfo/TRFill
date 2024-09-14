#include <map>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <bitset>
#include <algorithm>
#include <set>
#include <zlib.h>
#include <limits.h>
#include <time.h>
#include <boost/math/distributions/normal.hpp>
#include "thread_pool.hpp"
#include "uthash.h"
#include "statistic_test.hpp"
#include "stringtools.hpp"

using namespace std;

typedef struct my_struct {
    uint64_t key;                    /* key */
    //char name;
    UT_hash_handle hh;         /* makes this structure hashable */
} hash_node;

class PAF_info{
public:
    uint32_t qry_start;
    uint32_t qry_end;
    string strand;
    uint32_t ref_start;
    uint32_t ref_end;

    PAF_info(vector<string> infos){
        qry_start = atoi(infos[2].c_str());
        qry_end = atoi(infos[3].c_str());
        strand = infos[4];
        ref_start = atoi(infos[7].c_str());
        ref_end = atoi(infos[8].c_str());
    }
};

class hash_set{
public:
    hash_node *dict=NULL;
    void insert(uint64_t *key){
        hash_node *s;
        HASH_FIND_INT(dict, key, s);
        if (s == NULL)
        {
            s = (hash_node*)malloc(sizeof *s);
            s->key = *key;
            HASH_ADD_INT(dict, key, s);
        }
    }
    hash_node *find(uint64_t *key)
    {
        hash_node *tmp=NULL;
        HASH_FIND_INT(dict, key, tmp);
        return tmp;
    }
    int size() {
        return HASH_COUNT(dict);
    }
};

mutex m;
vector<hash_set> kmer;
vector<uint32_t> ref_pos;
vector<uint64_t> ref_kmers;
map<string, bool> read_writen;
map<string, vector<PAF_info>> PAF_infos;
double delta;
double sigma;
int k_size; // k-mer length


pair<bool, uint32_t> BinarySearch(uint32_t key){
    int low = 0, mid = 0, high = ref_pos.size() - 1;
    while(low <= high){
        mid = (low + high) / 2;
        if(ref_pos[mid] == key)
            return make_pair(true,mid);
        else if(ref_pos[mid] < key)
            low = mid + 1;
        else
            high = mid - 1;
    }
    return make_pair(false,low);
}

// 读取k-mer为二进制整数
uint64_t ctoi(char c)
{
    switch (c)
    {
    case 'A':
        return 0ull;
    case 'a':
        return 0ull;
    case 'T':
        return 3ull;
    case 't':
        return 3ull;
    case 'C':
        return 1ull;
    case 'c':
        return 1ull;
    case 'G':
        return 2ull;
    case 'g':
        return 2ull;
    default:
        return 4ull;
    }
}

//
uint64_t tobin(string *s)
{
    uint64_t k_value = 0;
    int i = (k_size - 1) * 2;
    for (auto c = s->begin(); c < s->end(); c++)
    {
        uint64_t v = ctoi(*c);
        if (v != 4)
        {
            k_value = k_value | (v << i);
            i -= 2;
        }
    }
    return k_value;
}
// reverse the bin of kmer
uint64_t reversebin(uint64_t x)
{
    uint64_t mid = x & 0x300000ull;
    x = (x & 0x3ffffc00000ull) >> 22 | ((x & 0xfffffull) << 22);
    x = (x & 0x3ff000ffc00ull) >> 10 | (x & 0xffc003ffull) << 10;
    uint64_t quarter = x & 0x300c00c030ull;
    x = (x & 0x3C0F00F03C0ull) >> 6 | (x & 0xF03C03C0Full) << 6;
    x = (x & 0x30CC30C330Cull) >> 2 | (x & 0xC330C30CC3ull) << 2;
    x = x | mid | quarter;
    return (~x) & 0x3ffffffffffull;
}

std::string reverse_complement(const std::string& dnaSequence) {
    std::string complement = dnaSequence;

    std::transform(complement.begin(), complement.end(), complement.begin(), [](char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'a':
                return 't';
            case 't':
                return 'a';
            case 'c':
                return 'g';
            case 'g':
                return 'c';
            default:
                return base;
        }
    });

    std::reverse(complement.begin(), complement.end());

    return complement;
}

int creat_dict(const char *k_path, uint64_t chunk_line_num, int idx)
{
    ifstream k_mer(k_path);
    hash_set dict;
    uint64_t k_value = 0;
    string s,temp_kmer;
    uint64_t line_count = 0;
    while(line_count < idx * chunk_line_num){
        getline(k_mer,s);
        line_count += 1;
    }
    line_count = 0;
    while(line_count < chunk_line_num && !k_mer.eof()){
        getline(k_mer,s);
        if(s == ""){
            continue;
        }
        temp_kmer = s.substr(0,k_size);
        k_value = tobin(&temp_kmer);
        // k_value = stoull(s.substr(0,(k_size + 1) / 2).c_str(),0,16);
        dict.insert(&k_value);
        line_count += 1;
    }
    cout<< idx << "\t" << dict.size() << "\n";
    m.lock();
    kmer.push_back(dict);
    m.unlock();
    k_mer.close();
    return 1;
}

int read_kmer(const char *k_path, int t ,int k)
{
    /*
    将文件分割为与线程数一样多的字符串块；
    每块文件获得一个线程，对其进行二进制转换并存入unorderset集合
    */
    ThreadPool pool(t);
    pool.init();
    ifstream k_mer(k_path);
    k_size = k;

    // 分块
    string s;
    uint64_t line_count = 0;
    while(!k_mer.eof()){
        getline(k_mer,s);
        line_count += 1;
    }
    k_mer.close();
    uint64_t chunk_line_num = (line_count + t - 1) / t;
    for (int i = 0; i < t; i++)
    {
        pool.submit(creat_dict, k_path, chunk_line_num, i);
    }

    pool.shutdown();
    return 1;
}

/* store the sbin(bin of kmer),rbin(reverse bin of kmer), pos at contig and contig name*/
class KMER
{
public:
    uint64_t sbin;
    uint32_t pos;
    set<uint64_t> kmers;
    int is_ukmer_reference(int t)
    {
        for(int i = 0;i < t;i ++){
            if (kmer[i].find(&sbin) != NULL)
            {
                ref_pos.push_back(pos);
                ref_kmers.push_back(sbin);
                break;
            }
        }
        return 0;
    }
    int is_ukmer_reads(int t)
    {
        for(int i = 0;i < t;i ++){
            if (kmer[i].find(&sbin) != NULL)
            {
                kmers.emplace(sbin);
                break;
            }
        }
        return 0;
    }
};

int search_kmer_reference(string line, uint32_t start, uint32_t end, long bin_mask, int t1, int n)
{
    KMER k;
    uint64_t base;
    string qkmer = line.substr(start, k_size);
    k.sbin = tobin(&qkmer);  
    for (uint32_t i = start + k_size; i < end + 1; i++)
    {
        k.pos = i - k_size;

        k.is_ukmer_reference(n);

        base = ctoi(line[i]);
        k.sbin = (k.sbin << 2 | base) & bin_mask;
    }
    k.pos++;
    k.is_ukmer_reference(n);

    return 1;
}

int process_reads_line(FILE *fp, int idx, string line, uint64_t line_count, mutex qry_m[], bool *flag, string *name, vector<PAF_info> *temp_infos, int round, long bin_mask, double z_value, int n){
    vector<string> ids;    
    string strand_line;

    if(line_count % round == 0){
        ids = split_find(line, " ");
        *name = ids[0].substr(1);
        if(PAF_infos.find(*name) != PAF_infos.end()){
            *temp_infos = PAF_infos[*name];
            read_writen[*name] = false;
            *flag = true;
        }
    }
    else if(line_count % round == 1 && *flag){
        for (auto iter = (*temp_infos).begin(); iter != (*temp_infos).end(); iter++){
            cout << idx << "\t" << line_count << "\t" << *name << "\t" << (*iter).ref_start << "\t" << (*iter).ref_end << "\t" << (*iter).strand << "\t" << (*iter).qry_start << "\t" << (*iter).qry_end << endl;
            if((*iter).strand == "-")
                strand_line = reverse_complement(line);
            else
                strand_line = line;
            KMER k;
            uint64_t base;
            string qkmer = strand_line.substr((*iter).qry_start, k_size);
            k.sbin = tobin(&qkmer);
            for (uint32_t i = (*iter).qry_start + k_size; i < (*iter).qry_end + 1; i++)
            {
                k.is_ukmer_reads(n);

                base = ctoi(strand_line[i]);
                k.sbin = (k.sbin << 2 | base) & bin_mask;
            }
            k.is_ukmer_reads(n);
            cout << idx << "\t" << line_count << "\t" << *name << "\tk.count:" << k.kmers.size() << endl;

            pair<bool,uint32_t> result_start = BinarySearch((*iter).ref_start);
            cout << idx << "\t" << line_count << "\t" << *name << "\tresult_start:" << result_start.first << "\t" << result_start.second << endl;
            pair<bool,uint32_t> result_end = BinarySearch((*iter).ref_end - k_size + 1);
            cout << idx << "\t" << line_count << "\t" << *name << "\tresult_end:" << result_end.first << "\t" << result_end.second << endl;
            result_end.second = result_end.first ? result_end.second : result_end.second - 1;
            set<uint64_t> ref_kmer;
            for(int i = result_start.second; i <= result_end.second; i++){
                ref_kmer.emplace(ref_kmers[i]);
            }
            set<uint64_t> common_qry_kmers;
            set_intersection(k.kmers.begin(),k.kmers.end(),ref_kmer.begin(),ref_kmer.end(),inserter(common_qry_kmers,common_qry_kmers.begin()));
            cout << idx << "\t" << line_count << "\t" << *name << "\t" << common_qry_kmers.size() << "\t" << ref_kmer.size() << "\t" << ref_kmer.size() * delta - z_value * sigma << "\t" << ref_kmer.size() * delta + z_value * sigma << endl;
            if(ref_kmer.size() * delta - z_value * sigma <= common_qry_kmers.size() && common_qry_kmers.size() <= ref_kmer.size() * delta + z_value * sigma){
                cout << idx << "\t" << line_count << "\t" << *name << read_writen[*name] << endl;
                if(!read_writen[*name]){
                    cout << idx << "\t" << line_count << "\t" << *name << "\tline_length:" << strlen(line.c_str()) << endl;
                    qry_m[idx].lock();
                    fprintf(fp, ">%s\n", (*name).c_str());
                    fprintf(fp, "%s\n", strand_line.c_str());
                    qry_m[idx].unlock();
                    read_writen[*name] = true;
                    break;
                }
            }
        }
        *flag = false;
    }
    return 1;
}

int search_kmer_reads(const char *qry_file, uint64_t chunk_line_num, int idx, int round, long bin_mask, double z_value, FILE *file[], mutex qry_m[], int n)
{
    FILE *fp = file[idx];
    string line,name;
    vector<PAF_info> temp_infos;
    
    if(qry_file[strlen(qry_file) - 1] == 'z'){
        uint64_t line_count = 0;
        bool flag = false;
        gzFile gzfile = gzopen(qry_file, "rb");
        if(gzfile == nullptr){
            std::cerr << "Failed to open GZ file: " << qry_file << std::endl;
            return 1;
        }

        const int BUFFER_SIZE = 1024;
        char buffer[BUFFER_SIZE];

        ostringstream lineStream;

        while(true){
            if(line_count >= (idx + 1) * chunk_line_num){
                break;
            }

            int bytesRead = gzread(gzfile, buffer, BUFFER_SIZE);

            if(bytesRead <= 0){
                if(gzeof(gzfile)){
                    // handle the last line without '\n'
                    line = lineStream.str();
                    if(!line.empty() && line_count >= idx * chunk_line_num){
                        process_reads_line(fp, idx, line, line_count, &qry_m[0], &flag, &name, &temp_infos, round, bin_mask, z_value, n);
                        line_count += 1;
                    }
                }
                else{
                    std::cerr << "Error reading GZ file: " << qry_file << std::endl;
                }

                break;
            }

            lineStream.write(buffer,bytesRead);

            while(true){
                size_t pos = lineStream.str().find('\n');
                if(pos == std::string::npos){
                    // no '\n'
                    break;
                }

                line = lineStream.str().substr(0, pos);
                // process line
                if(line_count >= idx * chunk_line_num){
                    process_reads_line(fp, idx, line, line_count, &qry_m[0], &flag, &name, &temp_infos, round, bin_mask, z_value, n);
                }
                line_count += 1;

                lineStream.str(lineStream.str().substr(pos + 1));
                lineStream.seekp(0,ios_base::end);
            }
        }

        gzclose(gzfile);
        return 1;
    }
    else if(qry_file[strlen(qry_file) - 1] == 'a' || qry_file[strlen(qry_file) - 1] == 'q'){
        ifstream qry(qry_file);

        uint64_t line_count = 0;
        while(line_count < idx * chunk_line_num){
            getline(qry,line);
            line_count += 1;
        }
        cout << idx << "\t" << line_count << endl;   

        bool flag = false;
        line_count = 0;
        while(line_count < chunk_line_num && !qry.eof()){
            getline(qry,line);
            if(line == "")
                break;
            process_reads_line(fp, idx, line, line_count, &qry_m[0], &flag, &name, &temp_infos, round, bin_mask, z_value, n);
            line_count += 1;
        }

        qry.close();

        return 1;             
    }
    else{
        cout << "False format!\n";
        exit(EXIT_FAILURE);  
    }
}

int build_fasta(const char *ref_fasta_file, const char *chr, uint32_t gap_start, uint32_t gap_end, vector<string> qry_files, const char *paf_file, string out_path, int k, double significance, double my_delta, double my_sigma, int t1, int t2)
{
    cout << "Number of threads: " << t2 << "\n"
         << "Reading the kmer to dicts\n";
    ThreadPool fpool(t2);
    fpool.init();
    FILE *file[t2];
    string line, name;
    vector<string> ids;
    k_size = k;
    delta = my_delta;
    sigma = my_sigma;
    int round = 0;
    uint64_t line_count;

	boost::math::normal_distribution<> normal(0.0, 1.0);
	double z_value = quantile(normal, 1 - significance / 2);

    /*
        Creat mutiple files for mutiple threads;
        File pointers are store in a FILE * array.
    */
    char outname[15];
    strcpy(outname, "/tempAA.fasta");
    string mk = "mkdir -p " + out_path;
    system(mk.c_str());

    for (int i = 0; i < t2; i++)
    {
        char *outfile = (char *)malloc((out_path.size() + 25) * sizeof(char));
        strcpy(outfile, out_path.c_str());
        strcat(outfile, outname);
        file[i] = fopen(outfile, "w");
        if(outname[6] == 'Z'){
            outname[5]++;
            outname[6] = 'A';            
        }
        else
            outname[6]++;
    }
    
    long bin_mask = 0x3;
    for(int i = 0;i < k_size - 1;i++){
        bin_mask = (bin_mask << 2) | bin_mask;
    }


    uint32_t left_most = INT_MAX;
    uint32_t right_most = 0;
    cout << "Reading paf file!\n";
    ifstream paf(paf_file);
    while(getline(paf,line)){
        vector<string> infos = split_find(line, "\t");
        if(strcmp(infos[5].c_str(),chr) == 0 && max(gap_start, uint32_t(atoi(infos[7].c_str()))) <= min(gap_end, uint32_t(atoi(infos[8].c_str())))){
            PAF_info info(infos);
            PAF_infos[infos[0]].push_back(info);
            left_most = min(left_most,uint32_t(atoi(infos[7].c_str())));
            right_most = max(right_most,uint32_t(atoi(infos[8].c_str())));
        }
    }
    cout << "Reading paf file done!\n";
    paf.close();
    left_most = min(left_most,gap_start);
    right_most = max(right_most,gap_end);
    // cout << left_most << "\t" << right_most << endl;


    cout << "Searching kmer in reference file!\n";
    ifstream ref_fa(ref_fasta_file); // TODO: multi line fa, fq, gz
    if(ref_fasta_file[strlen(ref_fasta_file) - 1] == 'a')
        round = 2;
    else if(ref_fasta_file[strlen(ref_fasta_file) - 1] == 'z'){
        line_count = 0;
        gzFile gzfile = gzopen(ref_fasta_file, "rb");
        if(gzfile == nullptr){
            std::cerr << "Failed to open GZ file: " << ref_fasta_file << std::endl;
            return 1;
        }

        const int BUFFER_SIZE = 1024;
        char buffer[BUFFER_SIZE];

        ostringstream lineStream;

        while(true){
            int bytesRead = gzread(gzfile, buffer, BUFFER_SIZE);

            if(bytesRead <= 0){
                if(gzeof(gzfile)){
                    // handle the last line without '\n'
                    line = lineStream.str();
                    if(!line.empty()){
                        if (line_count % 2 == 0)
                        {
                            ids = split_find(line, " ");
                            name = ids[0].substr(1);
                        }
                        else if(line_count % 2 == 1 && strcmp(name.c_str(), chr) == 0){
                            search_kmer_reference(line, left_most, right_most, bin_mask, t1, t2);
                            break;
                        }
                        line_count += 1;
                    }
                }
                else{
                    std::cerr << "Error reading GZ file: " << ref_fasta_file << std::endl;
                }

                break;
            }

            lineStream.write(buffer,bytesRead);

            while(true){
                size_t pos = lineStream.str().find('\n');
                if(pos == std::string::npos){
                    // no '\n'
                    break;
                }

                line = lineStream.str().substr(0, pos);
                // process line
                if (line_count % 2 == 0)
                {
                    ids = split_find(line, " ");
                    name = ids[0].substr(1);
                }
                else if(line_count % 2 == 1 && strcmp(name.c_str(), chr) == 0){
                    search_kmer_reference(line, left_most, right_most, bin_mask, t1, t2);
                    break;
                }
                line_count += 1;

                lineStream.str(lineStream.str().substr(pos + 1));
                lineStream.seekp(0,ios_base::end);
            }
        }

        gzclose(gzfile);
    }
    else{
        cout << "False format!\n";
        exit(EXIT_FAILURE);  
    }
    if(round != 0){
        line_count = 0;
        while (getline(ref_fa, line))
        {
            if (line_count % round == 0)
            {
                ids = split_find(line, " ");
                name = ids[0].substr(1);
            }
            else if(line_count % round == 1 && strcmp(name.c_str(), chr) == 0){
                search_kmer_reference(line, left_most, right_most, bin_mask, t1, t2);
                break;
            }
            line_count += 1;
        }        
    }
    cout << "Searching kmer in reference done!\n";
    ref_fa.close();


    cout << "Searching kmer in query file!\n"; // TODO: multi line fa, fq, gz
    mutex qry_m[t2];
    for(const auto& qry_file : qry_files){
        if(qry_file[strlen(qry_file.c_str()) - 1] == 'z'){
            if(qry_file[strlen(qry_file.c_str()) - 4] == 'a')
                round = 2;
            else if(qry_file[strlen(qry_file.c_str()) - 4] == 'q')
                round = 4;
            else{
                cout << "False format!\n";
                exit(EXIT_FAILURE);  
            }

            line_count = 0;
            gzFile gzfile = gzopen(qry_file.c_str(), "rb");
            if(gzfile == nullptr){
                std::cerr << "Failed to open GZ file: " << qry_file << std::endl;
                return 1;
            }

            const int BUFFER_SIZE = 1024;
            char buffer[BUFFER_SIZE];

            ostringstream lineStream;

            while(true){
                int bytesRead = gzread(gzfile, buffer, BUFFER_SIZE);

                if(bytesRead <= 0){
                    if(gzeof(gzfile)){
                        // handle the last line without '\n'
                        line = lineStream.str();
                        if(!line.empty()){
                            line_count += 1;
                        }
                    }
                    else{
                        std::cerr << "Error reading GZ file: " << qry_file << std::endl;
                    }

                    break;
                }

                lineStream.write(buffer,bytesRead);

                while(true){
                    size_t pos = lineStream.str().find('\n');
                    if(pos == std::string::npos){
                        // no '\n'
                        break;
                    }

                    line = lineStream.str().substr(0, pos);
                    line_count += 1;

                    lineStream.str(lineStream.str().substr(pos + 1));
                    lineStream.seekp(0,ios_base::end);
                }
            }

            gzclose(gzfile);
        }
        else if(qry_file[strlen(qry_file.c_str()) - 1] == 'a' || qry_file[strlen(qry_file.c_str()) - 1] == 'q'){
            if(qry_file[strlen(qry_file.c_str()) - 1] == 'a')
                round = 2;
            else if(qry_file[strlen(qry_file.c_str()) - 1] == 'q')
                round = 4;
            else{
                cout << "False format!\n";
                exit(EXIT_FAILURE);  
            }          

            ifstream qry(qry_file);
            line_count = 0;
            while(getline(qry,line)){
                line_count += 1;
            }
            qry.close();              
        }
        else{
            cout << "False format!\n";
            exit(EXIT_FAILURE); 
        }
        cout << line_count << endl;

        uint64_t chunk_line_num =  ((line_count / round + t2 - 1) / t2) * round;
        for (int i = 0; i < t2; i++)
        {
            fpool.submit(search_kmer_reads, qry_file.c_str(), chunk_line_num, i, round, bin_mask, z_value, &file[0], &qry_m[0], t2);
        }        
    }
    cout << "Searching kmer in query file done!\n";


    fpool.shutdown();
    cout << "Shutting down!\n";

    for (int i = 0; i < t2; i++)
    {
        fclose(file[i]);
    }
    string cmd = "cat " + out_path + "/temp*.fasta > " + out_path + "/read_" + chr + "_" + to_string(gap_start) + "_" + to_string(gap_end) + "_" + to_string(significance) + ".fasta;" + "rm " + out_path + "/temp*.fasta";
    system(cmd.c_str());
#if 0
    string cmd,cmd2;
	if(type){
		cmd="cat " + out_path +"/*_r.pos > "+ out_path + "/ref.pos";
		cmd2="rm "+out_path + "/*_r.pos";
	}
	else{
		cmd="cat " + out_path +"/*_q.pos > "+ out_path + "/query.pos";
		cmd2="rm "+out_path + "/*_q.pos";
	}
	system(cmd.c_str());
	system(cmd2.c_str());
#endif

    return 1;
}
