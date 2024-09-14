#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <mutex>
#include <bitset>
#include <algorithm>
#include <unordered_set>
#include <zlib.h>
#include "thread_pool.hpp"
#include "uthash.h"
#include "kmerpos.hpp"
#include "stringtools.hpp"

using namespace std;

// vector<phmap::flat_hash_set<uint64_t>> kmer; //存储kmer的hash>
mutex m;
int k_size;          // k-mer length

typedef struct my_struct {
    uint64_t key;                    /* key */
    //char name;
    UT_hash_handle hh;         /* makes this structure hashable */
} hash_node;

class hash_set{
public:
    hash_node *dict=NULL;
    void insert(uint64_t *key){
        hash_node *s;
        HASH_FIND_INT(dict, key, s);
        if (s == NULL)
        {
            s = (struct my_struct*)malloc(sizeof *s);
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
    uint64_t size() {
        return HASH_COUNT(dict);
    }
};

vector<hash_set> kmer;


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

int creat_dict(const char *k_path, uint64_t chunk_line_num, int idx)
{
    ifstream k_mer(k_path);
    // unordered_set<uint64_t> dict;
    hash_set dict;
    uint64_t k_value = 0;
    // printf("%s",s);
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
    t = 1;
    k_mer.seekg(0,k_mer.end);
    uint64_t total_len=k_mer.tellg();
    k_mer.seekg(0,k_mer.beg);
    string s;
    getline(k_mer, s);
    uint64_t line_len = s.size()+1;
    k_mer.close();
    uint64_t chunk_line_num = (total_len / line_len + t - 1) / t;

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
    uint64_t sbin, rbin;
    uint32_t pos;
    string nameid;
    int flag;
    int is_ukmer_strand(int t, FILE *fp)
    {        
        for(int i = 0;i < t;i ++){
            if (kmer[i].find(&sbin) != NULL)
            {
                if(!flag){
                    flag = 1;
                    fprintf(fp, "@%s\n", nameid.c_str());
                }   
                fprintf(fp, "%-8x\t%0*lx\t+\n", pos, (k_size / 2) + 1, sbin);
            }
            else if (kmer[i].find(&rbin) != NULL)
            {
                if(!flag){
                    flag = 1;
                    fprintf(fp, "@%s\n", nameid.c_str());
                }   
                fprintf(fp, "%-8x\t%0*lx\t-\n", pos, (k_size / 2) + 1, rbin);
            }
        }

        
        #if 0
        for(int i = 0;i < t;i ++){
            if (kmer[i].find(&sbin) != NULL)
            {
                fprintf(fp, "%-8x\t%0*lx\t+\n", pos, (k_size / 2) + 1, sbin);
            }
            else if (kmer[i].find(&rbin) != NULL)
            {
                fprintf(fp, "%-8x\t%0*lx\t-\n", pos, (k_size / 2) + 1, rbin);
            }
        }
        if (kmer[0].find(sbin) != kmer[0].end())
        {
            if(!flag){
                flag = 1;
                fprintf(fp, "@%s\n", nameid.c_str());
            }   
            fprintf(fp, "%-8x\t%0*lx\n", pos, (k_size / 2) + 1, sbin);
        }
        else if (kmer[0].find(rbin) != kmer[0].end())
        {
            if(!flag){
                flag = 1;
                fprintf(fp, "@%s\n", nameid.c_str());
            }
            fprintf(fp, "%-8x\t%0*lx\n", pos, (k_size / 2) + 1, rbin);
        }
        #endif
        return 0;
    }
    int is_ukmer_not_strand(int t, FILE *fp)
    {
        for(int i = 0;i < t;i ++){
            if (kmer[i].find(&sbin) != NULL)
            {
                fprintf(fp, "%-8x\t%0*lx\n", pos, (k_size / 2) + 1, sbin);
            }
        }
        return 0;
    }
};

int search_kmer_strand(string line, int t1, string name, long bin_mask, FILE *file[], bool mask[], int n)
{
    // Find a idle file
    FILE *fp;
    int index = -1;
    m.lock();
    for (int i = 0; i < n; i++)
    {
        if (mask[i] == 0)
        {
            index = i;
            fp = file[i];
            mask[i] = 1;
            break;
        }
    }
    if(index == -1){
        cout << "Wild pointer!\n";
        exit(EXIT_FAILURE);        
    }
    m.unlock();

    KMER k;
    // k.flag = 0;
    vector<string> id = split_find(name, " ");
    k.nameid = id[0];
    k.flag = 0;
    // fprintf(fp, "@%s\n", id[0].c_str());
    uint64_t base;
    uint32_t len = line.length();
    string qkmer = line.substr(0, k_size);
    k.sbin = tobin(&qkmer);
    reverse(qkmer.begin(),qkmer.end());
    k.rbin = (tobin(&qkmer) ^ bin_mask);
    // cout<<k.rbin<<" "<<k.sbin<<"\n";   
    for (uint32_t i = k_size; i < len; i++)
    {
        k.pos = i - (k_size-1);

        k.is_ukmer_strand(t1, fp);

        base = ctoi(line[i]);
        k.sbin = (k.sbin << 2 | base) & bin_mask;
        k.rbin = (k.rbin >> 2 | ((3ull - base) << ((k_size - 1) * 2)));
    }
    k.pos++;
    k.is_ukmer_strand(t1, fp);

    m.lock();
    mask[index] = 0;
    m.unlock();
    return 1;
}

int search_kmer_not_strand(string line, int t1, string name, long bin_mask, FILE *file[], bool mask[], int n)
{
    // Find a idle file
    FILE *fp;
    int index;
    m.lock();
    for (int i = 0; i < n; i++)
    {
        if (mask[i] == 0)
        {
            index = i;
            fp = file[i];
            mask[i] = 1;
            break;
        }
    }
    if(index == -1){
        cout << "Wild pointer!\n";
        exit(EXIT_FAILURE);        
    }
    m.unlock();

    KMER k;
    vector<string> id = split_find(name, " ");
    fprintf(fp, "@%s\n", id[0].c_str());
    uint64_t base;
    uint32_t len = line.length();
    string qkmer = line.substr(0, k_size);
    k.sbin = tobin(&qkmer);
    // cout<<k.rbin<<" "<<k.sbin<<"\n";
    for (uint32_t i = k_size; i < len; i++)
    {
        k.pos = i - (k_size-1);

        k.is_ukmer_not_strand(t1, fp);

        base = ctoi(line[i]);
        k.sbin = (k.sbin << 2 | base) & bin_mask;
    }
    k.pos++;
    k.is_ukmer_not_strand(t1, fp);

    m.lock();
    mask[index] = 0;
    m.unlock();
    return 1;
}

int build_pos(vector<string> target_files, string out_path,int k, bool strand, int t1, int t2)
{
    cout << "Number of threads: " << t2 << "\n"
         << "Reading the kmer to dicts\n";
    ThreadPool fpool(t2);
    fpool.init();
    string line, name;
    int round = 0;
    k_size = k;
    FILE *file[t2];

    /*
        Creat mutiple files for mutiple threads;
        File pointers are store in a FILE * array.
    */
    char outname[15];
    strcpy(outname, "/kmerAA.pos");
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
    bool mask[t2] = {0};            // In order to search a idle file to write
    
    long bin_mask = 0x3;
    for(int i = 0;i < k_size - 1;i++){
        bin_mask = (bin_mask << 2) | bin_mask;
    }

    cout << "Searching kmer in query file!\n";
    for(const auto& target_file : target_files){
        if(target_file[strlen(target_file.c_str()) - 1] == 'z'){
            if(target_file[strlen(target_file.c_str()) - 4] == 'a')
                round = 2;
            else if(target_file[strlen(target_file.c_str()) - 4] == 'q')
                round = 4;
            else{
                cout << "False format!\n";
                exit(EXIT_FAILURE);  
            }

            uint64_t count = 0;
            gzFile gzfile = gzopen(target_file.c_str(), "rb");
            if(gzfile == nullptr){
                std::cerr << "Failed to open GZ file: " << target_file << std::endl;
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
                            if(count % round == 0)
                                name = line.substr(1, line.size());
                            else if(count % round == 1){
                                if(strand){
                                    fpool.submit(search_kmer_strand, line, t1, name, bin_mask, &file[0], &mask[0], t2);
                                }
                                else{
                                    fpool.submit(search_kmer_not_strand, line, t1, name, bin_mask, &file[0], &mask[0], t2);
                                }
                            }
                            count += 1;
                        }
                    }
                    else{
                        std::cerr << "Error reading GZ file: " << target_file << std::endl;
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

                    if(count % round == 0)
                        name = line.substr(1, line.size());
                    else if(count % round == 1){
                        if(strand){
                            fpool.submit(search_kmer_strand, line, t1, name, bin_mask, &file[0], &mask[0], t2);
                        }
                        else{
                            fpool.submit(search_kmer_not_strand, line, t1, name, bin_mask, &file[0], &mask[0], t2);
                        }
                    }
                    count += 1;

                    lineStream.str(lineStream.str().substr(pos + 1));
                    lineStream.seekp(0,ios_base::end);
                }
            }

            gzclose(gzfile);
        }
        else if(target_file[strlen(target_file.c_str()) - 1] == 'a' || target_file[strlen(target_file.c_str()) - 1] == 'q'){
            if(target_file[strlen(target_file.c_str()) - 1] == 'a')
                round = 2;
            else if(target_file[strlen(target_file.c_str()) - 1] == 'q')
                round = 4;
            else{
                cout << "False format!\n";
                exit(EXIT_FAILURE);              
            }

            ifstream fl(target_file);
            uint64_t count = 0;
            while (getline(fl, line))
            {
                if (count % round == 0)
                {
                    name = line.substr(1, line.size());
                }
                else if(count % round == 1){
                    if(strand){
                        fpool.submit(search_kmer_strand, line, t1, name, bin_mask, &file[0], &mask[0], t2);
                    }
                    else{
                        fpool.submit(search_kmer_not_strand, line, t1, name, bin_mask, &file[0], &mask[0], t2);
                    }
                }
                count += 1;
            }
            fl.close();
        }
        else{
            cout << "False format!\n";
            exit(EXIT_FAILURE);  
        }        
    }


    cout << "Searching kmer done!\n";
    fpool.shutdown();
    for (int i = 0; i < t2; i++)
    {
        fclose(file[i]);
    }
    
    cout << "Shutting down!\n";
    string cmd = "cat " + out_path + "/kmer*.pos > " + out_path + "/ref.pos;" + "rm " + out_path + "/kmer*.pos";
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
