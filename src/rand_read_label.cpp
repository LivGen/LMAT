#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include <cassert>
#include <vector>
#include <cmath>
#include <list>
#include <stdlib.h>
#include <omp.h>
#include <random>
#include <functional>
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"
#include <gzstream.h>

#define MMAP_SIZE 0
#define TID_T uint32_t

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;

//using namespace kencode_ns;
using namespace metag;

#include "rkmer.hpp"

bool verbose=false;

bool add_root_on_kmer_drop = true;
bool gPERMISSIVE_MATCH = false;
map<TID_T,string> gRank_table;


size_t perm_bytes_allocd;

typedef pair<TID_T,float> ufpair_t;
typedef pair<TID_T,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<int16_t,tax_data_t> label_info_t;
typedef map<TID_T,TID_T> hmap_t;
typedef map<TID_T,float> ufmap_t;
typedef map<uint16_t, ufmap_t> u_ufmap_t;

vector <int> read_len_vec;
vector <int> read_len_avgs;

#define _USE_KPATH_IDS 0


my_map tid_rank_map;
id_convback_map_t conv_map;


// this determines which style of tid reduction we want, default is
// the "comphrehensive reduction; requires the tid to rank mapping
// table.  -w option enables strain to species only reduction.
bool tid_map_is_strain_species = false;

typedef std::pair<TaxNode<TID_T>*,TID_T> tncpair_t;

struct SimpleCmp {
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
      return a.second > b.second;
   }
};

struct CmpDepth {
   CmpDepth(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
      const int adepth = (*_imap.find(a.first)).second;
      const int bdepth = (*_imap.find(b.first)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};

void genRandRead(string& read_buff, int beg, int end) {
   const unsigned rl = read_buff.length();
   const int range = (end-beg)+1;
   const int gc_draw = (rand() % range) + beg;
   const float gc_pcnt = gc_draw/100.0;
   const unsigned num_gc = static_cast<unsigned>(gc_pcnt*rl);
   for(unsigned i = 0; i < num_gc; ++i) {
      const int coin_flip = rand() % 100;
      const char base = coin_flip < 50 ? 'g' : 'c';
      read_buff[i] = base;
   }
   for(unsigned i = num_gc; i < rl; ++i) {
      const int coin_flip = rand() % 100;
      const char base = coin_flip < 50 ? 'a' : 't';
      read_buff[i] = base;
   }
   random_shuffle(read_buff.begin(),read_buff.end());
   //cout<<"debug: "<<beg<<" "<<end<<" "<<gc_pcnt<<" "<<num_gc<<" "<<read_buff<<endl;
}

void
fill_in_labels(const TaxTree<TID_T>& taxtree, vector<TID_T>& row, const tax_data_t& taxids, const hmap_t& idx2taxid, const hmap_t& taxid2idx) {
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const TID_T col_tax_val = (*mtch).second;
      const TaxTree<TID_T>::const_iterator col_tax_val_it = taxtree.find(col_tax_val); 
      if(col_tax_val_it == taxtree.end()) {
         cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
         continue;
      }
      TaxNode<TID_T>* col_tax_val_node = (*col_tax_val_it).second;
      // see if this id's genome count needs to be update
      //  based on its children   - can skip if its a leaf
      if( !col_tax_val_node->isLeaf() && row[tax_idx] == 0 ) {
         tax_data_t::const_iterator it = taxids.begin();
         const tax_data_t::const_iterator is = taxids.end();
         unsigned has_cnt_sum = 0;
         for(; it != is; ++it) {   
            const TID_T tax_id = (*it).first;
            const TaxTree<TID_T>::const_iterator tax_val_it = taxtree.find(tax_id); 
            if(tax_val_it == taxtree.end()) {
               cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
               continue;
            }
            TaxNode<TID_T>* tax_val_node = (*tax_val_it).second;
            // test if tax_id is a descendant of tax_val
            // save max value
            const hmap_t::const_iterator idx_mtch = taxid2idx.find(tax_id);
            assert(idx_mtch != taxid2idx.end());
            const TID_T tax_val_idx = (*idx_mtch).second;
            
            if( col_tax_val != tax_id && row[tax_val_idx] > 0 && tax_val_node->isAncestor(col_tax_val)) {
               // want to capture the descendant tax nodes closest to col_tax_val_node 
               has_cnt_sum = 1;
               break;
            }
         }
         if( has_cnt_sum > 0 ) {
            row[tax_idx] = has_cnt_sum;
            if(verbose) cout<<"Final Tax Node Val "<<col_tax_val<<" final transfer cnt for taxid="<<row[tax_idx]<<endl;
         }
      }
   }
   if(verbose) cout<<"fill_row: "; 
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const TID_T tax_val = (*mtch).second;
      if(verbose) cout<<"["<<tax_val<<","<<row[tax_idx]<<"]";
   }
   if(verbose) cout<<endl;
}

struct TCmp { 
   TCmp(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
	if( fabs(a.second- b.second) < 0.001 ) {
		const int adepth = (*_imap.find(a.first)).second;
		const int bdepth = (*_imap.find(b.first)).second;
		return adepth < bdepth;
        } 
	return a.second < b.second; }
   const hmap_t& _imap;
};

struct ScoreOptions {
   ScoreOptions(hmap_t& imap) : _equal_kmer_vote(false), _strict_kmer_match(false), _prn_all(false), _imap(imap), _diff_thresh(1.0), _diff_thresh2(3.0) {}
   bool _equal_kmer_vote;
   bool _strict_kmer_match;
   bool _prn_all;
   hmap_t& _imap;
   float _diff_thresh, _diff_thresh2;
   u_ufmap_t _rand_hits; 
   bool _comp_rand_hits;
};

typedef map<TID_T,vector<float> > match_map_t;
typedef map<TID_T,vector<int> > match_cnt_map_t;
void
construct_labels(const TaxTree<TID_T>& tax_tree, 
                match_map_t& max_match, match_cnt_map_t& match_cnt, const ScoreOptions& sopt,
                int cand_kmer_cnt, const map<TID_T,int>& cnt_tids, unsigned gcbucket, unsigned num_gcbuckets) {
   //cout<<"check: "<<cnt_tids.size()<<endl;
   map<TID_T,int>::const_iterator it = cnt_tids.begin();
   const map<TID_T,int>::const_iterator is = cnt_tids.end();
   for(; it != is; ++it) {
      const TID_T taxid = (*it).first;
      const int found_genome_cnt = (*it).second;
      const float label_prob = (float)found_genome_cnt / (float) cand_kmer_cnt;
      if(verbose) cout<<"chk: "<<taxid<<" "<<found_genome_cnt<<" "<<label_prob<<endl;
      if( found_genome_cnt > cand_kmer_cnt ) cout<<" How does thishappen: "<<found_genome_cnt<<" "<<cand_kmer_cnt<<" "<<taxid<<endl;
      assert(label_prob <= 1);
      if( max_match.find(taxid) == max_match.end() ) {
         //max_match.insert(make_pair(taxid,label_prob));
         //match_cnt.insert(make_pair(taxid,1));
         vector<float> buff(num_gcbuckets,0);
         buff[gcbucket] = label_prob;
         max_match.insert(make_pair(taxid,buff));
         vector<int> buff2(num_gcbuckets,0);
         buff2[gcbucket] = 1;
         match_cnt.insert(make_pair(taxid,buff2));
      } else {
         const float curr_max = max_match[taxid][gcbucket];
         max_match[taxid][gcbucket] = (curr_max < label_prob) ? label_prob : curr_max;
         match_cnt[taxid][gcbucket] += 1;
      }
   }
}

#define ENCODE(t, c, k) \
switch (c) { \
case 'a': case 'A': t = 0; break; \
case 'c': case 'C': t = 1; break; \
case 'g': case 'G': t =2; break; \
case 't': case 'T': t = 3; break; \
default: k = 0; continue; \
}

#if 0
/*

   Given a sequence of nucleotide characters,

   break it into canonical k-mers in one pass.

   Nucleotides are encoded with two bits in

   the k-mer. Any k-mers with ambiguous characters

   are skipped.

   str:  DNA sequence (read)

   slen: DNA sequence length in nucleotides

   klen: k-mer length in nucleotides

*/

static

int retrieve_kmer_labels(INDEXDB<DBTID_T>* table, const char* str, const int slen, const kmer_t klen, 
                          const hmap_t& dmap, const TaxTree<TID_T>& tax_tree, uint16_t max_count, map<TID_T,int>& cnt_tids) 
   {
    int j; /* position of last nucleotide in sequence */
    int k = 0; /* count of contiguous valid characters */
    int highbits = (klen-1)*2; /* shift needed to reach highest k-mer bits */
    kmer_t mask = ((kmer_t)1 << klen*2)-1; /* bits covering encoded k-mer */
    kmer_t forward = 0; /* forward k-mer */
    kmer_t reverse = 0; /* reverse k-mer */
    kmer_t kmer_id; /* canonical k-mer */
    set<kmer_t> no_dups;
    int valid_kmers=0;
    set<TID_T> track_tid;
    for (j = 0; j < slen; j++) {
        register int t;
        ENCODE(t, str[j], k);
        forward = ((forward << 2) | t) & mask;
        reverse = ((kmer_t)(t^3) << highbits) | (reverse >> 2);
        if (++k >= (signed)klen) {
           valid_kmers++;
           kmer_id = (forward < reverse) ? forward : reverse;
            /* zero based position of forward k-mer is (j-klen+1) */
            /* zero based position of reverse k-mer is (slen-j-1) */
           const int pos = j-klen+1;
           /* kmer_lookup(kmer); do k-mer lookup here... */
           if(verbose) cout<<"debug valid: "<<j<<" "<<pos<<endl;
           if( no_dups.find(kmer_id) != no_dups.end() ) continue;
           no_dups.insert(kmer_id);

           TaxNodeStat<DBTID_T> *h = new TaxNodeStat<DBTID_T>(*table);
           if(verbose) cout<<"lookup kmer at posit:"<<j<<endl;
#if (TID_SIZE == 16)
           h->begin(kmer_id, tid_rank_map,  max_count, tid_map_is_strain_species, &conv_map);
#else
           h->begin(kmer_id, tid_rank_map,  max_count, tid_map_is_strain_species);
#endif
           unsigned dcnt = 0, mtch = 0;
           list<TID_T> obs_tids;
           set<TID_T> once;
           while( h->next() ) {
              const TID_T tid = h->taxid();
              if( tid == 20999999 ) continue;
              if( once.find(tid) != once.end() ) continue;

              once.insert(tid); 
              if( track_tid.find(tid) == track_tid.end()) {
                 track_tid.insert(tid);
                 cnt_tids[tid] = 1;
              } else {
                 cnt_tids[tid] += 1;
                 if( cnt_tids[tid] > valid_kmers ) {
                    cout<<"happens here: "<<tid<<" "<<cnt_tids[tid]<<" "<<valid_kmers<<" "<<kmer_id<<endl;
                    TaxNodeStat<DBTID_T> *hdebug = new TaxNodeStat<DBTID_T>(*table);
                    hdebug->begin(kmer_id, tid_rank_map,  max_count, tid_map_is_strain_species);
                    while(hdebug->next()) {
                        const TID_T tiddebug = hdebug->taxid();
                        cout<<"strange: "<<tiddebug<<endl;
                    }
                    exit(0);
                 }
              } 
              //cout<<"comeone: "<<pos<<" "<<tid<<" "<<cnt_tids[tid];
              obs_tids.push_back(tid);
           }
           vector<TID_T> obs_tids_vec(obs_tids.size());
           list<TID_T>::const_iterator it = obs_tids.begin();
           const list<TID_T>::const_iterator is = obs_tids.end();
           for(unsigned i =0; it != is; ++it, ++i) {
              obs_tids_vec[i] = *it;
           }
           CmpDepth1 cd(dmap);
           sort(obs_tids_vec.begin(),obs_tids_vec.end(),cd);
           int last_depth = -1;
           set<TID_T> sanity_check;
           for(unsigned i = 0; i < obs_tids_vec.size(); ++i) {
             const TID_T tid=obs_tids_vec[i];
             const int depth = (*dmap.find(tid)).second;
             if( depth == 0 ) {
               break;
             }
             if( last_depth == depth || last_depth == -1 ) {
                  vector<TID_T> path;
                  TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
                  tax_tree_tmp.getPathToRoot(tid,path);
                  for(unsigned p = 0; p < path.size(); ++p) {
                     const TID_T ptid = path[p];
                  // There seems to be alot of wasted visits that need to be fixed here
                     if( sanity_check.find(ptid) != sanity_check.end() ) {
                        cout<<"Not what I was expecting (should see just once no?): "<<pos<<" "<<ptid<<" "<<tid<<endl;
                        continue;
                     }
                     if( once.find(ptid) != once.end() ) {
                        //cout<<"We already saw this one?: "<<pos<<" "<<ptid<<" "<<tid<<endl;
                        break;
                     }
                     once.insert(ptid);
                     if( track_tid.find(ptid) == track_tid.end()) {
                       track_tid.insert(ptid);
                       cnt_tids[ptid] = 1;
                     } else {
                       cnt_tids[ptid] += 1;
                       if( cnt_tids[tid] > valid_kmers ) {
                          cout<<"happens here2: "<<tid<<" "<<cnt_tids[tid]<<" "<<valid_kmers<<endl;
                          exit(0);
                       }
                     }
                  }
             } else {
               break;
             }
          }
          if(dcnt > 0 && verbose ) cout<<" end k-mer lookup"<<endl;
          if(mtch == 0 && verbose ) cout<<" no k-mer matches "<<endl;
          delete h;
       }
   }
   return valid_kmers;
}
#endif

void proc_line(const TaxTree<TID_T>& tax_tree, 
               const string& line,
               int k_size, int max_count, INDEXDB<DBTID_T> *table, const ScoreOptions& sopt, 
               match_map_t& max_match, match_cnt_map_t& cnt_match, unsigned gcbucket, unsigned num_gcbuckets ) {
     const unsigned ri_len = line.length();
     if( (signed)ri_len < k_size ) {
         return;
     } else {
        map<TID_T,int> cnt_tids;
        vector<label_info_t> label_vec(ri_len-k_size+1,make_pair(-1,tax_data_t()));
        list<TID_T> taxid_lst;
        hmap_t tax2idx, idx2tax;
        const pair<int,int> res = retrieve_kmer_labels(table, line.c_str(), ri_len, k_size,label_vec,taxid_lst,tax2idx,idx2tax, sopt._imap, tax_tree, max_count);
        const int valid_kmers = res.first;
        if( valid_kmers > 0 ) {
           for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
              tax_data_t::const_iterator it = label_vec[pos].second.begin();
              const tax_data_t::const_iterator is = label_vec[pos].second.end();
              for(; it != is; ++it) {
                 const TID_T tax_id = (*it).first;
                 if( cnt_tids.find(tax_id) == cnt_tids.end() ) {
                    cnt_tids.insert(make_pair(tax_id,1));
                 } else {
                    cnt_tids[tax_id] += 1;
                 } 
              } 
            }
            construct_labels(tax_tree,max_match,cnt_match,sopt,valid_kmers, cnt_tids, gcbucket, num_gcbuckets); 
        }
     }
  }


void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <input db file (list)> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:restore persistent db] [-h:max_tid_count\n";
  cout << "[-r <rank/tid-map-file>] [-h <tid-cutoff>] [-w <with-strain-species-map> (affects -r option)]\n";
  cout << "note: -r makes -k unnecessary\n"; 
}


int main(int argc, char* argv[]) 
{
   std::srand ( unsigned ( std::time(0) ) );
   char c = '\0';
   int k_size=-1;
   unsigned n_threads = 0;

#include <random>
#include <functional>

/*
May want to consider this snippet, HOWEVER, since we are using 80 threads to s 
std::uniform_int_distribution<int> dice_distribution(1, 2);
std::mt19937 random_number_engine; // pseudorandom number generator
auto dice_roller = std::bind(dice_distribution, random_number_engine);
int random_roll = dice_roller();  
*/
   bool restore = true;

   float threshold = 0.0, min_score = 0.0;
   int min_kmer = 35;

   string rank_map_file, rank_ids, kmer_db_fn, ofname, ofbase, tax_tree_fn, tax_tree_options, depth_file, rand_hits_file, rank_table_file, id_bit_conv_fn;
   hmap_t imap;	
   ScoreOptions sopt(imap);

   size_t mmap_size = 0;
   uint16_t max_count = ~0;
   bool prn_read = true;
   unsigned num_reads = 0, read_len = 0;

   while ((c = getopt(argc, argv, "u:ah:n:j:b:ye:w:mpk:c:v:k:i:d:l:t:s:r:o:x:f:g:z:q:")) != -1) {
      switch(c) {
      case 'f':
        id_bit_conv_fn = optarg;
        break;
      case 'e':
         depth_file = optarg;
         break;
      case 'j':
         min_kmer = atoi(optarg);
         break;
      case 'u':
         rank_ids = optarg;
         break;
      case 'x':
         min_score = atof(optarg);
         break;
      case 'a':
         prn_read=false;
         break;
      case 'w':
         rank_map_file = optarg;
	      break;
      case 'h' :
        max_count = atoi(optarg);
        break;
      case 's':
         mmap_size = atoi(optarg);
         mmap_size = mmap_size * (1<<30);
         cout << "Input heap size: " << mmap_size << endl;
         break;
      case 'n':
         rand_hits_file=optarg;
         break;
      case 'y':
         verbose = true;
         break;
      case 'p':
         sopt._prn_all = true;
         break;   
      case 'r':
	      rank_table_file = optarg;
        break;
      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
        break;
      case 'v':
         threshold = atof(optarg);
         break;
      case 'c':
         tax_tree_fn = optarg;
	      break;
      case 'k':
         k_size = atoi(optarg);
         break;
      case 'g':
         num_reads = atoi(optarg);
         break;
      case 'i':
         read_len = atoi(optarg);
         break;
      case 'd':
         kmer_db_fn = optarg;
         break;
      case 'o':
         ofbase = optarg;
         break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }
   cout<<"Total reads to evaluate: "<<num_reads*n_threads<<endl;
   if (ofbase == "") cout << "ofbase\n";
   if (n_threads == 0) cout << "n_threads\n";
   if (kmer_db_fn == "") cout << "kmer_db_fn\n";

   if( rank_map_file.size() > 0) {
      ifstream ifs1(rank_map_file.c_str());
      TID_T tid;
      string rank;
      while(ifs1>>tid>>rank) {
         gRank_table.insert(make_pair(tid,rank));
      }
   }


#if TID_SIZE == 16
   if (id_bit_conv_fn == "") {
     usage(argv[0]);
     return -1;
   }
#endif

   if (ofbase == "" || n_threads == 0 || kmer_db_fn == "" ) {
     cout<<ofbase<<" "<<n_threads<<" "<<kmer_db_fn<<" "<<depth_file<<endl; 
     usage(argv[0]);
     return -1;

   }
   if (!restore && k_size == -1)  {
     cout << "missing kmer size!\n";
     usage(argv[0]);
     return -1;

   }
   if (id_bit_conv_fn.length() > 0) {
     cout << "Loading map file: "<<id_bit_conv_fn<<endl;
     FILE * tfp = fopen(id_bit_conv_fn.c_str(), "r");
     if( !tfp ) {
         cerr<<"Unable to read: "<<id_bit_conv_fn<<endl; 
         return -1;
     }
     uint32_t src;
     uint16_t dest;

     while (fscanf(tfp,"%d%hd", &src, &dest) > 0) {

       conv_map[dest] = src;
     }
     fclose(tfp);
   }


   cout << "Start kmer DB load\n";
   INDEXDB<DBTID_T> *taxtable;


#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<INDEXDB<DBTID_T>*, std::size_t> ret = mfile.find<INDEXDB<DBTID_T>>("KmerDB");

     taxtable = ret.first;
     k_size = taxtable->get_kmer_length();
     cout << "k size:  " << k_size  <<   endl ;
     taxtable->conv_ptrs();
#else


#if WITH_PJMALLOC == 1

   if (restore) {
      
     perm(&taxtable, sizeof(taxtable));
     if( mopen(kmer_db_fn.c_str(), "r", mmap_size) != 0 ) {
         cerr<<"Unable to open ["<<kmer_db_fn<<"]"<<endl;
         return -1;
     }
     if (k_size < 1)
       k_size = taxtable->get_kmer_length();

     cout << "num kmers: " << taxtable->size() << " - " << k_size  <<   endl ;


   } else 
    
#endif


#if (USE_SORTED_DB == 0)
     taxtable = new INDEXDB<DBTID_T>;

     ifstream qifs(kmer_db_fn.c_str());
     if( !qifs ) {
       cerr<<"Unable to open: "<<kmer_db_fn<<endl;
       return -1;


     }



     string fname;
     
     while(qifs>>fname) {
       cout<<"register file: "<<fname<<endl;
       taxtable->registerFile(fname.c_str());
     }
     taxtable->ingest();
#endif


   //}

#endif
   //cout << "End kmer DB load\n";
   //cout << "DB size is " << table->size() << endl;
   if( k_size <= 0 ) {
      cerr<<"Unable to read database, k_size="<<k_size<<endl;
      return -1;
   }
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */

   ofstream ofs;

   if (rank_table_file.length() > 0) {
     if(max_count == ~0) {
       cout << "Need to set -h <tid-cutoff> to use rank file map!\n";
     } else {
       FILE * rmfp = fopen(rank_table_file.c_str(), "r");
     
       uint32_t src, dest;
       
       while (fscanf(rmfp,"%d%d", &src, &dest) > 0) {
	      tid_rank_map[src] = dest;
       }
       
       fclose(rmfp);
       
     }
	 
   }
   cout<<"Read taxonomy tree: "<<tax_tree_fn<<endl;
   TaxTree<TID_T> tax_tree(tax_tree_fn.c_str());
   cout<<"Read taxonomy depth: "<<depth_file<<endl;
   ifstream ifs1(depth_file.c_str());
   if(!ifs1) {
   cerr<<"unable to open: "<<depth_file<<endl;
   return -1;
   }
   TID_T taxid,depth;
   while(ifs1>>taxid>>depth) {
   sopt._imap[taxid] = depth;
   }


   StopWatch clock;
   clock.start();
   vector<match_map_t> max_match_all(n_threads); 
   vector<match_cnt_map_t> match_cnt_all(n_threads); 
   const int min_gc=0;
   const int num_bins = 10;
   float width = 100.0/(float)num_bins;
   vector<pair<int,int> > gc_range(num_bins);
   float lval=min_gc;
   const unsigned num_gcbuckets = gc_range.size();
   for(unsigned i = 0; i < gc_range.size(); ++i) {
      const int i_lval = static_cast<float>(lval);
      const int rval = static_cast<float>(lval+width-1);
      gc_range[i] = make_pair(i_lval,rval);
      cout<<"gc check "<<i<<" "<<gc_range[i].first<<" "<<gc_range[i].second<<endl;
      lval += width;
   }
   

#pragma omp parallel shared(max_count,read_len, num_reads, gc_range,k_size, ofbase, taxtable, tax_tree, sopt, max_match_all,match_cnt_all,n_threads)  
   {
   const int thread = omp_get_thread_num();
   string read_buff(read_len,'\0');
   match_map_t& max_match = max_match_all[thread]; 
   match_cnt_map_t& match_cnt = match_cnt_all[thread]; 
   for(unsigned i = 0; i < num_reads; ++i) {
      const int gc_bucket = i % num_gcbuckets;
      const int beg = gc_range[gc_bucket].first; 
      const int end = gc_range[gc_bucket].second; 
      genRandRead(read_buff,beg,end);
	   proc_line(tax_tree, read_buff, k_size, max_count, taxtable, sopt, max_match, match_cnt, gc_bucket, num_gcbuckets );
      if( i % 10000 == 0 ) cout<<"progress "<<thread<<" "<<i<<endl;
   }
}
   cout<<"Merge phase"<<endl;
   match_cnt_map_t  merge_count;
   match_map_t  merge_score;
   for(unsigned thread = 0; thread < max_match_all.size();  ++thread) {
      match_map_t& gt = max_match_all[thread];
      match_map_t::const_iterator it = gt.begin();
      const match_map_t::const_iterator is = gt.end();
      for(; it != is; ++it) {
         TID_T tid = (*it).first;
         const vector<float>& score = (*it).second;
         if( merge_score.find(tid) == merge_score.end() ) {
            merge_score.insert( make_pair(tid,score));
         } else {
            const vector<float>& curr_score = merge_score[tid];
            for(unsigned ti = 0; ti < curr_score.size();++ti) {
               merge_score[tid][ti] = (score[ti] > curr_score[ti]) ? score[ti] : curr_score[ti];
            }
         }
      }
      match_cnt_map_t& gtcnt = match_cnt_all[thread];
      match_cnt_map_t::const_iterator it1 = gtcnt.begin();
      const match_cnt_map_t::const_iterator is1 = gtcnt.end();
      for(; it1 != is1; ++it1) {
         const TID_T tid = (*it1).first;
         const vector<int>& cnt = (*it1).second;
         if( merge_count.find(tid) == merge_count.end() ) {
            merge_count.insert( make_pair(tid,cnt));
         } else {
            for(unsigned ti = 0; ti < cnt.size();++ti) {
               merge_count[tid][ti] += cnt[ti];
            }
         }
      }
   }
   ostringstream ostrm;
   ostrm<<ofbase<<".rand_lst";
   ofstream sum_ofs(ostrm.str().c_str());
   if( !sum_ofs ) {
      cout<<"Could not open for writing "<<ostrm.str()<<endl;
      return -1;
   }
   match_map_t::const_iterator it = merge_score.begin();
   match_map_t::const_iterator is = merge_score.end();
   for(unsigned i = 0; it != is; ++it, ++i) {
      const TID_T tid = (*it).first;
      const vector<float>& max_score = (*it).second;
      assert( merge_count.find( tid ) != merge_count.end());
      const vector<int>& cnt = (*merge_count.find(tid)).second;
      sum_ofs<<tid;
      for(unsigned val = 0; val < num_gcbuckets; ++val ) {
         sum_ofs<<" "<<max_score[val]<<" "<<cnt[val];
      }
      sum_ofs<<endl;
   }
   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
