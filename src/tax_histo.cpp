#include <unistd.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <set>
#include "all_headers.hpp"
#include "TaxTree.hpp"
#include "metag_typedefs.hpp"
#include <version.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::set;
using std::map;
using std::fstream;
using std::ofstream;
using std::ostream;

using namespace metag;

template<class tid_T>
static bool
isHuman(tid_T id) {
   bool res=false;
   switch(id) {
   case 10000348:
   case 10000349:
   case 10000350:
   case 10000351:
   case 10000352:
   case 10000353:
   case 10000354:
   case 10000355:
   case 10000356:
   case 10000357:
   case 10000358:
   case 10000359:
   case 10000360:
   case 10000361:
   case 10000362:
   case 10000363:
   case 10000364:
   case 10000365:
   case 10000366:
   case 10000367:
   case 10000368:
   case 10000369:
   case 10000370:
   case 10000371:
   case 10000372:
   case 10000373:
   case 9606:
   case 63221: //neanderthal:
      res=true;
      break;
   }
   return res; 
}

template<class tid_T>
static bool
hasHuman(const set<tid_T>& tax_ids, TaxTree<tid_T>& tax_tree) {
   bool human=false;
   for (typename set<tid_T>::const_iterator t = tax_ids.begin(); t != tax_ids.end(); t++) {
      TaxNode<tid_T> *tn = tax_tree[*t];
      if(tn) {
         if(isHuman(tn->id())) {
            //cout<<"found human skip "<<tn->id()<<endl;
            human=true;
            break;
         }
      } else {
        cout << "failed to find node in TaxTree for tid: " << *t << endl;
      }
   }
   return human;
}

template<class tid_T>
void doit(string &taxtree_fn, string &kmer_db_fn, string &outfile, string &ranks_fn, size_t quit_early); 

void
usage(const string& prog_name, const string& optstr);

int main(int argc, char* argv[]) {

  cout << "invocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;

  char c = '\0';
  bool prn_help = false;
  string kmer_db_fn, outfile;
  string taxtree_fn, ranks_fn;
  size_t quit_early = 0;
  int count = 0;
  int which = 0;

  const string opt_string="q:o:t:d:h:r:f:V";
  while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
    switch (c) {
    case 'f' :
      ++count;
      which = atoi(optarg);
      break;
    case 'r' :
      ranks_fn = optarg;
      break;
    case 'q':
      quit_early = atoi(optarg);
      break;
    case 'o':
      ++count;
      outfile = optarg;
      break;
    case 'd':
      ++count;
      kmer_db_fn = optarg;
      break;
    case 't':
      ++count;
      taxtree_fn = optarg;
      break;
    case 'h':
      prn_help = true;
      break;
    case 'V':
      cout << "LMAT version " << LMAT_VERSION  << "\n";
	exit(0);
    default:
      cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      prn_help = true;
      break;
    }
  }
  if ( prn_help || count != 4 || !(which == 16 || which == 32)) {
    usage(argv[0],opt_string);
    exit(0);
  }

  if (which == 16) {
    //doit<uint16_t>(taxtree_fn, kmer_db_fn, outfile, ranks_fn, quit_early);
    cout<<"not working yet with this hack"<<endl;
    return -1;
  } else {
    doit<uint32_t>(taxtree_fn, kmer_db_fn, outfile, ranks_fn, quit_early);
  }
}

template<class tid_T>
void doit(string &taxtree_fn, string &kmer_db_fn, string &outfile, string &ranks_fn, size_t quit_early) {
  cout << "info: starting tax tree load from filename: " << taxtree_fn << endl;
  TaxTree<tid_T> tax_tree(taxtree_fn.c_str());
  if (ranks_fn.size()) {
    tax_tree.setRanks(ranks_fn.c_str());
  }
  cout << "info: tax tree size: " << tax_tree.size() << endl;

  cout << "info: start kmer DB load\n";
  FILE *fp = Utils::openReadFile(kmer_db_fn.c_str());
  KmerFileMetaData metadata;
  metadata.read(fp);




  uint64_t kmer_count = metadata.size();

  uint64_t test, sanity = ~0;
  std::set<tid_T> tax_ids;

  cerr << "opening for writing: " << outfile.c_str() << endl;
  FILE *out_bin = fopen(outfile.c_str(), "wb");
  assert(out_bin);

  //write metadata
  metadata.setVersion(TAX_HISTO_VERSION);
  metadata.write(out_bin);

  set<tid_T> tids_that_were_already_processed;
  KmerNode<tid_T> w;
  StopWatch c2;
  c2.start();
  uint64_t j;
  uint64_t kmer;
  size_t count = 0;

  set<tid_T> bad_tid;
  string mer;
  uint64_t tid_ct = 0;
  cout << "starting; kmer count: " << kmer_count << endl;

  map<tid_T, set<tid_T> > species_test;

  map<int, int> all_species_test;

  uint64_t total_tid = 0;

  size_t singletons = 0;

  int ignore_kmer_cnt = 0; 

  for (j=0; j<kmer_count; j++) {

    if (quit_early && j == quit_early) {
      cout << "quit_early: " << j << endl;
      break;
    }

    w.read(fp);

    const set<tid_T> tax_ids = w.getTaxIDs();

    //allen99 quick hack to remove human k-mers
    bool doWrite=true;

    __gnu_cxx::hash_map<tid_T, set<tid_T> > tid_set;

    /*

    if( hasHuman( tax_ids, tax_tree ) ) {
      doWrite=false;
      ++ignore_kmer_cnt;
    } else {
    */

    tax_tree.getLcaMap(tax_ids, tid_set);




      if (tid_set.size() == 0) {
	cout << "\nfrom tax_histo_new_fmt: WARNING: tid_set is empty; no entry will be written for kmer " << w.getKmer() << endl;
	cout << "    this is for kmer #" << j+1 << " of " << kmer_count << endl;
	cout << "    entries in tax_id set: ";
	for (typename set<tid_T>::iterator t = tax_ids.begin(); t != tax_ids.end(); t++) {
	  cout << *t << " ";
	}
	cout <<  endl << endl;
	doWrite = false;

      } else  {

	++count;
 
 
	kmer= w.getKmer();
    

	if(doWrite)assert(fwrite(&kmer,8, 1, out_bin) == 1);  //write kmer
	uint16_t sz = tid_set.size();
	if(doWrite) assert(fwrite(&sz,2, 1, out_bin) == 1); //write tid count
	tid_ct += tid_set.size();



	if (sz == 1)
	  singletons++;

	//write the tax IDs

	for (typename __gnu_cxx::hash_map<tid_T, set<tid_T> >::const_iterator t = tid_set.begin(); t != tid_set.end(); t++) {	  
	  tid_T tid = t->first;

	  if(doWrite) assert(fwrite(&tid,4, 1, out_bin) == 1);

	  total_tid++;
	}
      
	if ((count) % TAX_HISTO_SANITY_COUNT ==0) {


	  assert(fwrite(&sanity, 8, 1, out_bin) == 1);
	}


      }


      if ((j+1) % KMER_SANITY_COUNT == 0) {
	assert(fread(&test, 8, 1, fp) == 1);
	assert(test == sanity);
      }
      

      //    }  // else for has human

  }

  cout << "total taxids: " << total_tid << "\nrem kmer cnt: "<<ignore_kmer_cnt<<endl;
  cout << "singletons: " << singletons << endl;

  
  fseek(out_bin, 0, SEEK_SET);


  
  metadata.setSize(count);
  metadata.write(out_bin);

  fclose(fp);
  fclose(out_bin);

  double tm = c2.stop();
  cout << endl << "total annotate time: " << tm << endl;
  cout << "num mapping kmers processed: " << count << endl;
  cout << "kmers per second (not counting startup time): " << (double)kmer_count/tm << endl;
}


void
usage(const string& prog_name, const string& optstr) {
  cout<<"Usage: \n"
      "  -f 16|32      - determines if 16 or 32 bits are used for tax IDs [required]\n"
      "  -o <string>   - output fn           [required]\n"
      "  -d <string>   - kmer data filename  [required]  \n"
      "  -t <string>   - tax tree data file  [required]\n"
      "  -q <int>      - quit after processing q kmers [optional]\n"
      "  -r <string>   - file containing taxid to rank mapping [optional; used to collect species statistics]\n"
      "  -h            - print help and exit\n"
      "\n";
}
