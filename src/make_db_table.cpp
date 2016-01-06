/**
Create a mmap kmer db file of the desired data structure
*/
#ifndef USE_SORTED_DB
#define USE_SORTED_DB 1
#endif

#include "all_headers.hpp"
#include "SortedDb.hpp"

#include "version.h"

#define DATA_DIR "../../src/kmerdb/examples/tests/data/"

#include <fstream>
#include <omp.h>
#include <getopt.h>

using namespace std;
using namespace metag;

#ifdef USE_BOOST
#include <boost/interprocess/managed_mapped_file.hpp>
#endif

#define ILLU_TAXID  32630


bool use_tax_histo_format = true;


//typedef pair<size_t, int> kmer_info_t; 

void usage() {
  cout << "LMAT version " << LMAT_VERSION  << "\n";
  cout<<"Usage: \n"
    //      "  -f 16|32  - determines if 16 or 32 bit tax IDS are used [required]\n"
    "  -i <fn>   - input file -or- filename of file that contains a listing of the input binary files\n"
    
    "              output from tax_histo_fast_limited [required]\n"
    "  -l        - is file a list?        [default=no, optional]\n"       
    "  -o <fn>   - output filename                    [required]\n"
    "  -k <int>  - k-mer length used in input file    [required]\n"
    "  -h        - if given, assumes input files are from kmerPrefixCounter, instead of tax_histo [optional]\n"
    "  -s <int>  - size to reserve for the memory-mapped file [optional; default is -s 500, or 500G]\n"
    "  -q <int>  - stop at N k-mers (per input file) - helpful for troubleshooting \n"
    "  -g <int>  - taxid list cutoff \n"
    "  -m <fn>   - pruning taxid table file  \n"    
    "  -f <fn>   - use 32 to 16-it id conversion using <fn>  \n"    
    "  -j <fin>  - file containing human kmers \n"
    "  -c <int>  - count of human kmers \n"
    "  -w        - pruing uses strain-to-species mapping\n" 
    "  -u        - file for illumina adaptor k-mers\n"
    "  -V        - print version and exit\n"
    "guidance for setting -s: from our paper, our full reference DB required 619G;\n"
    "this was for a fasta file that was ~19G\n";
}


size_t  single_f_count(const char *fname)
{
  FILE *fp = fopen(fname, "r");

  if (!fp)
    {
      cout << "Input file " << fname << " not found!\n";
      exit(1);

    }

      KmerFileMetaData metadata;
      metadata.read(fp);

      if (use_tax_histo_format)
	assert(metadata.version() == TAX_HISTO_VERSION);


      return metadata.size();
      fclose(fp);
}


size_t get_kmer_count_preproc(const char *fname, bool list)
{
  if (list) {
    ifstream input_st(fname);

    string line;
    size_t kmer_count = 0;

    while (input_st >> line) {

      kmer_count += single_f_count(line.c_str()); 
    }
    return kmer_count;
  }
    else {
      return single_f_count(fname);
   }    

}



int main(int argc, char *argv[]) {

  string inputfn, outputfn;

  int c, n_threads;

  inputfn = DATA_DIR "test.fa.int.bin.no_loc.12";
  outputfn = DATA_DIR "test.fa.mmap";
  

  bool restore = false;

  // default database size in GB
  size_t mmap_size = 255;

  bool list = false;

  int kmer_len = 0;

  size_t hash_size = 1;
  size_t storage_size = 1;

  unsigned long long stopper = 0;
  n_threads = -1;

  int count = 0;  

  int tid_cut = 0;

  size_t xtra_kmers = 0;

  string species_map_fn, id_bit_conv_fn, human_kmer_fn, illu_kmer_fn;

  bool strainspecies = false;
  //  while ((c = getopt(argc, argv, "t:g:q:k:i:o:s:t:h:r l a ")) != -1) {

  cout << "invocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;



 
 while ((c = getopt(argc, argv, "g:q:k:i:o:s: l h m:f:wj:c:u:V")) != -1) {
    switch(c) {
    case 'j':
      human_kmer_fn=optarg;
      break;
    case 'u':
      illu_kmer_fn=optarg;
      break;
    case 'w':
      strainspecies = true;
      break;
    case 'f':
      id_bit_conv_fn = optarg;
      break;
    case 'q':
      stopper = strtoull(optarg, NULL, 10);
      break;
      /*    case 'h':
      hash_size = strtoull(optarg, NULL, 10);
      break; */
    case 'k':
      ++count;
      kmer_len = atoi(optarg);
      break;
    case 'l':
      list = true;
      break;
      /*    case 'r':
      restore = true;
      break;*/
    case 'i':
      ++count;
      inputfn = optarg;
      break;
    case 'c':
      xtra_kmers = strtoull(optarg, NULL, 10);
      break;
    case 'o':
      ++count;
      outputfn = optarg;
      break;
      //    case 't':
      //      n_threads = atoi(optarg);
      // omp_set_num_threads(n_threads);
      //    break;
    case 's':
      mmap_size = atoi(optarg);
      break;
    case 'h':
      use_tax_histo_format = false;
      break;
    case 'g':
      tid_cut = atoi(optarg);
      break;
    case 'm':
      species_map_fn = optarg;
      break;
    case 'V':
      cout << "LMAT version " << LMAT_VERSION  << "\n" ; 
	exit(0);
    default:
      usage();
      exit(1);
    }
  }

  if (count != 3)    {
    usage();
    exit(1);
  }
  
  if (strainspecies && species_map_fn.length() == 0) {
    cout << "missing pruning map file!\n";
    exit(1);
    
  }
  

  // setup kmerdb
  SortedDb<DBTID_T> *ttable;
  


  mmap_size = mmap_size * (1<<30);

  cout << "size requested: " << mmap_size << endl;


  hash_size = get_kmer_count_preproc(inputfn.c_str(), list);

  cout << hash_size << " k-mers found in input files.\n";

  // calculate how much space is left over for tid list storage 
  // for now keep 1 extra GB for keeping 
  //(TODO: narrow that down to necessary amount)


  size_t space = ((TT_BLOCK_COUNT * 8) + (hash_size * 8));

  storage_size = (size_t)((mmap_size - space) * (double).95) - (size_t)(xtra_kmers*8);

 
  cout << "storage requested: " << storage_size << endl;

  bitreduce_map_t  br_map;

  my_map species_map;


  if (id_bit_conv_fn.length() > 0) { 
     FILE * tfp = fopen(id_bit_conv_fn.c_str(), "r");
     
     uint32_t src; 
     uint16_t dest;
     
     while (fscanf(tfp,"%d%hd", &src, &dest) > 0) {

       //       assert(dest < 37369);       
       br_map[src] = dest;

     }

     fclose(tfp);
  }
  
  FILE *human_fp = NULL;

  if (human_kmer_fn.length() > 0) {

    
    human_fp = fopen(human_kmer_fn.c_str()  ,"r");
    
    if (human_fp)
      cout << "opened human k-kmer file at " << human_kmer_fn << "\n" ;
	
    
  }

  if (! human_fp)
    cout << "No human k-mer file.\n";


  FILE *illu_fp = NULL;

  if (illu_kmer_fn.length() > 0) 
    
    illu_fp = fopen(illu_kmer_fn.c_str()  ,"r");
    

  if (! illu_fp)
    cout << "No Illumina k-mer file.\n";

   
   if (  (tid_cut > 0) && species_map_fn.length() > 0) {
     
     FILE * smfp = fopen(species_map_fn.c_str(), "r");
     
     uint32_t src, dest;
     
     while (fscanf(smfp,"%d%d", &src, &dest) > 0) {
       
       species_map[src] = dest;
     }

     fclose(smfp);
   }
   
   


#if USE_BOOST == 1
  bip::managed_mapped_file mfile(bip::create_only, (const char*) outputfn.c_str(
										) , mmap_size );  
  ttable = mfile.construct<SortedDb<uint32_t> >("KmerDB")(hash_size, storage_size, mfile);
  // ,n_threads);
  

#else


  perm(&ttable, sizeof(ttable));
  
 

  if (restore) 
    mopen(outputfn.c_str(), "r+", mmap_size);
 
  else {
    mopen(outputfn.c_str(), "w+", mmap_size);
    
    if (kmer_len < 16)
      hash_size=(1<<(kmer_len*2));

    ttable = PERM_NEW(SortedDb<DBTID_T>)(hash_size+xtra_kmers, storage_size);

  }

#endif

  
  ttable->set_kmer_length(kmer_len);

  

  StopWatch clock;
  clock.start();




  if (list) {

    std::vector<string> input_files;
  
    ifstream ifs(inputfn.c_str());
    
    string line;
  
    while (ifs>>line) {
      input_files.push_back(line);
    }

    int i;

    if (n_threads < 2)
      {
	for (i = 0 ; i < input_files.size(); i++) {

	StopWatch clock2;

	clock2.start();
	
	bitreduce_map_t *p_map = NULL;
	if (br_map.size() > 0)
	  p_map = & br_map;
	ttable->add_data(input_files[i].c_str() , stopper, use_tax_histo_format, p_map,  species_map, tid_cut, strainspecies, human_fp, illu_fp, ILLU_TAXID);
	cout << "elapsed time: " << clock2.stop() << endl;
	} // end for

      }
    /*
    else


      
#pragma omp parallel shared(ttable, input_files, stopper) private(i)

    {

      for (i = 0 ; i < input_files.size()/ n_threads ; i++) {

	StopWatch clock2;

	clock2.start();
      
	ttable->add_data(input_files[i*n_threads+omp_get_thread_num()].c_str() , stopper, n_threads, omp_get_thread_num());
	cout << "elapsed time: " << clock2.stop() << endl;
      } // end for
    } // end pragma
    */    
  } 

  else {
    cout << "input fn: " << inputfn << endl;
    bitreduce_map_t *p_map = NULL;
    if (br_map.size() > 0)
      p_map = & br_map;
    ttable->add_data(inputfn.c_str(), stopper, use_tax_histo_format, p_map,  species_map, tid_cut, strainspecies, human_fp, illu_fp, ILLU_TAXID); 
  }
  




  cout << "KmerDB load time: " << clock.stop() << endl;


#if USE_BOOST != 1

  mclose();
#endif

  return(0);
}
