#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>

// the "raw" prefix table
using td_t = std::set<uint32_t>;
using db_t = std::map<uint64_t, td_t>;

#define DEFAULT_KMER_LEN 20

#define SANITY_INTERVAL 1000

#define LMAT_ERR(e) { \
  std::cerr << "ERROR! file: " << __FILE__ << " line: " << __LINE__  \
            << "; error: " << e << std::endl; \
  exit(9); \
}

#define FWRITE_STATUS(b) { \
  if (b != 1) { \
    std::cerr << "ERROR; file: "<<__FILE__<<"; line: "<<__LINE__<<"; err_msg: fwrite failed\n"; \
    exit(9); \
  }}   

#define FREAD_STATUS(b) { \
  if (b != 1) { \
    std::cerr << "ERROR; file: "<<__FILE__<<"; line: "<<__LINE__<<"; err_msg: fread failed\n"; \
    exit(9); \
  }}   

// reads a list of filenames from 'fn'
void get_filenames(const std::string fn, std::vector<std::string> &filenames); 

// Formats number with commas for easier comprehe;nsion
std::string commify(size_t n);

std::string get_binary_filename(int thread_id, uint64_t prefix, std::string base_output_dir); 

namespace metag {
    /*
     * static functions that aren't particular to any specific class
     *
     */
    
void print_invocation(int argc, char **argv);

    // a string of 64 1s after every SANITY_FREQUENCY kmers
#define SANITY_FREQUENCY 1000
    
#define CHECK_HELP {for (int j=0; j<argc; j++) { \
if (strcmp(argv[j], "-h") == 0) { \
usage(argv); \
} \
} \
exit(0); }
    
    
    class Utils {
        public :
        
        //! the input file should be a kmer data file;
        //! counts the number of TaxNode IDs required
        //! for all KmerNodes
        //  static uint64_t countTaxNodeStorage(const char *fn, GenomeIdToTaxId &lookup);
        
        //! returns the number of bytes in the file;
        //! on entry file offset can be anywhere;
        //! on return, file offset is zero.
        static size_t getFileSize(FILE *fp);
        
        //! returns the number of bytes in the file;
        static size_t getFileSize(const char *fn);
        
        //! reads the contents of the file into a buffer;
        //! caller is responsible for freeing the returned buffer;
        //! on return, fp will be at offset zero.
        static char * readFile(FILE *fp);
        
        //! reads the contents of the file into a buffer;
        //! caller is responsible for freeing the returned buffer
        static char * readFile(const char *fn);
        
        //! opens a file for binary reading
        static  FILE * openReadFile(const char *fn);
        
        //! opens a file for binary writing
        static  FILE * openWriteFile(const char *fn);
        
        //! closes a file
        static void closeFile(FILE *fp);
        
        //! opens a file; this is used internally: end users
        //! should probably call openReadFile or openWriteFile instead
        static FILE * openFile(const char *fp, const char *mode);
        
        //! returns, in 'out,' all tax IDs associated with kmers
        //! in 'fn,' which is a kmer data filename
        //  static void getTids(const char *fn, GenomeIdToTaxId &genome_to_taxid, std::set<uint32_t> &out);
        
        static uint64_t countLines(const char *fn, bool fn_is_file_containing_list_of_files=true);
        
        static void skipHeaderLines(std::ifstream &in);
        
        static void readTaxHisto_ascii(uint64_t &kmer, uint16_t &genome_count, uint16_t &tid_count, std::map<uint32_t, uint16_t> &tids, std::ifstream &in);
        
        static void readTaxHisto_bin(uint64_t &kmer, uint16_t &genome_count, uint16_t &tid_count, std::map<uint32_t, uint16_t> &tids, std::ifstream &in);
        
        static void readTaxHisto_bin(uint64_t &kmer, uint16_t &genome_count, uint16_t &tid_count, std::map<uint32_t, uint16_t> &tids, FILE *in);
        static uint64_t encode(std::string &s); 
        
        static void decode(uint64_t j, int kmer_len, std::string &mer);
        
        static void showBits(uint64_t kmer, uint64_t num_to_print);
        
    };

    
} // namespace metag 

#endif
