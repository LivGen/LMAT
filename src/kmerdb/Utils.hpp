#ifndef __UTILS_HPP__
#define __UTILS_HPP__

//#include <cstdio>
//#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <stdint.h>

namespace metag {
/*
 * static methods that aren't particular to any specific class
 *
 */

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

}
#endif
