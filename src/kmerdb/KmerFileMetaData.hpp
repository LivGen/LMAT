#ifndef __KMERFILEMETADATA__
#define __KMERFILEMETADATA__

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <cstdlib>

namespace metag {

/**
container class for metadata at the beginning of binary kmer data files
output from jellylist_bin.  Includes read, write, and print methods.
*/

class KmerFileMetaData {

public:

  //! default ctor
  KmerFileMetaData() : m_has_locations(false), m_was_initialized(false)  {}

  void setWasInitialized() {
    m_was_initialized = true;
  }

  //! returns the length of the kmers in the file
  uint32_t kmerLength() const {
    checkmate();
    return m_kmer_len;
  }

  void setKmerLength(int s) {
    m_kmer_len = s;
  }

  //! returns the number of kmers in the file
  uint64_t size() {
    checkmate();
    return m_kmer_count;
  }

  //! sets the number of kmers in the file
  void setSize(uint64_t s) {
    m_kmer_count = s;
  }

  //! 
  void setVersion(uint32_t v) {
    m_version = v;
  }

  //! returns the file's version number
  uint32_t version() {
    checkmate();
    return m_version;
  }

  //! returns true if the data file contains genome location information
  bool hasLocations() {
    checkmate();
    return m_has_locations;
  }

  //! returns the offset in the file where data for the first kmer begins
  uint32_t dataStart() {
    checkmate();
    return m_data_start;
  }

  //! prints class variables
  void write(std::ostream &out = std::cout);

  //! set the number of kmers; this is needed by
  //! partition_kmer_data_file.cpp
  void setCount(uint64_t c) {
    m_kmer_count = c;
  }

  void setTidCount(uint32_t ct) {
    m_data_start = ct;
  }

  uint32_t tidCount() {
    return m_data_start;
  }

  void setDefaultDataStart() {
    m_data_start = 29;
  }

  //! opens file, reads meta-data, closes file; exits with 
  //! error message if the kmer data file is found to be corrupt.
  //! The caller is responsible for deleting the KmerMetaData structure.
  void read(const char *fn); 

  //! reads meta-data; exits with error message if the 
  //!kmer data file is found to be corrupt; on return,
  //! fp will be at the offset where data for the 1st kmer begins.
  //! The caller is responsible for deleting the KmerMetaData structure.
  void read(FILE *fp);

  //! writes meta-data to file
  void write(FILE *fp);

  //! opens file, writes metadata; on return, fp will be
  //! at the offset where data for the firest kmer should
  //! be written
  void write (const char *fn);

private :
  uint32_t m_kmer_len;
  uint64_t m_kmer_count;
  uint32_t m_version;
  bool m_has_locations;   //true, if the file contains genome location
  uint32_t m_data_start;  //file offset where data for the first kmer begins;
                          //when used with tax_histo files, stores the number
                          //of tax IDs in the file
                         

  bool m_was_initialized;

  void checkmate() const {
    if (! m_was_initialized) {
      std::cerr << "you must call readMetaData before invoking this method\n";
      exit(-1);
    }
  }
  
};

}
#endif
