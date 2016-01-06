#include "KmerFileMetaData.hpp"
#include "Utils.hpp"
#include <iostream>
#include <cassert>
#include <fstream>

using namespace std;
using namespace metag;

void KmerFileMetaData::write(const char *fn) {
  FILE *fp = Utils::openReadFile(fn);
  write(fp);
  Utils::closeFile(fp);
}

void KmerFileMetaData::write(FILE *fp) {
  //cout << "starting KmerFileMetaData::write\n";
  assert(fwrite(&m_data_start, 4, 1, fp) == 1);

  assert(fwrite(&m_kmer_count, 8, 1, fp) == 1);

  uint64_t sanity = ~0;
  assert(fwrite(&sanity, sizeof(uint64_t), 1, fp) == 1);

  assert(fwrite(&m_version, sizeof(uint32_t), 1, fp) == 1);

  char c = m_has_locations? 'Y' : 'N';
  assert(fwrite(&c, 1, 1, fp) == 1);

  assert(fwrite(&m_kmer_len, sizeof(uint32_t), 1, fp) == 1);

  //uint32_t mark = ftell(fp);
  //assert(mark == m_data_start);
}

void KmerFileMetaData::write(ostream &out) {
  cout << "K:             " << kmerLength() << endl;
  cout << "kmer count:    " << size() << endl;
  cout << "version:       " << version() << endl;
  cout << "data start:    " << dataStart() << endl;
  cout << "has locations? " << hasLocations() << endl;
}

void KmerFileMetaData::read(FILE *fp) {
  uint64_t sanity = ~0;
  uint64_t test;

  assert(fread(&m_data_start, sizeof(uint32_t), 1, fp) == 1);
  assert(fread(&m_kmer_count, sizeof(uint64_t), 1, fp) == 1);
  assert(fread(&test, sizeof(uint64_t), 1, fp) == 1);

  if (test != sanity) {
    cerr << __FILE__<<" "<<__LINE__ ;
    cerr << " kmer data file is invalid; should have read 64 1s, but didn't\n";
    cerr << "from file: " << test << " should be: " << sanity << endl;
    exit(-1);
  }

  //todo: check that the version is valid??
  assert(fread(&m_version, sizeof(uint32_t), 1, fp) == 1);

  char c;
  assert(fread(&c, 1, 1, fp) == 1);
  if (! (c == 'Y' || c == 'N')) {
    cerr << "invalid location flag; should be 'Y' or 'N', but got " << c << endl;
    exit(-1);
  }

  m_has_locations = false;
  if (c == 'Y') {
    m_has_locations = true;
  }


  #if USE_GENOME_LOCATIONS == 1
    if (m_has_locations == false) {
      cerr << "your data file was generated without genome locations,\n"
           << "but the metag library was compiled using USE_GENOME_LOCATIONS=1\n";
      cerr << "please recompile the library or use a different kmer data file\n";
      exit(0);
    }
  #else
    if (m_has_locations == 1) {
      cerr << "your data file was generated with genome locations,\n"
           << "but the metag library was compiled using USE_GENOME_LOCATIONS=0\n";
      cerr << "please recompile the library or use a different kmer data file\n";
      exit(0);
    }  
  #endif

  assert(fread(&m_kmer_len, sizeof(uint32_t), 1, fp) == 1);
  assert(ftell(fp) == m_data_start);
  m_was_initialized = true;
}

void KmerFileMetaData::read(const char *fn) {
  FILE *fp = Utils::openReadFile(fn);
  read(fp);
  Utils::closeFile(fp);
}

