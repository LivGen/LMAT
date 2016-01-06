#include "Utils.hpp"
#include <iostream>
#include <cassert>
#include <string>
#include <fstream>
#include <set>
#include "all_headers.hpp"

using namespace std;
using namespace metag;


size_t Utils::getFileSize(FILE *fp) {
  fseek(fp, 0, SEEK_END);
  long kmer_file_size = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  return (size_t)kmer_file_size;
}

size_t Utils::getFileSize(const char *fn) {
  FILE *fp = openReadFile(fn);
  size_t sz = getFileSize(fp);
  closeFile(fp);
  return sz;
}


FILE * Utils::openFile(const char *fn, const char *mode) {
  FILE *fp = fopen(fn, mode);
  if (fp == NULL) {
    cerr << "failed to open " << fn << " for mode \"" << mode << "\" at " << __FILE__ << " " << __LINE__ << endl;
    exit(-1);
  }
  return fp;

}
   
FILE * Utils::openReadFile(const char *fn) {
  return openFile(fn, "r");
}
        
FILE * Utils::openWriteFile(const char *fn) {
  return openFile(fn, "w");
}

void Utils::closeFile(FILE *fp) {
  if (fclose(fp)){
    cerr << "fclose failed\n";
    exit(-1);
  }
}

char * Utils::readFile(FILE *fp) {
  size_t sz = getFileSize(fp);
  char *retval = new char[sz];
  assert(fread(retval, 1, sz, fp) == sz);
  fseek(fp, 0, SEEK_SET);
  return retval;
}

char * Utils::readFile(const char *fn) {
  FILE *fp = openReadFile(fn);
  char *retval = readFile(fp);
  closeFile(fp);
  return retval;
}

uint64_t countLinesInOnefile(const char *fn) {
  cout << "Utils::countLinesInOnefile(): " << fn << endl;
  FILE *fp = fopen(fn, "r");
  if (!fp) {
    cerr << "Utils::countLinesInOnefile(), failed to open " << fn << " for reading\n";
    exit(1);
  }
  uint64_t ct = 0;
  int c;
  while(true) {
    c = fgetc(fp);
    if (c == EOF) break;
    if (c == '\n') {
      ++ct;
    }  
  }
  fclose(fp);
  return ct;
}

uint64_t Utils::countLines(const char *fn, bool fn_list) {
  if (!fn_list) {
    return countLinesInOnefile(fn);
  }
  string line;
  ifstream in(fn);
  assert(in);
  uint64_t ct = 0;
  while (true) {
    line = "";
    getline(in,line);
    if (line.size() < 2) break;
    ct += countLinesInOnefile(line.c_str());
  }
  return ct;
}

void Utils::skipHeaderLines(std::ifstream &in) {
  int c;
  while (true) {
    c = in.peek();
    string line;
    if (c == '#') {
      getline(in, line);
    } else {
      break;
    }
  }
}

void Utils::readTaxHisto_bin(uint64_t &kmer, uint16_t &genome_count, uint16_t &tid_count, std::map<uint32_t, uint16_t> &tids, std::ifstream &in) {
  tids.clear();
  uint16_t tuple_count;
  in.read((char*)&kmer, 8);
  in.read((char*)&genome_count, 2);
  in.read((char*)&tid_count, 2);
  in.read((char*)&tuple_count, 2);
  uint32_t tid;
  uint16_t present;
  for (size_t j=0; j<tuple_count; j++) {
    in.read((char*)&tid, 4);
    in.read((char*)&present, 2);
    tids[tid] = present;
  }
}

void Utils::readTaxHisto_bin(uint64_t &kmer, uint16_t &genome_count, uint16_t &tid_count, std::map<uint32_t, uint16_t> &tids, FILE *in) {
  tids.clear();
  uint16_t tuple_count;
  int x = fread(&kmer, 8, 1, in);
  if (x != 1) cout << ">>>>>>> " << ftell(in) << endl;
  assert(fread(&genome_count, 2, 1, in) == 1);
  assert(fread(&tid_count, 2, 1, in) == 1);
  assert(fread(&tuple_count, 2, 1, in) == 1);
  uint32_t tid;
  uint16_t present;
  for (size_t j=0; j<tuple_count; j++) {
    assert(fread(&tid, 4, 1, in) == 1);
    assert(fread(&present, 2, 1, in) == 1);
    tids[tid] = present;
  }
}

void Utils::readTaxHisto_ascii(uint64_t &kmer, uint16_t &genome_count, uint16_t &tid_count, std::map<uint32_t, uint16_t> &tids, std::ifstream &in) {
  tids.clear();
  uint16_t tuple_count;
  in >> kmer;
  in >> genome_count;
  in >> tid_count;
  in >> tuple_count;
  uint32_t tid;
  uint16_t present;
  for (size_t j=0; j<tuple_count; j++) {
    in >> tid;
    in >> present;
    tids[tid] = present;
  }
}

uint64_t Utils::encode(string &s) {
  uint64_t out = 0;
  int len = s.size();
  for (int j=0; j<len; j++) {
    switch (s[j]) {
      case 'A' :
      case 'a' : out = out << 2;
                 break;
      case 'C' :
      case 'c' : out = out << 2;
                 out = out | 1;
                 break;
      case 'G' :
      case 'g' : out = out << 2;
                 out = out | 2;
                 break;
      case 'T' :
      case 't' : out = out << 2;
                 out = out | 3;
                 break;
      default : cout << "ERROR!\n"; exit(1);
    }
  }
  return out;
}

void Utils::decode(uint64_t j, int kmer_len, string &mer) {
  static char *x = 0;
  if (!x) x = new char[kmer_len+1];
  static uint64_t mask = 3, tmp;
  static char c;
  x[kmer_len] = '\0';
  for (int k=0; k<kmer_len; k++) {
    tmp = j & mask;
    switch (tmp) {
      case 0 : c = 'a'; break;
      case 1 : c = 'c'; break;
      case 2 : c = 'g'; break;
      case 3 : c = 't'; break;
      default : cerr << "ERROR!!!!!!!!!\n"; exit(1);
    }
    x[kmer_len-k-1] = c;
    
    j = j >> 2;
  }
  mer = x;
}

void Utils::showBits(uint64_t kmer, uint64_t num_to_print) {
  vector<bool> b;
  uint64_t mask = 1;
  for (int j=0; j<num_to_print; j++) {
    int i = kmer & mask;
    b.push_back(i);
    kmer = kmer >> 1;
  }

  int k = 0;
  for (int j = b.size()-1; j>=0; j--, ++k) {
    cout << b[j] << " ";
    if ((k+1) % 4 == 0) cout << " ";
  }
  cout << endl;

}

