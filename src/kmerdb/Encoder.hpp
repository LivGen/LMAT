#ifndef __ENCODER_HPP__
#define __ENCODER_HPP__

#include <string>
#include <cassert>
#include <stdint.h>
#include <cstdlib>
#include <iostream>

namespace metag {

class Encoder {
public :
  Encoder(std::string &seq, int k) : m_seq(seq), m_k(k), m_idx(0), m_cur(0) {
    m_mask = ~0;
    int n = 64 - 2*(m_k);
    m_mask = m_mask >>n;
  } 

  static uint64_t encode(std::string &kmer) {
    uint64_t out = 0;
    for (size_t j=0; j<kmer.size(); j++) {
      switch (tolower(kmer[j])) {
        case 'a' : out = out << 2; 
                   break;
        case 'c' : out = out << 2;
                   out = out | 1;
                   break;
        case 'g' : out = out << 2;
                   out = out | 2;
                   break;
        case 't' : out = out << 2;
                   out = out | 3;
                   break;
        default : std::cout << "ERROR!\n"; exit(1);
      }
    }
    return out;
  }

  static void decode(uint64_t kmer, int kmer_len, std::string &mer) {
    mer.resize(kmer_len);
    uint64_t mask = 3, tmp;
    char c;
    for (int k=0; k<kmer_len; k++) {
      tmp = kmer & mask;
      switch (tmp) {
        case 0 : c = 'a'; break;
        case 1 : c = 'c'; break;
        case 2 : c = 'g'; break;
        case 3 : c = 't'; break;
        default : std::cerr << "ERROR!!!!!!!!!\n"; exit(1);
      }
      mer[kmer_len-k-1] = c;
      kmer = kmer >> 2;
    }
  }

  //! returns the reverse complement
  static uint64_t rc(uint64_t kmer, int kmer_len) {
    uint64_t rev = 0, work, tmp;
    uint64_t mask = 3;
    for (int j=0; j<kmer_len; j++) {
      work = kmer & mask;
      switch (work)  {
        case 0 : tmp = 3; break;
        case 1 : tmp = 2; break;
        case 2 : tmp = 1; break;
        case 3 : tmp = 0; break;
        default : assert(false);
      }
      rev = rev << 2;
      rev = rev | tmp;
      kmer = kmer >> 2;
    }
    return rev;
  }

  static bool rc(const std::string &fwd, std::string &rc) {
    rc.resize(fwd.size());
    char c;
    for (size_t j=0; j<fwd.size(); j++) {
      c = tolower(fwd[j]);
      switch(c) {
        case 'a' : c = 't'; break;
        case 't' : c = 'a'; break;
        case 'c' : c = 'g'; break;
        case 'g' : c = 'c'; break;
        default : assert(false);
      }
      rc[rc.size()-1-j] = c;
    }
    return true;
  }

  bool next(uint64_t &output) {
   //case 3: no kmers remain
      if (m_idx == m_seq.size()) {
          return false;
      }
    //case 1: this executes the first time next() is called
    if (m_idx == 0) {
      bool status = start_over();
      if (!status) {
        return false;
      }
      output = m_cur;
      return true;
    }

    //case 2: this executes when we encounter a degenerate bp
    char x = tolower(m_seq[m_idx]);
    if (! (x=='a'||x=='c'||x=='g'||x=='t')) {
      //++m_idx;
      //skip over all sequential degenerate bp
      while (true) {
        if (m_idx == m_seq.size()) {
          return false;
        }
        x = tolower(m_seq[m_idx++]);
        if (x=='a'||x=='c'||x=='g'||x=='t') {
          --m_idx;
          break;
        }
      }
      //at this point m_seq[m_idx] should not be degenerate!
      bool status = start_over();
      if (!status) {
        return false;
      }
      output = m_cur;
      return true;
    }

    //case 4: we're good to go!

        uint64_t p;
        switch(m_seq[m_idx++]) {
          case 'a' :
          case 'A' : p = 0; 
                     break;
          case 'c' :
          case 'C' : p = 1; 
                     break;
          case 'g' :
          case 'G' : p = 2; 
                     break;
          case 't' :
          case 'T' : p = 3; break;
          default : assert(false);
        }
        m_cur = m_cur << 2;
        m_cur = m_cur & m_mask;
        m_cur = m_cur | p;
        output = m_cur;
      return true;
  }

private:
  std::string m_seq;
  unsigned short m_k;
  uint64_t m_idx;
  uint64_t m_cur;
  uint64_t m_mask;

  bool start_over() {
    m_cur = 0;
    bool good = false;
    unsigned short c = 0;
    bool kontinue = true;
    while (kontinue) {
      for (; m_idx<m_seq.size(); ++m_idx) {
        char x = tolower(m_seq[m_idx]);
        if (! (x=='a'||x=='c'||x=='g'||x=='t')) {
          c = 0;
        } else {
          ++c;
          if (c == m_k) {
            good = true;
            kontinue = false;
            break;
          }
        }
      }
      break;
    }

    if (!good) {
      return false;
    }

    ++m_idx;
    m_cur = 0;
    int z = 0;
    for (uint64_t j=m_idx-m_k; j<m_idx; j++, ++z) {
          uint64_t p;
          switch(m_seq[j]) {
            case 'a' :
            case 'A' : p = 0; break;
            case 'c' :
            case 'C' : p = 1; break;
            case 'g' :
            case 'G' : p = 2; break;
            case 't' :
            case 'T' : p = 3; break;
            default : assert(false);
          }

          m_cur = m_cur | p;
          if (z < m_k-1) {
            m_cur = m_cur << 2;
          }
        }
        return true;
      }
};

}

#endif
