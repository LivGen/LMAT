#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <cassert>
#include <string>
#include <string.h>
#include <stdint.h>
#include "all_headers.hpp"
#include "Encoder.hpp"
#include <ext/hash_set>
#include <version.h>

using namespace std;
using namespace metag;

void usage(char *argv[]) {
  cerr << "Usage:\n"
       << " -i <string>  - input fasta_fn\n"
       << " -k <int>     - kmer length\n"
       << " -o <string>  - output filename\n"
       << " -l <int>     - prefix length (as in: the string representation of the prefix)\n"
       << " -f <int>     - prefix\n"
       << " -m <int>     - initial memory; optional; default: 100,000,000\n"
       << " -h           - print help and exit\n\n";
}


uint64_t write_binary(map<uint64_t, set<uint32_t> > &hash, int mer_len, const char *fn);

int main(int argc, char *argv[]) {

  cout << "invocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;

  string input_fn, output_fn;
  int mer_len = 0;
  bool prn_help = false;
  uint64_t prefix = 0;
  int p_len = 0;
  int quit_early = 0;
  uint64_t mem = 100000000;

  int count = 0;
  const string opt_string="i:k:h o:p:l:q:f:m:V";
  char c;
  while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
    switch (c) {
    case 'm' :
      mem = atol(optarg);
      break;
    case 'f':
      ++count;
      prefix = atoi(optarg);
      break;
    case 'q':
      quit_early = atoi(optarg);
      break;
    case 'i':
      input_fn = optarg;
      ++count;
      break;
    case 'k':
      mer_len = atoi(optarg);
      ++count;
      break;
    case 'h':
      prn_help = false;
      break;
    case 'o':
      ++count;
      output_fn = optarg;
      break;
      break;
    case 'l':
      ++count;
      p_len = atoi(optarg);
      break;
    case 'V':
      cout << "LMAT version " << LMAT_VERSION  << "\n";
	exit(0);
    default:
      cerr << "Unrecognized option: "<<c<<", ignore."<<endl;
      prn_help = true;
      break;
    }
  }

  if (prn_help || count != 5) {
    usage(argv);
    exit(-1);
  }

  uint64_t kmerid, kmerid2, rc, kmer;
  string header, seq, mer, prefix_str;
  int moveme = (mer_len*2) - (p_len*2);
  map<uint64_t, set<uint32_t> > hash;

  ifstream in;
  cout << "opening: " << input_fn << endl;
  in.open(input_fn.c_str());
  if (!in) {
    cerr << "failed to open " << input_fn << " for reading\n";
    exit(1);
  }
  int s_ct = 0;
  int gid = 0; //stop compiler complaints
  StopWatch clk;
  clk.start();
  while (true) {
      header = "";
      in>>header;
      if (header.size() == 0) {
        break;
      } else {
        ++s_ct;
        if (header[0] != '>') {
          cerr << "header[0] != '>' for line number " << s_ct*2 << endl;
        }

        cerr << s_ct << " header: " << header << endl;
        assert(header[0] == '>');
        string id_s = header.substr(1);
        gid = atoi(id_s.c_str());
        in >> seq;
      }

      //exit condition, for development and testing
      if (quit_early && s_ct == quit_early) {
        in.close();
      }  

    //hash kmers and genome IDs
    Encoder e(seq, mer_len);
    while (e.next(kmerid)) {
      rc = Encoder::rc(kmerid, mer_len);
      kmer = kmerid < rc ? kmerid : rc;
      kmerid2 = kmer >> moveme;
      if (prefix == kmerid2) {
        hash[kmer].insert(gid);
      }
    }
  }
  double compute_time = clk.stop();

  clk.reset();
  clk.start();
  char buf[1024];
  sprintf(buf, "%s.%d", output_fn.c_str(), (int)prefix);
  cout << "writing to: " << buf << endl;
  write_binary(hash,  mer_len, buf);
  double write_time = clk.stop();

  cerr << endl << "time to write kmers: " << write_time << endl;

  cout << "sizeof(uint32_t): " << sizeof(uint32_t) << endl;
  cout << "invocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;
  uint64_t kmer_ct = hash.size();

  cout << "proc/[pid]/status:\n";
  pid_t p = getpid();
  sprintf(buf, "cat /proc/%d/status", p);
  system(buf);

  cout << "kmer_ct: " << kmer_ct << " compute time: " << compute_time << " write_time: " << write_time << endl;
  cout << "SUCCESS!\n";
}



uint64_t write_binary(map<uint64_t, set<uint32_t> > &hash, int mer_len, const char *fn) {
    KmerFileMetaData metadata;
    metadata.setSize(hash.size());
    metadata.setKmerLength(mer_len);
    metadata.setDefaultDataStart();
    cout << "opening: " << fn << endl;
    FILE *out = fopen(fn, "wb");
    assert(out);
    metadata.write(out);
    uint64_t sanity = ~0;
    uint32_t  id;
    int  count;
    uint64_t kmer_ct = 0;
    uint64_t kmer;
    string mer;
    //for (__gnu_cxx::hash_map<uint64_t, set<uint32_t> >::iterator t = hash.begin(); t != hash.end(); t++) {
    for (map<uint64_t, set<uint32_t> >::iterator t = hash.begin(); t != hash.end(); t++) {
      ++kmer_ct;
      kmer = t->first;
      assert(fwrite(&kmer, 8, 1, out) == 1);
      count = t->second.size();
      assert(fwrite(&count, 4, 1, out) == 1);
      for (set<uint32_t>::iterator t2 = t->second.begin(); t2 != t->second.end(); t2++) {
        id = *t2;
        assert(fwrite(&id, 4, 1, out) == 1);
      }
      if (kmer_ct % 1000 == 0) {
        assert(fwrite(&sanity, 8, 1, out) == 1);
      }
    }
    fclose(out);

    return kmer_ct; 
}

