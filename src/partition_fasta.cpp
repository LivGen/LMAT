#include <iostream>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <vector>
#include <string>
#include "Utils.hpp"

using namespace std;

void usage(char *argv[]) {
    cerr << "Usage:\n"
     " -i <string>  - input fasta_fn (required)\n"
     " -o <string>  - base output directory (required)\n"
     " -n <int>     - num files\n"
     "function: breaks the input fasta file(s) into multiple files\n"
     "using a round-robin approach.\n";
}


int main(int argc, char *argv[]) {
    cout << "invocation: ";
    for (int j=0; j<argc; j++) {
        cout << argv[j] << " ";
    }
    cout << endl;

    string inf;
    string output_directory;
    int n_files = 0;

    int count = 0;
    const string opt_string="i:o:n:";

    signed char c;
    while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
        switch (c) {
            case 'i' :
                inf = optarg;
                ++count;
                break;
            case 'o':
                ++count;
                output_directory = optarg;
                break;
            case 'n':
                n_files = atoi(optarg);
                ++count;
                break;
            default:
                cerr << "Unrecognized option: " << c << endl;
                exit(9);
                break;
        }
    }
  if (count != 3) {
    cout << "usage: " << argv[0] 
         << " -i input_filename -o output_base_dir -n num_output_files\n";
    exit(9);
  }

  stringstream err;
    
  // Open input file
  ifstream in(inf.c_str());
  if (!in) {
    err << "failed to open " << inf << " for reading";
    LMAT_ERR(err.str());
  }

  // Make output directory; assume we're not making a crazy number
  // of partition blocks, so all fit in a single directory nicely
  stringstream d;
  d << "mkdir -p " << output_directory;
  int r = system(d.str().c_str());
  if (r == -1) {
    err << "system call to mkdir failed\n";
    LMAT_ERR(err.str());
  }

  // Open output files
  vector<ofstream> out(n_files);
  for (size_t j=0; j<n_files; j++) {
    stringstream ss;
    ss << output_directory << "/" << j;
    out[j].open(ss.str().c_str());
    if (!out[j]) {
      err << "failed to open output file: " << ss.str();
      LMAT_ERR(err.str());
    }
  }

  // Use round-robin approach to partition the sequences
  string header;
  string seq;
  int cur_file = 0;
  while (in >> header >> seq) {
    out[cur_file] << header << '\n' << seq << '\n';
    ++cur_file;
    if (cur_file == out.size()) {
      cur_file = 0;
    }
  }

  // Close output files
  for (size_t j=0; j<out.size(); j++) {
    out[j].close();
  }
}

