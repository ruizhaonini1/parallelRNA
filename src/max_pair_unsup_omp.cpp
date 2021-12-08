#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <omp.h>
#include <chrono>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::max;

int max_pair(const vector< vector<uint> > &S, const vector< vector<bool> > &pairing, 
    const uint i, const uint j){
  uint max_pair = S[i+1][j];
  for (int k = i+4; k <= j; k++){
    if(pairing[i][k]){
        max_pair = max(max_pair,1+S[i+1][k-1]+S[k+1][j]);
    }
  }
  return max_pair;
}

//Change RNA bases to numbers
void base_to_num(const string &line, vector<uint>& seq) {
  for (char c : line) {
    if (c == 'A') 
      seq.push_back(0);
    else if (c == 'C') 
      seq.push_back(1);
    else if (c == 'G') 
      seq.push_back(2);
    else if (c == 'U') 
      seq.push_back(3);
  }
  return;
}

// return true only for {A,U}, {C,G}, and {G,U} pairs 
bool can_pair(const int a, const int b) {
  return ((a+b==3) || (a+b==5));
}

int
main(int argc, const char **argv) {
  //Get input  sequence
  string filename = argv[1];
  uint numthreads = atoi(argv[2]);

  omp_set_num_threads(numthreads);
  std::ifstream infile(filename);

  string line;
  std::getline(infile,line);
  int N = line.size();
  //Convert the sequence into array of integers
  vector<uint> seq;
  base_to_num(line, seq);

  //Build pairing matrix
  // pairing[i][j]=true if i and j can pair
  vector<vector<bool> > pairing(N, vector<bool>(N,false));
  
  #pragma omp parallel for
  for (uint i = 0; i < N; i++) {
    for (uint j = i+4; j < N; j++) {
      if (can_pair(seq[i], seq[j]))
        pairing[i][j] = true;
    }
  }

  //Define dynamic programming matrix
  vector<vector<uint> > S (N+1, vector<uint>(N+1,0));

  auto start = std::chrono::steady_clock::now();
  // main, cubic time loop
  for (uint d =  1; d < N; d++){
    uint dlen = N - d;  // length of the diagonal wavefront
    
    #pragma omp parallel for
    for (uint i = 0; i < dlen; i++){
      uint j = i+d;
      S[i][j] = max_pair(S, pairing, i, j);
    }
  }
  
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  // ouput the total number of secondary structures
  cout << S[0][N-1] << endl;
  cout << "Time to fill matrix: " << elapsed_seconds.count() << "s" << endl;
  return EXIT_SUCCESS;

}



