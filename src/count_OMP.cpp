#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <omp.h>

using std::string;
using std::cout;
using std::endl;
using std::vector;

typedef boost::multiprecision::int128_t long_num;

#define NUMTHREADS 2

int dlen;
int chunklen;

//Change RNA bases to numbers
void base_to_num(const string line, vector<int>& seq) {
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
  std::ifstream infile(filename);
  string line;
  std::getline(infile,line);
  int N = line.size();
  
  omp_set_num_threads(NUMTHREADS);

  //Convert the sequence into array of integers
  vector<int> seq;
  base_to_num(line, seq);

  //Build pairing matrix
  // pairing[i][j]=true if i and j can pair
  vector<vector<bool> > pairing(N, vector<bool>(N,false));
  for (int i = 0; i < N; i++) {
    for (int j = i+4; j < N; j++) {
      if (can_pair(seq[i], seq[j]))
        pairing[i][j] = true;
    }
  }

  //Define dynamic programming matrix
  vector<vector<long_num> > S (N+1, vector<long_num>(N+1,0));

  // initialize with base cases
  for (int i =0; i < N+1; i++) {
    S[i][i] = 1;
  }
  for (int i =1; i < N+1; i++) {
    S[i][i-1] = 1;
  }
  
  // put this omp pragma into a for loop to go through all of the diagonals
  // Move main for loop into this omp
  // determine chunk of matrix to work on based on tid
  
  for (int d = 0 ; d < N - 2; d++){
    dlen = N - d;
    chunklen = ceil(dlen / NUMTHREADS);
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int starting_pos = (chunklen * tid) + d;
      int end_pos = starting_pos + chunklen;
      
      if (tid == (NUMTHREADS-1)){
          end_pos = N;
      }
      
      for (int j = starting_pos; j < end_pos; j++){
          for (int i = d*tid; i < starting_pos; i++){
              S[i][j] = S[i+1][j];
              for (int k = i+4; k <=j; k++){
                S[i][j] += (int)pairing[i][k]*S[i+1][k-1]*S[k+1][j];
              }
          }
      }
      
    }
  }
  
  // main loop
  for (int j = N - 2 ; j < N ; j++) {
    for (int i = j-1; i >= 0; i--) {
      S[i][j] = S[i+1][j];
      for (int k = i+4; k <=j; k++){
        S[i][j] += (int)pairing[i][k]*S[i+1][k-1]*S[k+1][j];
      }
    }
  }
    
  // ouput the total number of secondary structures
  cout << S[0][N-1] << endl;

  return EXIT_SUCCESS;

}



