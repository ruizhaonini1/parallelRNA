#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <chrono>
#include <algorithm>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::max;

typedef boost::multiprecision::int128_t long_num;

//#define NUMTHREADS 10

int dlen;
int chunklen;
int i,j,k;

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
  int numthreads = atoi(argv[2]);
  std::ifstream infile(filename);
  string line;
  std::getline(infile,line);
  int N = line.size();
  
  //omp_set_num_threads(NUMTHREADS);
  omp_set_num_threads(numthreads);
  //Convert the sequence into array of integers
  vector<int> seq;
  base_to_num(line, seq);

  //Build pairing matrix
  // pairing[i][j]=true if i and j can pair
  // COULD WE DO THIS IN PARALLEL TOO?
  vector<vector<bool> > pairing(N, vector<bool>(N,false));
  for (int i = 0; i < N; i++) {
    for (int j = i+4; j < N; j++) {
      if (can_pair(seq[i], seq[j]))
        pairing[i][j] = true;
    }
  }

  //Define dynamic programming matrix
  vector<vector<int> > S (N+1, vector<int>(N+1,0));

  
  auto start = std::chrono::steady_clock::now();
  // Dynamically change how many threads are being used as we get closer to the end of the matrix?
  // Like check if we can divide the diagonal nicely each loop, if we can't just get 1 or 2 threads to finish 
  for (int d =  1; d < N -2 ; d++){
    // How long is the diagonal
    dlen = N - d;
    
    if(dlen < numthreads){
        numthreads = ceil(numthreads / 4);
        omp_set_num_threads(numthreads);
    }
    
    // How big are the chunks for the different threads to work on?
    chunklen = ceil(dlen / numthreads);
    
    #pragma omp parallel private(i,j,k)
    {
      int tid = omp_get_thread_num();
      int starting_pos = (chunklen * tid) + d;
      int end_pos = starting_pos + chunklen;
      
      if (tid == (numthreads-1) || end_pos > N){
      //if(end_pos > N) {
	    end_pos = N;
      }

      i = chunklen*tid;
      for (j = starting_pos; j < end_pos; j++){
            int max_pair = S[i+1][j];
            for (k = i+4; k <=j; k++){
                if(pairing[i][k]){
                    max_pair = max(max_pair,  1+S[i+1][k-1]+S[k+1][j]);
                }
            }
            S[i][j] = max_pair;
            i++;
      }
    }
  }
  
//  for (int a = 0; a < N; a++){
//      for (int b = 0; b < N; b++){
//          cout << S[a][b] << " ";
//      }
//      cout << endl;
//  }
  
  // main loop
  for (int j = N - 2 ; j < N ; j++) {
    for (int i = j-1; i >= 0; i--) {
      int max_pair = S[i+1][j] ;
      for (int k = i+4; k <=j; k++){
        if (pairing[i][k]){
          max_pair = max(max_pair,  1+S[i+1][k-1] +S[k+1][j]);
        }
      }
      S[i][j] = max_pair;
    }
  }
  
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  // ouput the total number of secondary structures
  cout << S[0][N-1] << endl;
  cout << "Time to fill matrix: " << elapsed_seconds.count() << "s" << endl;
  return EXIT_SUCCESS;

}



