#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::max;



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
  vector<vector<int> > S (N+1, vector<int>(N+1,0));

  // initialize with base cases
  //for (int i =0; i < N+1; i++) {
    //S[i][i] = 0;
  //}
  //for (int i =1; i < N+1; i++) {
    //S[i][i-1] = 1;
  //}

  // main loop
  for (int j = 0 ; j <N ; j++) {
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
    
  // ouput the total number of secondary structures
  cout << S[0][N-1] << endl;

  return EXIT_SUCCESS;

}



