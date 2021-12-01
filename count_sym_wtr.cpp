#include "read_constants_nupack.hpp"
#include "constants.hpp"
#include "utils.hpp"
#include "parameters.hpp"
#include <boost/multiprecision/cpp_int.hpp>
// #include "boost/multiprecision/cpp_int.hpp"

using std::string;
using std::stringstream;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;
using std::make_pair;
using std::runtime_error;
using std::to_string;
using std::map;
using std::tuple;
using std::make_tuple;
using std::get;
using nlohmann::json;
using std::unordered_map;
using std::reference_wrapper;
using namespace ss_names;


struct
count_matrices{
  s_matrix Cb;
  s_matrix Ce;
  s_matrix C;
  s_matrix Cbx;
  s_matrix Crb;
  s_matrix Cre;
  s_matrix Cr;
  s_matrix Crb1;
  s_matrix Crb2;
  s_matrix S;
  s_matrix Se;
  s_matrix Sr;
  s_matrix Sb;


  void
  initialize(size_t m_size){
    Cb = s_matrix(m_size, vector<ull>(m_size, 0));
    Cbx = s_matrix(m_size, vector<ull>(m_size, 0));
    Ce = s_matrix(m_size, vector<ull>(m_size, 0));
    C = s_matrix(m_size, vector<ull>(m_size, 0));
    Crb = s_matrix(m_size, vector<ull>(m_size, 0));
    Crb1 = s_matrix(m_size, vector<ull>(m_size, 0));
    Crb2 = s_matrix(m_size, vector<ull>(m_size, 0));
    Cre = s_matrix(m_size, vector<ull>(m_size, 0));
    Cr = s_matrix(m_size, vector<ull>(m_size, 0));

    S = s_matrix(m_size, vector<ull>(m_size, 0));
    Se = s_matrix(m_size, vector<ull>(m_size, 0));
    Sr = s_matrix(m_size, vector<ull>(m_size, 0));
    Sb = s_matrix(m_size, vector<ull>(m_size, 0));

    for (pos i = 0; i < m_size; i++) {
      S[i][i] = 1;
    }
    for (pos i = 1; i < m_size; i++) {
      S[i][i-1] = 1;
    }


  }

};




void
update_S(const vector<base> &seq,
         const vector<pos> &nicks,
         const vector<vector<bool> > &eta,
         const pos i, const pos j,
         count_matrices &cm){

  ull count = eta[i][i+1] * cm.S[i+1][j];
  for (pos k = i+1; k <= j; k++){
    ull factor = (k == j) ? 1 : eta[k][k+1] * cm.S[k+1][j];
    count += bases_can_pair(seq, eta, i, k)
      * (eta[i][i+1] * cm.S[i+1][k-1] * eta[k-1][k] + cm.Se[i][k])
      * factor;
  }
  cm.S[i][j] = count;
}



void
update_Se(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          count_matrices &cm){
  if (eta[i][i+1] || eta[j-1][j] || (j == i + 1) ){
    ull count = 0;
    if (!(eta[i][i+1]) || !(eta[j-1][j])) {
      count = cm.S[i+1][j-1];
    }
    else {
      for (pos d : nicks) {
        if ((i <= d) && (d < j)) {
          count += cm.S[i+1][d]*cm.S[d+1][j-1];
        }
      }
    }
    cm.Se[i][j] = count;
  }
}

void
update_Sr(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          count_matrices &cm){

  ull count = eta[i][i+1] * cm.Sr[i+1][j];
  for (pos k = i+1; k < j - 1; k++){
    count += bases_can_pair(seq, eta, i, k)
      * (eta[i][i+1] * cm.S[i+1][k-1] * eta[k-1][k] + cm.Se[i][k])
      * eta[k][k+1] * cm.Sr[k+1][j];
  }
  for (pos k = i+1; k <= j; k++){
    if (bases_can_pair(seq[i], seq[k])) {
      ull factor = (k == j) ? 1 : eta[k][k+1] * cm.S[k+1][j];
      ull sub_count = 0;
      sub_count += (eta[i][i+1] * cm.S[i+1][k-1] * eta[k-1][k])
        * factor;
    }
  }

  cm.S[i][j] = count;
}



void
update_Cb(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          count_matrices &cm){
  // printf("Cb(%d, %d)\n", i, j);
  if ((!eta[i][i+1]) && (!eta[j-1][j]) && (i + 1 != j ))
    return;

  if (!bases_can_pair((size_t)seq[i], (size_t)seq[j])) {
    return;
  }
  if ( (eta[i][j]) &&
       (j <= i + hairpin_property::min_loop_size_int) )
    return;

  cm.Cb[i][j] = cm.Ce[i][j] + eta[i][i+1]*eta[j-1][j]*cm.C[i+1][j-1];
}

void
update_Cbx(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          count_matrices &cm){
  ull count = 0;
  for (pos k = j; (k > i) && (eta[k][j]); k--) {
    count += cm.Cb[i][k];
  }

  cm.Cbx[i][j] = count;
}


void
update_Ce(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          count_matrices &cm){
  // printf("Ce(%d, %d)", i, j);
  if ((!eta[i][i+1]) && (!eta[j-1][j]) && (i + 1 != j ))
    return;
  ull count = 0;
  if ((!eta[i][i+1]) || (!eta[j-1][j])) {
    count += cm.C[i+1][j-1];
  }
  else {
    for (pos n : nicks) {
      if ((i <= n) && (n < j)) {
        count += cm.C[i+1][n]*cm.C[n+1][j-1];
      }
    }
  }

  cm.Ce[i][j] = count;
}


void
update_C(const vector<base> &seq,
         const vector<pos> &nicks,
         const vector<vector<bool> > &eta,
         const pos i, const pos j,
         count_matrices &cm){
  ull count = eta[i][j];
  // for (pos l = j; (l > i) && (eta[l][j]); l--){
  //   count += cm.Cb[i][l];
  // }
  for (pos k = i; (k < j); k++) {
    ull factor1 = (k == i) ? 1 : cm.C[i][k-1];
    bool factor2 = (k == i) ? true : eta[k-1][k];
    if (factor2){
        count += factor1*cm.Cbx[k][j];
    }
  }
  cm.C[i][j] = count;
}




void
update_Cre(const vector<base> &seq,
           const vector<pos> &nicks,
           const vector<vector<bool> > &eta,
           const pos i, const pos j,
           count_matrices &cm){
  if ((!eta[i][i+1]) && (!eta[j-1][j]) )
    return;
  ull count = 0;
  if (eta[i][i+1]) {
    for (pos k = i+2; k < j; k++ ) {
      count += cm.Crb2[i+1][k]*cm.Ce[k][j];
    }
  }
  if (eta[j-1][j]) {
    for (pos k = i+1; k < j-1; k++ ) {
      count += cm.Ce[i][k]*cm.Crb1[k][j-1];
    }
  }



  cm.Cre[i][j] = count;
}

void
update_Crb(const vector<base> &seq,
           const vector<pos> &nicks,
           const vector<vector<bool> > &eta,
           const pos i, const pos j,
           count_matrices &cm){
  // printf("Cb(%d, %d)\n", i, j);
  if ((!eta[i][i+1]) && (!eta[j-1][j]) )
    return;
  if (!bases_can_pair(seq[i], seq[j])) {
    return;
  }
  if ((!eta[i][i+1]) || (!eta[j-1][j]))
    cm.Crb[i][j] = cm.Cre[i][j];
  else
    cm.Crb[i][j] = cm.C[i+1][j-1] + cm.Cre[i][j] + cm.Cr[i+1][j-1];
}


void
update_Crb1(const vector<base> &seq,
            const vector<pos> &nicks,
            const vector<vector<bool> > &eta,
            const pos i, const pos j,
            count_matrices &cm){
  ull count = 0;
  for (pos k = i + 1; k <= j; k++) {
    ull factor = (k == j) ? 1: cm.C[k+1][j];
    bool factor2 = (k == j) ? true : eta[k][k+1];
    count += cm.Crb[i][k]*factor2*factor;
  }

  cm.Crb1[i][j] = count;
}

void
update_Crb2(const vector<base> &seq,
            const vector<pos> &nicks,
            const vector<vector<bool> > &eta,
            const pos i, const pos j,
            count_matrices &cm){
  ull count = 0;
  for (pos k = i; k < j ; k++) {
    ull factor = (k == i) ? 1 : cm.C[i][k-1];
    bool factor2 = (k == i) ? true: eta[k-1][k];
    count += factor*factor2*cm.Crb[k][j];
  }

  cm.Crb2[i][j] = count;
}


void
update_Cr(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          count_matrices &cm){
  ull count = 0;
  for (pos k = i; k < j ; k ++){
    ull factor1 = (k == i) ? 1 : eta[k-1][k];
    ull factor2 = (k == i) ? 1 : cm.C[i][k-1];
    count += factor1*factor2*cm.Crb1[k][j];
  }


  cm.Cr[i][j] = count;
}


void
update_Cy(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &C, const s_matrix &Cb,
          s_matrix &Cy){
  ull count = 0;
  for (pos k = j-2; k-1 > i+1; k--) {
    bool proceed = true;
    if (!eta[k-1][k]) proceed = false;
    if (proceed){
      count += Cb[i][k-1]*C[k][j];
    }
  }
  Cy[i][j] += count;
}


void
update_C(const vector<base> &seq,
         const vector<pos> &nicks,
         const vector<vector<bool> > &eta,
         const pos i, const pos j,
         const s_matrix &Cx,
         s_matrix &C){
  ull count = eta[i][j];
  for (pos k = i; k <= j - 2; k++) {
    bool proceed = true;

    if ( (k > 0) && (!eta[k-1][k])) proceed = false;
    if (proceed) {
      if (k == i) count += Cx[k][j];
      else count += C[i][k-1]*Cx[k][j];
    }
  }
  C[i][j] += count;
}

void
update_Cx(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cb,
          s_matrix &Cx){
  if ( (eta[i][j]) &&
       (j <= i + hairpin_property::min_loop_size_int) )
    return;
  ull count = 0;
  for (pos k = j; eta[k][j] && (i < k); k--){
    count += Cb[i][k];
  }
  Cx[i][j] += count;
}

void
update_Cm(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cx,
          s_matrix &Cm){
  ull count = 0;
  for (pos k = i; eta[i][k] && (k < j); k++){
    count += Cx[k][j];
  }
  for (pos k = j-1; (k > i+1); k--){
    bool proceed = true;
    if (!eta[k-1][k]) proceed = false;
    if (proceed)
      count += Cm[i][k-1]*Cx[k][j];
  }
  Cm[i][j] += count;
}

void
update_Cf(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cb,
          s_matrix &Cf){
  ull count = 0;
  for (pos k = j; eta[k][j] && (k > i); k--){
    // printf ("adding Cb[%d][%d] = %ld to Cf[%d][%d]\n", i,k, Cb[i][k], i, j);
    count += Cb[i][k];
  }
  Cf[i][j] += count;


}

void
update_Ce(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cm,
          s_matrix &Ce){
  if ((i == 0) || (j == seq.size()-1)) return;
  ull count = 0;
  for (pos n : nicks) {
    if ((i <= n) && (n <= j)) {
      ull count1 = (eta[i-1][n-1] == 0) ? 1 : 0;
      count1 += Cm[i][n-1];
      ull count2 = (eta[n+1][j+1] == 0) ? 1 : 0;
      count2 += Cm[n+1][j];
      count += count1*count2;
    }
  }
  Ce[i][j] = count;
}



void
update_Cbs(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cbp,
          const s_matrix &Cm,
          const s_matrix &Cep,
          const s_matrix &Cp,
          const s_matrix &Cmp,
          const s_matrix &Ce,
          s_matrix &Cbs){
  if (!(eta[i][i+1]) && !(eta[j-1][j]) && (i + 1 != j -1))
    return;

  if (!bases_can_pair((size_t)seq[i], (size_t)seq[j])) {
    return;
  }

  ull count = 0;
  // center interior loop or multiloop
  if (eta[i][j] == 0) count += 1;

  // interior loops
  for (pos k = i+1 ; (eta[i][k]==0) && (k<j-1); k++) {
    count += Cbp[k][j-1];
  }

  // exterior loops
  for (pos k = i + 1; k < j; k++) {
    count += ((eta[i][k-1] == 0) + Cm[i+1][k-1])*Cep[k][j-1];
    count += Ce[i+1][k-1]*Cp[k][j-1];
  }

  // center nontrivial multiloop
  count += Cm[i+1][j-1];

  // non-center multiloops
  for (pos k = i + 1; k < j; k++) {
    count += ((eta[i][k-1] == 0) + Cm[i+1][k-1])*Cmp[k][j-1];
    count += Cm[i+1][k-1]*Cp[k][j-1];
  }

  Cbs[i][j] = count;

}


void
update_Cbp(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cbs,
          s_matrix &Cbp){
  ull count = 0;
  for (pos k = j; (eta[k][j] == 0) && (k > i); k--) {
    count += Cbs[i][k];
  }
  Cbp[i][j] = count;
}

void
update_Cep(const vector<base> &seq,
           const vector<pos> &nicks,
           const vector<vector<bool> > &eta,
           const pos i, const pos j,
           const s_matrix &Cbs,
           const s_matrix &Ce,
           s_matrix &Cep){
  ull count = 0;
  for (pos k = j-1; k > i; k--) {
    count += Cbs[i][k]*Ce[k+1][j];
  }
  Cep[i][j] = count;
}

void
update_Cp(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cbs,
          const s_matrix &Cm,
          s_matrix &Cp){
  ull count = 0;
  count += Cbs[i][j];
  for (pos k = j-1; k > i; k--) {
    count += Cbs[i][k]*((eta[k+1][j]==0)+Cm[k+1][j]);
  }
  Cp[i][j] = count;
}

void
update_Cmp(const vector<base> &seq,
          const vector<pos> &nicks,
          const vector<vector<bool> > &eta,
          const pos i, const pos j,
          const s_matrix &Cbs,
          const s_matrix &Cm,
          s_matrix &Cmp){
  ull count = 0;

  for (pos k = j-1; k > i; k--) {
    count += Cbs[i][k]*Cm[k+1][j];
  }
  Cmp[i][j] = count;
}

int
main(int argc, const char **argv) {
  try {

    string outfile;
    string param_name = "rna1995";
    string param_dir;
    // run mode flags
    bool VERBOSE = false;
    double T = 37.0;

    const string description =
      "Calculate free energy of secondary structure for a particular\
      sequence. The input file must contain the sequence in the first line\
      and the secondary structure in the dot-parenthesis notation in the\
      seconda line.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           description, "<sequence-structure-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("params", 'p', "parameter name (default: rna1995)",
                      false, param_name);
    opt_parse.add_opt("paramdir", 'd', "parameter directory",
                      false, param_dir);
    opt_parse.add_opt("temp", 't', "temeprature in celcius (default: 37)",
                      false, T);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    
    
    uint num_strands;
    vector<vector<base> > strands;
    vector<uint> strand_ids;
    parse_pfunc_input_file(input_file, num_strands, strands, strand_ids);
    vector<uint> sym({1});
    uint n_strands = strand_ids.size();
    get_sym(strand_ids, sym);
    uint max_sym = sym.back();
    vector<base> tot_seq;
    vector<pos> nicks;
    vector<uint> gam;
    get_seq(strands, strand_ids, tot_seq, nicks, gam);
    unsigned size = tot_seq.size();

    vector<vector<bool> > eta(size,  vector<bool>(size,0));
    get_eta(tot_seq, nicks, eta);
    vector<vector<bool> > p(size,  vector<bool>(size,0));
    get_p(tot_seq, nicks, eta, p);



    size_t seq_size = size;
    const vector<base>::const_iterator first = tot_seq.begin();
    const vector<base>::const_iterator last = tot_seq.begin() + seq_size;
    vector<base> seq(first, last);
    count_matrices cm;
    // for (pos n : nicks) printf("%d\n", n);

    cm.initialize(seq_size);
    // for (pos i = 0; i < seq_size; i++){
    //   for (pos j = i; j < seq_size; j++) {
    //     printf("eta[%d][%d] = %d\n", i, j, (int) eta[i][j]);
    //   }
    // }

    for (pos n : nicks) cm.C[n+1][n] = 1;
    for (pos i = 0; i < seq_size-1; i++) {
      cm.C[i][i] = 1;
    }
    for (pos i = 1; i < seq_size-1; i++) {
      cm.C[i][i-1] = 1;
    }


    for (pos j = 0; j < (pos)seq_size ; j++) {
      for (pos i = j - 1; i != (unsigned)-1 ; i--) {
        update_Ce(seq, nicks, eta, i, j, cm);
        update_Cb(seq, nicks, eta, i, j, cm);
        update_Cbx(seq, nicks, eta, i, j, cm);
        update_C(seq, nicks, eta, i, j, cm);
        update_Cre(seq, nicks, eta, i, j, cm);
        update_Crb(seq, nicks, eta, i, j, cm);
        update_Crb1(seq, nicks, eta, i, j, cm);
        update_Crb2(seq, nicks, eta, i, j, cm);
        update_Cr(seq, nicks, eta, i, j, cm);

        update_Se(seq, nicks, eta, i, j, cm);
        update_S(seq, nicks, eta, i, j, cm);
        // printf("S[%d][%d] = ", i,j);
        // std::cout << cm.S[i][j] << std::endl;
        // printf("Se[%d][%d] = ", i,j);
        // std::cout << cm.Se[i][j] << std::endl;

        // printf("Cre[%d][%d] = %lld\n", i,j,cm.Cre[i][j]);
        // printf("Crb[%d][%d] = %lld\n", i,j,cm.Crb[i][j]);
        // printf("Cr[%d][%d] = %lld\n", i,j,cm.Cr[i][j]);

        // bool proceed = true;
        // if ((i > 0) && (eta[i-1][i] == 0)&&
        //     (j < (pos)seq_size - 1) && (seq[j+1] == base_K))
        //   proceed = false;
        // if ((i == 0) && (j < (pos)seq_size - 1) && (seq[j+1] == base_K))
        //   proceed = false;
        // if ((j == (pos) seq_size -1) && (i > 0) && (seq[i-1] == base_K))
        //   proceed = false;
        // if (proceed) {
        //   update_Cb(seq, nicks, eta, i, j, Cm, Cx, Cf, Ce, Cb);
        //   update_Cx(seq, nicks, eta, i, j, Cb, Cx);
        //   update_Cm(seq, nicks, eta, i, j, Cx, Cm);
        //   update_Cf(seq, nicks, eta, i, j, Cb, Cf);
        // }

      }
    }
    ull count = cm.C[0][seq_size-1];

    std::cout << cm.S[0][seq_size-1] << std::endl;
    // printf("tot: %lld\n", cm.S[0][seq_size-1]);

    vector<ull> sym_structures;
    sym_structures.push_back(count);
    for (uint k = 0; k < sym.size(); k++) {

      uint n_sym = sym[k];
      if (n_sym == 1) continue;
      size_t seq_size = size;
      seq_size = nicks[n_strands/n_sym - 1]+1;
      sym_structures.push_back(cm.Cr[0][seq_size-1]);
    }
    ull count_corrected = 0;
    for (pos i = sym.size()-1; i != (unsigned) -1; i--) {
      ull c_count = sym_structures[i];
      for (pos j = i+1; j < (pos)sym.size(); j++) {
        if (sym[j]%sym[i] == 0) {
          c_count -= sym_structures[j]*max_sym/sym[j];
        }
      }
      sym_structures[i] = sym[i]*c_count/max_sym;
      count_corrected += sym_structures[i];
    }



    // output
    if (outfile == "") {
      printf("Uncorrected number of structures\n");
      std::cout << "Total\t: " << count << std::endl;
      // printf("Total\t: %lld\n",  count);

      for (uint i = 0; i < sym.size(); i++) {
        std::cout << "sym = " << sym[i] << "\t: " <<
          sym_structures[i]*max_sym/sym[i] << std::endl;
      }
      printf("Corrected number of structures\n");
      std::cout << "Total\t: " << count_corrected << std::endl;

      for (uint i = 0; i < sym.size(); i++) {
        std::cout << "sym = " << sym[i] << "\t: " <<
          sym_structures[i] << std::endl;

      }

    }
    else {
      std::fstream fs;
      fs.open (outfile, std::fstream::out);
      fs << "Uncorrected number of structures\n";
      fs << "Total\t: " << count << std::endl;
      // printf("Total\t: %lld\n",  count);

      for (uint i = 0; i < sym.size(); i++) {
        fs << "sym = " << sym[i] << "\t: " <<
          sym_structures[i]*max_sym/sym[i] << std::endl;
      }
      fs <<"Corrected number of structures\n";
      fs << "Total\t: " << count_corrected << std::endl;

      for (uint i = 0; i < sym.size(); i++) {
        fs << "sym = " << sym[i] << "\t: " <<
          sym_structures[i] << std::endl;

      }
      fs.close();

      // FILE *fp;
      // fp = fopen(outfile.c_str(),"w");
      // for (uint i = 0; i < sym.size(); i++) {
      //   fprintf(fp, "sym = %d: %ld\n", sym[i], sym_structures[i]);
      // }
      // fclose(fp);
    }
  }
  catch (runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}



//TODO: Make function to check if the input is valid
