/* Use with somthing like
 *
 * g++ -o test_product -Wall -iquote ../include ../src/convolved_spike_train.cpp ../src/van_rossum_multiunit.cpp test_product.cpp
 */

#include "convolved_spike_train.hpp"
#include "van_rossum_multiunit.hpp"

#include<vector>
#include<iostream>
#include<limits>
#include<cstdlib>

using namespace std;

void print_train(ConvolvedSpikeTrain t){
  cout << t.tau << ' ' << t.size << ' ' << t.square_norm << endl;
  for (unsigned int i=0; i<t.size; ++i){
    cout << t.spikes[i] << ' ' << t.exp_pos[i] << ' ' << t.exp_neg[i] << ' ' << t.markage[i] << endl;
  }
  cout << endl;
}

int main(){
  double tau = 0.3;

  // create spike trains
  vector<double> spikes1;
  vector<double> spikes2;
  vector<double> spikes3;
  vector<double> spikes4;

  spikes1.push_back(0.2);
  spikes1.push_back(0.4);
  spikes1.push_back(0.5);
  spikes2.push_back(0.1);
  spikes2.push_back(0.3);
  spikes2.push_back(0.9);
  spikes3.push_back(0.2);
  spikes3.push_back(0.3);
  spikes4.push_back(0.4);
  spikes4.push_back(0.8);
  spikes4.push_back(0.95);

  vector< vector<ConvolvedSpikeTrain> > observations (2, vector<ConvolvedSpikeTrain>(3));
  observations[0][0] = ConvolvedSpikeTrain(spikes1, tau);
  observations[0][1] = ConvolvedSpikeTrain(spikes2, tau);
  observations[1][0] = ConvolvedSpikeTrain(spikes3, tau);
  observations[1][1] = ConvolvedSpikeTrain(spikes4, tau);

  print_train(observations[0][0]);
  print_train(observations[1][0]);

  // allocate distance matrix
  int n_obs = observations.size();
  double** g_matrix = new double*[n_obs];
  for (int i=0; i<n_obs; ++i) {
    g_matrix[i] = new double[n_obs];
  }
  
  // compute inner product
  inner_product(observations, 0.2, g_matrix);

  // print inner product matrix
  for (int i=0; i<n_obs; ++i){
    for (int j=0; j<n_obs; ++j){
      cout << g_matrix[i][j] << " ";
    }
    cout << endl;
  }
  return 0;
}
