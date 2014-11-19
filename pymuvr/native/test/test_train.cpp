/* Use with somthing like
 *
 * g++ -o test_train -Wall -iquote ../include ../src/convolved_spike_train.cpp test_train.cpp
 */

#include<vector>
#include<iostream>

#include "convolved_spike_train.hpp"

using namespace std;

void print_train(ConvolvedSpikeTrain t){
  cout << t.tau() << ' ' << t.size() << ' ' << t.square_norm() << endl;
  for (unsigned int i=0; i<t.size(); ++i){
    cout << t.spikes()[i] << ' ' << t.exp_pos()[i] << ' ' << t.exp_neg()[i] << ' ' << t.markage()[i] << endl;
  }
  cout << endl;
}

int main(){
  double tau = 0.1;

  vector<double> spikes;
  spikes.push_back(1.2);
  spikes.push_back(1.4);
  spikes.push_back(1.5);

  ConvolvedSpikeTrain ctrain = ConvolvedSpikeTrain(spikes, tau);
  print_train(ctrain);

  
  ctrain.UpdateConvolution(0.3);
  print_train(ctrain);

  vector<double> spikes2;
  ConvolvedSpikeTrain ctrain2 = ConvolvedSpikeTrain(vector<double>(), 0);
  print_train(ctrain2);

  vector< ConvolvedSpikeTrain > observation (4);
  print_train(observation[2]);

  observation[2] = ctrain;
  print_train(observation[2]);

  return 0;
}
