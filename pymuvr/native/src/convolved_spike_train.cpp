#include "convolved_spike_train.hpp"

#include<cstdlib>
#include<stdexcept>
#include<vector>
#include<cmath>
#include<iostream>

/* Check if we're compiling with Visual Studio */
#ifdef _MSC_VER
#include <float.h>
#endif

/* Visual Studio does not have isfinite; use _finite() instead */
#ifdef _MSC_VER
bool isfinite(long double x)
{
    return _finite(x);
}
#endif

using namespace std;

ConvolvedSpikeTrain::ConvolvedSpikeTrain(): tau_(), spikes_(), size_(), exp_pos_(), exp_neg_(), markage_(){}
    
ConvolvedSpikeTrain::ConvolvedSpikeTrain(vector<double> spikes, double tau) : tau_(tau), spikes_(spikes), size_(spikes.size()), exp_pos_(size_), exp_neg_(size_), markage_(size_){
  UpdateConvolution(tau);
}

void ConvolvedSpikeTrain::UpdateConvolution(double tau){
  tau_ = tau;
  UpdateExponentialVectors();
  UpdateMarkageVector();
}

void ConvolvedSpikeTrain::UpdateExponentialVectors(){
  if (size_ > 0){
    for (unsigned int i=0; i<size_; ++i){
      exp_pos_[i] = exp((long double)(spikes_[i]/tau_));
      exp_neg_[i] = exp((long double)(-spikes_[i]/tau_));
    }

    /*Check if any over/underflow occurred while calculating the
      exponentials*/
    if (exp_neg_.front()==0 || !isfinite(exp_pos_.back())){
      throw overflow_error("tau is too small compared to the spike times. Please use a larger value for tau, or shorter spike trains.");
    }
  }
}

void ConvolvedSpikeTrain::UpdateMarkageVector(){
  if (size_ > 0){
    markage_[0] = 0;
    for (unsigned int i=1; i<size_; ++i){
      markage_[i] = (1 + markage_[i-1]) * exp_pos_[i-1] * exp_neg_[i];
    }
  }
}
