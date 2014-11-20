#include "convolved_spike_train.hpp"

#include<stdexcept>
#include<vector>
#include<numeric>
#include<cmath>

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

using std::exp; using std::isfinite; // from cmath
using std::vector; // from vector
using std::accumulate; // from numeric
using std::overflow_error; // from stdexcept

ConvolvedSpikeTrain::ConvolvedSpikeTrain(): spikes(), tau(), size(), exp_pos(), exp_neg(), markage(), square_norm(){}
    
ConvolvedSpikeTrain::ConvolvedSpikeTrain(vector<double> spikes, double tau) : spikes(spikes), tau(tau), size(spikes.size()), exp_pos(size), exp_neg(size), markage(size){
  UpdateConvolution(tau);
}

void ConvolvedSpikeTrain::UpdateConvolution(double new_tau){
  tau = new_tau;
  UpdateExponentialVectors();
  UpdateMarkageVector();
  UpdateSquareNorm();
}

void ConvolvedSpikeTrain::UpdateExponentialVectors(){
  if (size > 0){
    for (unsigned int i=0; i<size; ++i){
      exp_pos[i] = exp((long double)(spikes[i]/tau));
      exp_neg[i] = exp((long double)(-spikes[i]/tau));
    }

    /* Check if any over/underflow occurred while calculating the
     * exponentials */
    if (tau!=0 && (exp_neg.front()==0 || !isfinite(exp_pos.back()))){
      throw overflow_error("tau is too small compared to the spike times. Please use a larger value for tau, or shorter spike trains.");
    }
  }
}

void ConvolvedSpikeTrain::UpdateMarkageVector(){
  if (size > 0){
      markage[0] = 0;
      for (unsigned int i=1; i<size; ++i){
	if (tau != 0){
	  markage[i] = (1 + markage[i-1]) * exp_pos[i-1] * exp_neg[i];
	} else {
	  /* If tau is zero, we have to work around the Inf * 0
	   * expression that would result from multiplyting exp_pos
	   * and exp_neg */
	  markage[i] = 0;
	}
      }
  }
}

void ConvolvedSpikeTrain::UpdateSquareNorm(){
  /* Note that this is still valid even if tau is zero or infinity */
  square_norm = size + 2 * accumulate(markage.begin(), markage.end(), 0.0);
}
