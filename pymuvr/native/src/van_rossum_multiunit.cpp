#include "van_rossum_multiunit.hpp"

#include "convolved_spike_train.hpp"

#include<stdexcept>
#include<vector>
#include<cmath>
/* Check if we're compiling with Visual Studio. */
#ifdef _MSC_VER
#include <float.h>
#endif

/* Visual Studio does not have isfinite; use _finite() instead */
#ifdef _MSC_VER
bool isfinite(long double x) {
    return _finite(x);
}
#endif

using std::sqrt; // from cmath 
using std::vector; // from vector
using std::invalid_argument; // from stdexcept

/*! 
 * Compute the inner product between two single-unit convolved spike
 * trains.
 */
double inner_product(ConvolvedSpikeTrain & u, ConvolvedSpikeTrain & v);

/*! 
 * Compute the inner product of two convolved single unit spike trains
 * for the limit case where tau=0.
 *
 * In practice, this counts the number of coincident spike pairs
 * between train_a and train_b.
 *
 * Note that this assumes that neither train_a nor train_b cointain
 * repeating spikes.
 */
int inner_product_zero_tau(ConvolvedSpikeTrain & u, ConvolvedSpikeTrain & v);


void distance(vector< vector<ConvolvedSpikeTrain> > & trains,
	      double c,
	      double **d_matrix)
{
  unsigned int big_n = trains.size(); // number of observations
  unsigned int big_p = trains.front().size(); // number of cells per observation

  /* make sure the diagonal of d_matrix is set to zero */
  for(unsigned int n=0; n<big_n; n++) {
    d_matrix[n][n] = 0;
  }

  /* compute helper vector */
  vector<double> squares(big_n, 0.0);
  for(unsigned int n=0; n<big_n; ++n) {
    /* same-cell part */
    for(unsigned int p=0; p<big_p; ++p) {
      squares[n] += trains[n][p].square_norm;
    }
    /* cross-cell part */
    for(unsigned int p=0; p<big_p-1; ++p) {
      for(unsigned int q=p+1; q<big_p; ++q) {
	squares[n] += 2 * c * inner_product(trains[n][p], trains[n][q]);
      }
    }
  }

  for(unsigned int n=0; n<big_n-1; ++n) {
    for(unsigned int m=n+1; m<big_n; ++m) {
      double d = squares[n] + squares[m];
      /* same-cell part */
      for(unsigned int p=0; p<big_p; ++p) {
	d -= 2 * inner_product(trains[n][p], trains[m][p]);
      }
      /* cross-cell part */
      for(unsigned int p=0; p<big_p-1; ++p) {
	for(unsigned int q=p+1; q<big_p; ++q) {
	  d -= 2 * c * inner_product(trains[n][p], trains[m][q]);
	  d -= 2 * c * inner_product(trains[m][p], trains[n][q]);
	}
      }
      /* d is actually the square distance, so we need to take the sqrt */
      if (d>0) {
	d_matrix[n][m] = d_matrix[m][n] = sqrt(d);
      } else {
	/* this could happen due to floating point inaccuracies */
	d_matrix[n][m] = d_matrix[m][n] = 0;
      }
    }
  }
}


void distance(vector< vector<ConvolvedSpikeTrain> > & trains1,
	      vector< vector<ConvolvedSpikeTrain> > & trains2,
	      double c,
	      double **d_matrix)
{
  unsigned int big_n=trains1.size(); // number of observations in first set
  unsigned int big_m=trains2.size(); // number of observations in second set
  unsigned int big_p=trains1.front().size(); // number of cells

  if (trains2.front().size() != big_p) {
    throw invalid_argument("trying to compare two observations with a different number of cells.");
  }

  for(unsigned int n=0; n<big_n; ++n) {
    for(unsigned int m=0; m<big_m; ++m) {
      /* compute same-cell ("labelled-line") term */
      double d_same = 0;
      for(unsigned int p=0; p<big_p; ++p) {
	d_same += trains1[n][p].square_norm;
	d_same += trains2[m][p].square_norm;
	d_same -= 2 * inner_product(trains1[n][p], trains2[m][p]);
      }
      
      /* compute cross-cell ("summed-population") term */
      double d_cross = 0;
      for(unsigned int p=0; p<big_p-1; ++p) {
	for(unsigned int q=p+1; q<big_p; ++q) {
	  /* same-observation, cross-cell */
	  d_cross += inner_product(trains1[n][p], trains1[n][q]);
	  d_cross += inner_product(trains2[m][p], trains2[m][q]);
	  /* cross-observation, cross-cell */
	  d_cross -= inner_product(trains1[n][p], trains2[m][q]);
	  d_cross -= inner_product(trains2[m][p], trains1[n][q]);
	}
      }
      
      /* Sum same-unit and cross-unit with the appropriate weighting
	 to give the desired interpolation between labelled-line and
	 summed-population (square) distance. */
      double d = d_same + 2 * c * d_cross;
      if (d>0) {
	d_matrix[n][m] = sqrt(d);
      } else {
	/* this could happen due to floating point inaccuracies */
	d_matrix[n][m] = 0;
      }
    }
  }
}


int inner_product_zero_tau(ConvolvedSpikeTrain & u,
			   ConvolvedSpikeTrain & v)
{
  if(u.size==0 || v.size==0) {
    return 0;
  }

  int result = 0;
  int j = u.size - 1;

  for (int i=v.size-1; i>=0; --i) {
    /* Look for the index of largest spike time in u which is smaller
	or equal than the ith spike time of v. Leave the index set to
	zero if you don't find any. */
    while(j>0 && u.spikes[j]>v.spikes[i]) {
      j--;
    }
    /* If the spike selected in u coincides with the one we're
       considering in v, add 1 to the result. */
    if (u.spikes[j]==v.spikes[i]) {
      result+=1;
    }
  }

  return result;
}


double inner_product(ConvolvedSpikeTrain & u, ConvolvedSpikeTrain & v)
{ 
  if(u.size==0 || v.size==0) {
    return 0;
  }

  if (u.tau==0) {
    if (u.tau!=v.tau) {
      throw invalid_argument("trying to compute the inner product of a spike train convolved with tau=0 and one with tau!=0");
    } else {
      return inner_product_zero_tau(u, v);
    }
  }
  double result = 0;

  // we'll use j as an index to iterate backwards on the spikes in u
  int j = u.size - 1;

  for(int i=v.size-1; i>=0; --i) {
    /* Look for the index of largest spike time in u which is smaller
       or equal than the ith spike time of v. Exit from the loop if
       you don't find any. Note that this, in general, is not the same
       thing as calculating J(i) in the paper (Houghton and Kreuz
       2012), because we allow for u[j] to be equal to v[i]. */
    while(j>=0 && u.spikes[j]>v.spikes[i]) {
      j--;
    }
    if(j<0) {
      break;
    }
    /* If the spike selected in u coincides with the one we're
       considering in v, add 1 to the quantity being
       calculated. Adjust the 'j' index to point at the previous
       spike, which now is guaranteed to correspond to the J(i) index
       in the paper mentioned above; unless it does not exist, in
       which case we exit from the loop. Note that this passage is not
       mirrored in the 'symmetric' loop below, as we only need to
       count each pair of coincident spikes once. */
    if (u.spikes[j]==v.spikes[i]) {
      result+=1;
      j--;
    }
    if(j<0) {
      break;
    }
    /* Compute the ith term of the main sum we're calculating, using
       the precomputed exponentials and the markage vector. */
    result += u.exp_pos[j] * v.exp_neg[i] * (u.markage[j] + 1);	
  }

  // reset j; we'll now use it as an index to iterate backwards on the spikes in v
  j = v.size - 1;

  for(int i=u.size-1; i>=0; --i) {
    /* Unlike above, calculate directly the index corresponding to the
       largest spike time in v which is _strictly_ smaller than the
       ith spike time of u. Exit from the loop if you don't find
       any. */
    while(j>=0 && v.spikes[j]>=u.spikes[i]) {
      j--;
    }
    if(j<0) {
      break;
    }
    /* Compute the ith term of the main sum we're calculating, using
       the precomputed exponentials and the markage vector. */
    result += v.exp_pos[j] * u.exp_neg[i] * (v.markage[j] + 1);
  }
  return result;
}
