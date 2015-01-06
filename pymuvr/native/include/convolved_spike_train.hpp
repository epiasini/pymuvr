/*!
 * \file convolved_spike_train.hpp
 * \author Eugenio Piasini
 * \date 2014
 * \copyright GPLv3+
 */

#ifndef __CONVOLVED_SPIKE_TRAIN_HPP_INCLUDED__
#define __CONVOLVED_SPIKE_TRAIN_HPP_INCLUDED__

#include<vector>

/* Visual Studio does not have isfinite; use _finite() instead */
#ifdef _MSC_VER
#include <float.h>
bool isfinite(long double x);
#endif

/*!
 * A class representing a single-unit spike train for which several
 * helper vectors relating to a Van Rossum-like exponential
 * convolution have been calculated.
 *
 */
class ConvolvedSpikeTrain {
  void UpdateExponentialVectors(); /*!< Compute and store exponential vectors. */
  void UpdateMarkageVector(); /*!< Compute and store markage vector. */
  void UpdateSquareNorm(); /*< Compute and store square norm of spike train. */
public:
  /* Data assigned upon instantiation. For ease of access, it's best
   * to keep all data members public. */
  std::vector<double> spikes; /*!< Vector of spike times */
  double tau; /*!< Time scale for the exponential kernel. */

  /* Data that gets computed internally */
  unsigned int size; /*!< Number of spikes. */
  std::vector<long double> exp_pos; /*!< Vector containing the value exp(t/tau) for each spike time. */
  std::vector<long double> exp_neg; /*!< Vector containing the value exp(-t/tau) for each spike time. */
  std::vector<double> markage; /*!< Markage vector */
  double square_norm; /*!< Square norm of the vector. */
  
  /*!
   * Default constructor.
   */
  ConvolvedSpikeTrain();
  
  /*!
   * Constructor accepting a spike train.
   *
   * \param[in] spikes A \a vector of spike times.
   *
   * \param[in] tau Time scale for the exponential kernel. There is a
   * limit to how small \a tau can be compared to the absolute value
   * of the spike times. An exception will be raised if this limit is
   * exceeded; its value is system-dependent, but as a rule of thumb
   * \a tau and the spike times should be within 4 orders of magnitude
   * of each other.
   *
   */
  ConvolvedSpikeTrain(std::vector<double> spikes, double tau);
  
  /*!
   * Update values of helper vectors using a new value for the time
   * constant of the exponential kernel.
   */
  void UpdateConvolution(double tau);  
};

#endif
