#ifndef __CONVOLVED_SPIKE_TRAIN_HPP_INCLUDED__
#define __CONVOLVED_SPIKE_TRAIN_HPP_INCLUDED__

#include<vector>

/*!
 * A class representing a single-unit spike train for which several
 * helper vectors relating to a Van Rossum-like exponential
 * convolution have been calculated.
 *
 */
class ConvolvedSpikeTrain {
  /* Data assigned upon instantiation */
  std::vector<double> spikes_;
  double tau_;
  /* Data that gets computed internally */
  unsigned int size_;
  std::vector<long double> exp_pos_;
  std::vector<long double> exp_neg_;
  std::vector<double> markage_;
  double square_norm_;
  
  void UpdateExponentialVectors(); /*!< Compute and store exponential vectors */
  void UpdateMarkageVector(); /*!< Compute and store markage vector */
  void UpdateSquareNorm(); /*< Compute and store square norm of spike train */
public:
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
  
  /* Attribute getters */
  unsigned int size() const { return size_; }
  double tau() const { return tau_; }
  double square_norm() const {return square_norm_; }
  std::vector<double> spikes() const { return spikes_; }
  std::vector<long double> exp_pos() const { return exp_pos_; }
  std::vector<long double> exp_neg() const { return exp_neg_; }
  std::vector<double> markage() const { return markage_; }

  /*!
   * Update values of helper vectors using a new value for the time
   * constant of the exponential kernel.
   */
  void UpdateConvolution(double tau);  
};

#endif