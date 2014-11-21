/*! \file van_rossum_multiunit.hpp
 *
 *  C++ implementation of a kernel-based algorithm for calculating
 *  multi-unit Van Rossum spike train distances. See Houghton and
 *  Kreuz, On the efficient calculation of van Rossum
 *  distances. Network: Computation in Neural Systems, 2012, 23,
 *  48-58.
 *
 *  Original implementation of multiunit distance algorithm by Conor
 *  Houghton (2012). Later work by Eugenio Piasini (2013, 2014).
 *
 */
#ifndef VAN_ROSSUM_H
#define VAN_ROSSUM_H

#include "convolved_spike_train.hpp"

#include<vector>


/*!
 * Square (all-to-all) distance matrix for a set of multi-unit spike
 * trains.
 *
 * 
 * \param[out] d_matrix The \a trains.size() x \a trains.size() matrix
 * where the results will be written.
 *
 * \param[in] trains The set of multi unit spike trains.
 *
 * \param[in] tau Time scale for the exponential kernel. There is a
 * limit to how small \a tau can be compared to the absolute value of
 * the spike times. An exception will be raised if this limit is
 * exceeded; its value is system-dependent, but as a rule of thumb \a
 * tau and the spike times should be within 4 orders of magnitude of
 * each other.
 *
 * \param[in] c Cosine of the multi-unit mixing angle. \a c=0 corresponds
 * to labelled-line code, \a c=1 to summed-population code.
 */
void d_exp_markage(double **d_matrix,
		   std::vector< std::vector<ConvolvedSpikeTrain> > & trains,
		   double c);


/*!
 * Rectangular (bipartite) distance matrix for a set of multi-unit
 * spike trains (observations).
 *
 * 
 * \param[out] d_matrix The \a trains1.size() x \a trains2.size()
 * matrix where the results will be written.
 *
 * \param[in] trains1,trains2 The sets of multi unit spike trains.
 *
 * \param[in] tau Time scale for the exponential kernel. There is a
 * limit to how small \a tau can be compared to the absolute value of
 * the spike times. An exception will be raised if this limit is
 * exceeded; its value is system-dependent, but as a rule of thumb \a
 * tau and the spike times should be within 4 orders of magnitude of
 * each other.
 *
 * \param[in] c Cosine of the multi-unit mixing angle. \a c=0 corresponds
 * to labelled-line code, \a c=1 to summed-population code.
 */
void d_exp_markage_rect(double **d_matrix,
			std::vector< std::vector<ConvolvedSpikeTrain> > & trains1,
			std::vector< std::vector<ConvolvedSpikeTrain> > & trains2,
			double c);

#endif
