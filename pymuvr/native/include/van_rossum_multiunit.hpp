/*! \file van_rossum_multiunit.hpp
 *
 *  C++ implementation of a kernel-based algorithm for calculating
 *  multi-unit Van Rossum spike train distances. See Houghton and
 *  Kreuz, On the efficient calculation of van Rossum
 *  distances. Network: Computation in Neural Systems, 2012, 23,
 *  48-58.
 *
 *  Credit to Conor Houghton (2012) for the first implementation of
 *  the algorithm published in the paper above. Later work by Eugenio
 *  Piasini (2013, 2014).
 *
 *  \date 2012-2014
 *  \copyright GPLv3+
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
 * \param[in] trains The set of multi unit spike trains.
 *
 * \param[in] c Cosine of the multi-unit mixing angle. \a c=0 corresponds
 * to labelled-line code, \a c=1 to summed-population code.
 * 
 * \param[out] d_matrix The \a trains.size() x \a trains.size() matrix
 * where the results will be written.
 * */
void distance(std::vector< std::vector<ConvolvedSpikeTrain> > & trains,
	      double c,
	      double **d_matrix);


/*!
 * Rectangular (bipartite) distance matrix for a set of multi-unit
 * spike trains (observations).
 *
 * 
 * \param[in] trains1,trains2 The sets of multi unit spike trains.
 *
 * \param[in] c Cosine of the multi-unit mixing angle. \a c=0 corresponds
 * to labelled-line code, \a c=1 to summed-population code.
 *
 * \param[out] d_matrix The \a trains1.size() x \a trains2.size()
 * matrix where the results will be written.
 */
void distance(std::vector< std::vector<ConvolvedSpikeTrain> > & trains1,
	      std::vector< std::vector<ConvolvedSpikeTrain> > & trains2,
	      double c,
	      double **d_matrix);


/*!
 * Square (all-to-all) inner product matrix for 
 */

#endif
