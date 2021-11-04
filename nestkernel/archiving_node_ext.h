/*
 *  archiving_node.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ARCHIVING_NODE_EXT_H
#define ARCHIVING_NODE_EXT_H

// C++ includes:
#include <algorithm>
#include <deque>

// Includes from nestkernel:
#include "histentry.h"
#include "nest_time.h"
#include "nest_types.h"
#include "node.h"
#include "structural_plasticity_node.h"
#include "archiving_node.h"

// Includes from sli:
#include "dictdatum.h"

#define DEBUG_ARCHIVER 1

namespace nest
{

/**
 * A node which archives spike history for the purposes of spike-timing
 * dependent plasticity (STDP)
 */
class ArchivingNodeExt : public ArchivingNode
{
protected:
  /**
   * \fn void set_spiketime(Time const & t_sp, double offset)
   * record spike history
   */
  void set_spiketime( Time const& t_sp, double offset = 0.0 );

//private:
//  // sum exp(-(t-ti)/tau_minus)
//  double Kminus_;
//
//  // sum exp(-(t-ti)/tau_minus_triplet)
//  double Kminus_triplet_;
//  double dAP_trace_;
//  double spike_trace_;
//
//  double tau_minus_;
//  double tau_minus_inv_;
//
//  // time constant for triplet low pass filtering of "post" spike train
//  double tau_minus_triplet_;
//  double tau_minus_triplet_inv_;
//
//  double max_delay_;
//  double trace_;
//
//  double last_spike_;
//
//  // spiking history needed by stdp synapses
//  std::deque< histentry > history_;
//
};

} // of namespace
#endif
