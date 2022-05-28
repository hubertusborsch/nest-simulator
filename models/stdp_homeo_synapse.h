/*
 *  stdp_homeo_synapse.h
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

#ifndef STDP_HOMEO_SYNAPSE_H
#define STDP_HOMEO_SYNAPSE_H

// C++ includes:
#include <cmath>
// #include <iostream>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

namespace nest
{

/** @BeginDocumentation
Name: stdp_homeo_synapse - Synapse type for spike-timing dependent
plasticity with homeostasis control.

Short description
+++++++++++++++++

Synapse type for spike-timing dependent plasticity and homeostasis control

Description
+++++++++++

stdp_homeo_synapse is a connector to create synapses with spike time
dependent plasticity and homeostasis controlled by the rate of dendritic
action potentials (dAPs) and firing rate of the postsynaptic neuron. 
A similar model (without firing rate homeostasis) is described in [1]. 

=======
Parameters
++++++++++

=========     =======  ======================================================
 tau_plus     ms       Time constant of STDP window, potentiation
                       (tau_minus defined in postsynaptic neuron)
 lambda_plus  real     Potentiation rate
 lambda_minus real     Depression rate
 lambda_h     real     Homeostasis rate for dAP rate homeostasis
 lambda_s     real     Homeostasis rate for firing rate homeostasis
 mu_plus      real     Weight dependence exponent, potentiation
 Wmax         real     Maximum allowed weight
=========     =======  ======================================================

Transmits
+++++++++

SpikeEvent

References
++++++++++

.. [1] Bouhadjar et al., (2021). Sequence learning, prediction, and replay in networks of spiking neurons.

See also
++++++++

tsodyks_synapse, static_synapse

EndUserDocs */

// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template

template < typename targetidentifierT >
class stdp_homeo_synapse : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  stdp_homeo_synapse();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  stdp_homeo_synapse( const stdp_homeo_synapse& ) = default;

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e, thread t, const CommonSynapseProperties& cp );


  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
  };

  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike_ - get_delay(), get_delay() );
  
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  double
  facilitate_( double w, double kplus )
  {
    double norm_w = ( w / Wmax_ ) + ( lambda_plus_ * std::pow( 1.0 - ( w / Wmax_ ), mu_plus_ ) * kplus);
    return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
  }

  // double
  // depress_( double w )
  // {
    //printf("# Depress #");
    // w = w - lambda_minus_ * Wmax_;
    // return w > Wmin_ ? w : Wmin_;
  // }

  double
  depress_( double w, double kminus )
  {
    double norm_w = ( w / Wmax_ ) - ( lambda_minus_ * std::pow( w / Wmax_, mu_minus_ ) * kminus );
    return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
  }

  double
  homeostasis_control_dAP( double weight, double z )
  {
    weight = weight + lambda_h_ * (zt_ - z) * Wmax_;
    return weight < Wmin_ ? Wmin_ : weight > Wmax_ ? Wmax_ : weight;
  }

  double
  homeostasis_control_soma( double weight, double s )
  {
    weight = weight + lambda_s_ * (st_ - s) * Wmax_;
    return weight < Wmin_ ? Wmin_ : weight > Wmax_ ? Wmax_ : weight;
  }

  // data members of each connection
  double weight_;
  double tau_plus_;
  double lambda_plus_;
  double lambda_minus_;
  double mu_plus_;
  double mu_minus_;
  double Wmax_;
  double Wmin_;
  double Kplus_;
  double zt_;
  double st_;
  double lambda_h_;
  double lambda_s_;
  double dt_min_;
  double dt_max_;

  double t_lastspike_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
stdp_homeo_synapse< targetidentifierT >::send( Event& e, thread t, const CommonSynapseProperties& )
{
  // synapse STDP depressing/facilitation dynamics
  double t_spike = e.get_stamp().get_ms();

  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  //bool reach_max_activity = target->get_reach_max_activity();

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< histentry >::iterator start;
  std::deque< histentry >::iterator finish;

  // For a new synapse, t_lastspike_ contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by ArchivingNode::register_stdp_connection(). See bug #218 for
  // details.
  target->get_history( t_lastspike_ - dendritic_delay,
    t_spike - dendritic_delay,
    &start,
    &finish );
  // facilitation due to post-synaptic spikes since last pre-synaptic spike
  double minus_dt;

  while ( start != finish )
  {
    minus_dt = t_lastspike_ - ( start->t_ + dendritic_delay );

    // get dendritic and somatic firing rates
    double z = start->dAP_trace_ ; 
    double s = start->spike_trace_;
    // std::cout << "this is z" << z << "\n";
    ++start;

    // get_history() should make sure that
    // start->t_ > t_lastspike - dendritic_delay, i.e. minus_dt < 0
    assert( minus_dt < -1.0 * kernel().connection_manager.get_stdp_eps() );
    if ( minus_dt > dt_max_ and  minus_dt < dt_min_ ){
    
        // Hebbian learning 
        // weight_ = facilitate_( weight_, Kplus_ * std::exp( minus_dt / tau_plus_ ));
        weight_ = facilitate_( weight_, 1.0 ); 

        // homeostasis control based on dAP firing rate
        weight_ = homeostasis_control_dAP( weight_, z );

        // homeostasis control based on somatic firing rate
        weight_ = homeostasis_control_soma( weight_, s );
 
    }

  }

  // depress 
  const double _K_value = target->get_K_value( t_spike - dendritic_delay );
  weight_ = depress_( weight_, _K_value ); 
  // weight_ = depress_( weight_, 1 ); 

  e.set_receiver( *target ); 
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  Kplus_ = Kplus_ * std::exp( ( t_lastspike_ - t_spike ) / tau_plus_ ) + 1.0;

  t_lastspike_ = t_spike;
}


template < typename targetidentifierT >
stdp_homeo_synapse< targetidentifierT >::stdp_homeo_synapse()
  : ConnectionBase()
  , weight_( 1.0 )
  , tau_plus_( 20.0 )
  , lambda_plus_( 0.01 )
  , lambda_minus_( 1.0 )
  , zt_( 1.0 )
  , st_( 10.0 )
  , lambda_h_( 0.01 )
  , lambda_s_( 0.01 )
  , mu_plus_( 0.005 )
  , mu_minus_( 1.0 )
  , Wmax_( 100.0 )
  , Wmin_( 0.0 )
  , dt_min_( 0.0 )
  , dt_max_( -50.0 )
  , Kplus_( 0.0 )
  , t_lastspike_( -100.0 )
{
}

template < typename targetidentifierT >
void
stdp_homeo_synapse< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::tau_plus, tau_plus_ );
  def< double >( d, names::lambda_plus, lambda_plus_ );
  def< double >( d, names::lambda_minus, lambda_minus_ );
  def< double >( d, names::zt, zt_ );
  def< double >( d, names::st, st_ );
  def< double >( d, names::lambda_h, lambda_h_ );
  def< double >( d, names::lambda_s, lambda_s_ );
  def< double >( d, names::mu_plus, mu_plus_ );
  def< double >( d, names::mu_minus, mu_minus_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::dt_min, dt_min_ );
  def< double >( d, names::dt_max, dt_max_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
stdp_homeo_synapse< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::tau_plus, tau_plus_ );
  updateValue< double >( d, names::lambda_plus, lambda_plus_ );
  updateValue< double >( d, names::lambda_minus, lambda_minus_ );
  updateValue< double >( d, names::zt, zt_ );
  updateValue< double >( d, names::st, st_ );
  updateValue< double >( d, names::lambda_h, lambda_h_ );
  updateValue< double >( d, names::lambda_s, lambda_s_ );
  updateValue< double >( d, names::mu_plus, mu_plus_ );
  updateValue< double >( d, names::mu_minus, mu_minus_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::dt_min, dt_min_ );
  updateValue< double >( d, names::dt_max, dt_max_ );

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }
}

} // of namespace nest

#endif // of #ifndef stdp_homeo_synapse_H
