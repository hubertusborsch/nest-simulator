/*
 *  stdsp_connection.h
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

#ifndef STDSP_SYNAPSE_H
#define STDSP_SYNAPSE_H

// C++ includes:
#include <cmath>

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
//TODO: change documentation of stdsp to stdsp ("spike-timing dependent structural plasticity")

/** @BeginDocumentation
Name: stdsp_synapse - Synapse type for spike-timing dependent structural plasticity.

Description:

stdsp_synapse is a connector to create synapses with spike-timing-dependent structural plasticity (as defined in [ref]).

Examples:

//TODO: change the following with exmples from stdsp
multiplicative STDP [2]  mu_plus = mu_minus = 1.0
additive STDP       [3]  mu_plus = mu_minus = 0.0
Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0

//TODO: change the following with stdsp parameters
Parameters:

tau_plus   double - Time constant of STDP window, potentiation in ms
                    (tau_minus defined in post-synaptic neuron)
lambda_plus     double - Step size
lambda_minus      double - Asymmetry parameter (scales depressing increments as
                    lambda_minus*lambda_plus)
mu_plus    double - Weight dependence exponent, potentiation
mu_minus   double - Weight dependence exponent, depression
Wmax       double - Maximum allowed weight

Transmits: SpikeEvent

References:

\\TODO add stdsp references

[1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
    Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

[2] Rubin, J., Lee, D. and Sompolinsky, H. (2001). Equilibrium
    properties of temporally asymmetric Hebbian plasticity, PRL
    86,364-367

[3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
    Hebbian learning through spike-timing-dependent synaptic
    plasticity,Nature Neuroscience 3:9,919--926

[4] van Rossum, M. C. W., Bi, G-Q and Turrigiano, G. G. (2000).
    Stable Hebbian learning from spike timing-dependent
    plasticity, Journal of Neuroscience, 20:23,8812--8821

FirstVersion: March 2006

Author: Moritz Helias, Abigail Morrison

Adapted by: Younes Bouhadjar

SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/
// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class stdsp_synapse : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  stdsp_synapse();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  stdsp_synapse( const stdsp_synapse& ) = default;

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

    t.register_stdp_connection( t_lastspike_ - get_delay(), get_delay());
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

  void set_permanence(double perm)
  {
    permanence_ = perm;
  }    

private:
  double
  facilitate_( double perm )
  {
    // printf("# Facilitate #");
    perm = perm + Delta_plus_;
    return perm < Pmax_ ? perm : Pmax_;
  }

  double
  facilitate_exp_( double perm, double kplus)
  {
    double mu = 0; 
    double norm_perm = ( perm / Pmax_ )
    + ( lambda_plus_ * std::pow( 1.0 - ( perm / Pmax_ ), mu ) * kplus );
    return norm_perm < 1.0 ? norm_perm * Pmax_ : Pmax_;
  }

  double
  depress_( double perm )
  {
    //printf("# Depress #");
    perm = perm - lambda_minus_ * Pmax_;
    return perm > init_perm_ ? perm : init_perm_;
  }
  // not used 
  double
  depress_exp_( double perm, double kminus )
  {
    double norm_perm = ( perm / Pmax_ )
      - ( lambda_minus_ * lambda_plus_ * std::pow( perm / Pmax_, mu_minus_ ) * kminus );
    return norm_perm > 0.0 ? norm_perm * Pmax_ : 0.0;
  }

  //RecordingDevice device_;

  // data members of each connection
  double weight_;
  double permanence_;
  double tau_plus_;
  double lambda_plus_;
  double lambda_minus_;
  double mu_plus_;
  double mu_minus_;
  double Wmax_;
  double Pmax_;
  double Kplus_;
  double Delta_plus_;
  double Delta_minus_;
  double lambda_h_;

  // remove t_mean_ and t_var_ 
  double t_mean_;
  double t_var_;

  double th_perm_;
  double t_lastspike_;
  double zt_;
  double max_dt_ = -100.;
  double min_dt_ = -4.;
  double init_perm_ = permanence_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param cp Common properties object, containing the stdsp parameters.
 */
template < typename targetidentifierT >
inline void
stdsp_synapse< targetidentifierT >::send( Event& e,
  thread t,
  const CommonSynapseProperties& )
{
  // synapse STDP depressing/facilitation dynamics
  double t_spike = e.get_stamp().get_ms();

  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< histentry >::iterator start;
  std::deque< histentry >::iterator finish;

  // For a new synapse, t_lastspike_ contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by Archiving_Node::register_stdsp_connection(). See bug #218 for
  // details.
  target->get_history( t_lastspike_ - dendritic_delay,
    t_spike - dendritic_delay,
    &start,
    &finish );
  
  // facilitation due to post-synaptic spikes since last pre-synaptic spike
  double minus_dt; 
  
  float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  
  // get the dendritic firing rate
  double z = target->get_dendritic_firing_rate();

  while ( start != finish )
  {
      minus_dt = t_lastspike_ - ( start->t_ + dendritic_delay ); 
      // printf("\n last spike %lf, start %lf, minus_dt %f, t %f", t_lastspike_ , start->t_, minus_dt, t_spike);
      ++start;

      // get_history() should make sure that
      // start->t_ > t_lastspike - dendritic_delay, i.e. minus_dt < 0
      assert( minus_dt < -1.0 * kernel().connection_manager.get_stdp_eps() );
      if ( minus_dt > max_dt_ and minus_dt < min_dt_ ){
          
          // hebbian learning     
          permanence_ = facilitate_exp_( permanence_, Kplus_ * std::exp( minus_dt / tau_plus_ ));
          
          // homeostasis control
          permanence_ += lambda_h_* (zt_ - z) * Pmax_; 
      }
   }
   
   // depression due to new pre-synaptic spike
   permanence_ = depress_( permanence_ );

  // update weight
  if ( permanence_ > th_perm_ )
  {
    weight_ = Wmax_;
    //permanence_ = depress_( permanence_ );
  }  
  else
  {
    weight_ = 0.001;
  }
      
  e.set_receiver( *target );
  e.set_weight( weight_ );
  set_permanence( permanence_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  Kplus_ = Kplus_ * std::exp( ( t_lastspike_ - t_spike ) / tau_plus_ ) + 1.0;

  t_lastspike_ = t_spike;
}


template < typename targetidentifierT >
stdsp_synapse< targetidentifierT >::stdsp_synapse()
  : ConnectionBase()
  , weight_( 0.1 )
  , permanence_(1.0)
  , tau_plus_( 80.0 )
  , lambda_plus_( 0.01 )
  , lambda_minus_( 1.0 )
  , mu_plus_( 1.0 )
  , mu_minus_( 1.0 )
  , Wmax_( 100.0 )
  , Pmax_( 100.0 )  
  , th_perm_( 10.0 )  
  , Kplus_( 0.0 ) 
  , Delta_plus_( 0.1 ) 
  , Delta_minus_( 0.0 )
  , t_lastspike_( 0.0 )
  , lambda_h_( 1.0 )
  , zt_( 1.0 )
  , t_mean_( 30.0 )
  , t_var_( 5.0 )
{
}

template < typename targetidentifierT >
void
stdsp_synapse< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ ); 
  def< double >( d, names::permanence, permanence_);
  def< double >( d, names::th_perm, th_perm_);
  def< double >( d, names::tau_plus, tau_plus_ );
  def< double >( d, names::lambda_plus, lambda_plus_ );
  def< double >( d, names::lambda_minus, lambda_minus_ );
  def< double >( d, names::mu_plus, mu_plus_ );
  def< double >( d, names::mu_minus, mu_minus_ );
  def< double >( d, names::Wmax, Wmax_ ); 
  def< double >( d, names::Pmax, Pmax_ ); 
  def< double >( d, names::lambda_h, lambda_h_ ); 
  def< double >( d, names::zt, zt_ ); 
  def< double >( d, names::t_mean, t_mean_ ); 
  def< double >( d, names::t_var, t_var_ ); 
  def< double >( d, names::Delta_plus, Delta_plus_ );
  def< double >( d, names::Delta_minus, Delta_minus_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
stdsp_synapse< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ ); 
  updateValue< double >( d, names::permanence, permanence_ );
  updateValue< double >( d, names::th_perm, th_perm_ ); 
  updateValue< double >( d, names::tau_plus, tau_plus_ );
  updateValue< double >( d, names::lambda_plus, lambda_plus_ );
  updateValue< double >( d, names::lambda_minus, lambda_minus_ );
  updateValue< double >( d, names::mu_plus, mu_plus_ );
  updateValue< double >( d, names::mu_minus, mu_minus_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::Pmax, Pmax_ );
  updateValue< double >( d, names::lambda_h, lambda_h_ );
  updateValue< double >( d, names::zt, zt_ );
  updateValue< double >( d, names::Delta_plus, Delta_plus_ );
  updateValue< double >( d, names::Delta_minus, Delta_minus_ );
  updateValue< double >( d, names::t_mean, t_mean_ );
  updateValue< double >( d, names::t_var, t_var_ );

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }
}

} // of namespace nest

#endif // of #ifndef STDSP_SYNAPSE_H
