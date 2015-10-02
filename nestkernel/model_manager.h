/*
 *  model_manager.h
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

#ifndef MODEL_MANAGER_H
#define MODEL_MANAGER_H

#include "manager_interface.h"
#include "nest_types.h"
#include "nest_time.h"
#include "nest_timeconverter.h"
#include "node.h"
#include "genericmodel.h"
#include "connector_model.h"

#include "dictutils.h"
#include "network.h"
#include "model.h"
#include "genericmodel.h"

namespace nest
{
class ModelManager : public ManagerInterface
{
public:
  ModelManager();
  ~ModelManager();

  /**
   *
   */
  void init();

  /**
   *
   */
  void reset();


  // TODO: register subnet, proxynode and siblingcontainer in the
  // constructor and provide three getters for them. See L130 and
  // following in network.cpp. Also make sure that the ModelManager is
  // initialized before the NodeManager is.

  /**
   * Resize the structures for the Connector objects if necessary.
   * This function should be called after number of threads, min_delay, max_delay,
   * and time representation have been changed in the scheduler.
   * The TimeConverter is used to convert times from the old to the new representation.
   * It is also forwarding the calibration
   * request to all ConnectorModel objects.
   */
  void calibrate( const TimeConverter& );


  /**
   *
   */
  void set_status( const DictionaryDatum& );

  /**
   *
   */
  void get_status( DictionaryDatum& );

  /**
   *
   */
  Model* get_subnet_model();

  /**
   *
   */
  Model* get_siblingcontainer_model();

  /**
   *
   */
  Node* get_proxy_node( thread, index );

  /**
   * Return pointer to protoype for given synapse id.
   * @throws UnknownSynapseType
   */

  //  TODO: make the return type const, after the increment of
  //  num_connections and the min_ and max_delay setting in
  //  ConnectorBase was moved out to the ConnectionManager
  ConnectorModel& get_synapse_prototype( synindex syn_id, thread t = 0 );

  /**
   * Register a node-model prototype.
   * This function must be called exactly once for each model class to make
   * it known in the simulator. The natural place for a call to this function
   * is in a *module.cpp file.
   * @param name of the new node model.
   * @param private_model if true, don't add model to modeldict.
   * @return ID of the new model object.
   * @see register_private_prototype_model, register_preconf_node_model, register_prototype_connection
   */
  template < class ModelT >
  index
  register_node_model( const Name& name, bool private_model = false );
  
  /**
   * Register a pre-configured model prototype with the network.
   * This function must be called exactly once for each model class to make
   * it known to the network. The natural place for a call to this function
   * is in a *module.cpp file.
   *
   * Pre-configured models are models based on the same class, as
   * another model, but have different parameter settings; e.g.,
   * voltmeter is a pre-configured multimeter.
   *
   * @param name of the new node model.
   * @param private_model if true, don't add model to modeldict.
   * @param dictionary to use to pre-configure model
   * @return ID of the new model object.
   * @see register_private_prototype_model, register_node_model, register_prototype_connection
   */
  template < class ModelT >
  index
    register_preconf_node_model( const Name& name, DictionaryDatum& conf, bool private_model = false );

  /**
   * Copy an existing model and register it as a new model.
   * This function allows users to create their own, cloned models.
   * @param old_name name of existing model.
   * @param new_name name of new model.
   * @param params default parameters of new model. 
   * @return model ID of new Model object.
   * @see copy_node_model_, copy_synapse_model_
   */
  index copy_model( Name old_name, Name new_name, DictionaryDatum params );

  /**
   * Set the default parameters of a model.
   * @param name of model.
   * @param params default parameters to be set. 
   * @see set_node_defaults_, set_synapse_defaults_
   */
  void
  set_model_defaults( Name name, DictionaryDatum params );

  /**
   * Register a synapse type.
   * @param cm ConnectorModel to be registered.
   * @return an ID for the synapse prototype.
   */
  synindex register_synapse_prototype( ConnectorModel* cf );


  /**
   * Register a synape with default Connector and without any common properties.
   */
  template < class ConnectionT >
  synindex
  register_connection_model( const std::string& name );


  /**
   * Copy an existing synapse type.
   * @see copy_model(), ModelManager::copy_synapse_prototype()
   * @param old_id ID of synapse model to copy.
   * @param new_name name of new synapse model.
   * @return ID of new synapse model.
   */
  synindex copy_synapse_prototype( synindex old_id, Name new_name );

  /**
   * @return The model id of a given model name
   */
  int get_model_id( const Name ) const;

  /**
   * @return The Model of a given model ID
   */
  Model* get_model( index ) const;

  DictionaryDatum get_connector_defaults( synindex syn_id ) const;
  void set_connector_defaults( synindex syn_id, const DictionaryDatum& d );

  /**
   * Check, if there are instances of a given model.
   * @param i the index of the model to check for.
   * @return True, if model is instantiated at least once.
   */
  bool is_model_in_use( index i );

  /**
   * Checks, whether connections of the given type were created
   */
  bool synapse_prototype_in_use( synindex syn_id );

  /**
   * Asserts validity of synapse index, otherwise throws exception.
   * @throws UnknownSynapseType
   */
  void assert_valid_syn_id( synindex syn_id, thread t = 0 ) const;

  /**
   * @return Reference to the model dictionary
   */
  const DictionaryDatum& get_modeldict();

  /**
   * @return Reference to the synapse dictionary
   */
  const DictionaryDatum& get_synapsedict() const;

  /**
   * Does the network contain copies of models created using CopyModel?
   */
  bool has_user_models() const;

  bool has_user_prototypes() const;

  bool are_model_defaults_modified() const;

  const std::vector< ConnectorModel* >& get_prototypes( const thread t ) const;

  size_t get_num_node_models() const;

  size_t get_num_synapse_prototypes() const;

private:
  /**
   * The list of clean node models. The first component of the pair is a
   * pointer to the actual Model, the second is a flag indicating if
   * the model is private. Private models are not entered into the
   * modeldict.
   */
  std::vector< std::pair< Model*, bool > > pristine_models_;

  std::vector< Model* > models_;  //!< List of available models


  /**
   * The list of clean synapse models. The first component of the pair is a
   * pointer to the actual Model, the second is a flag indicating if
   * the model is private. Private models are not entered into the
   * modeldict.
   */
  std::vector< ConnectorModel* > pristine_prototypes_; //!< The list of clean synapse prototypes
  std::vector< std::vector< ConnectorModel* > > prototypes_; //!< The list of available synapse
                                                             //!< prototypes: first dimension one
                                                             //!< entry per thread, second dimension
                                                             //!< for each synapse type


  /* BeginDocumentation
     Name: modeldict - dictionary containing all devices and models of NEST
     Description:
     'modeldict info' shows the contents of the dictionary
     SeeAlso: info, Device, RecordingDevice, iaf_neuron, subnet
  */
  DictionaryDatum modeldict_;    //!< DictionaryDatum of all models

  /* BeginDocumentation
     Name: synapsedict - DictionaryDatum containing all synapse models.
     Description:
     'synapsedict info' shows the contents of the dictionary
     FirstVersion: October 2005
     Author: Jochen Martin Eppler
     SeeAlso: info
  */
  DictionaryDatum synapsedict_;  //!< DictionaryDatum of all synapse models

  Model* subnet_model_;
  Model* siblingcontainer_model_;
  Model* proxynode_model_;

  std::vector< std::vector< Node* > >
    proxy_nodes_; //!< Placeholders for remote nodes, one per thread

  bool model_defaults_modified_;  //!< True if any model defaults have been modified

  /**  */
  void clear_models_( bool called_from_destructor = false );

  /**  */
  void clear_prototypes_();

  /**  */
  index
  register_node_model_( Model* model, bool private_model = false );

  /**
   * Copy an existing node model and register it as a new model.
   * @param old_id ID of existing model.
   * @param new_name name of new model.
   * @return model ID of new Model object.
   * @see copy_model(), copy_synapse_model_()
   */
  index
  copy_node_model_( index old_id, Name new_name );

  /**
   * Copy an existing synapse model and register it as a new model.
   * @param old_id ID of existing model.
   * @param new_name name of new model.
   * @return model ID of new Model object.
   * @see copy_model(), copy_node_model_()
   */
  index
  copy_synapse_model_( index old_id, Name new_name );

  /**
   * Set the default parameters of a model.
   * @param model_id of model.
   * @param params default parameters to be set. 
   * @see set_model_defaults, set_synapse_defaults_
   */
  void
  set_node_defaults_(index model_id, const DictionaryDatum& params );

  /**
   * Set the default parameters of a model.
   * @param name of model.
   * @param params default parameters to be set. 
   * @see set_model_defaults, set_node_defaults_
   */
  void
  set_synapse_defaults_( index model_id, const DictionaryDatum& params );
};

} // namespace nest

#endif /* MODEL_MANAGER_H */
