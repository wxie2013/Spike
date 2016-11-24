// Higgins STDP Class Header
// HigginsSTDP.h
//
//	Author: Nasir Ahmad
//	Date: 23/06/2016

#ifndef HIGGINS_STDP_H
#define HIGGINS_STDP_H

// Get Synapses & Neurons Class
#include "../Synapses/SpikingSynapses.hpp"
#include "../Neurons/SpikingNeurons.hpp"
#include "../STDP/STDP.hpp"

// stdlib allows random numbers
#include <stdlib.h>
// Input Output
#include <stdio.h>
// allows maths
#include <math.h>

namespace Backend {
  class HigginsSTDP : public virtual STDPCommon,
                      public STDP {
  public:
    virtual void prepare() {
      printf("TODO Backend::HigginsSTDP::prepare\n");
    }

    virtual void apply_ltp_to_synapse_weights(float* d_last_spike_time_of_each_neuron, float current_time_in_seconds) = 0;
    virtual void apply_ltd_to_synapse_weights(float* d_last_spike_time_of_each_neuron, float current_time_in_seconds) = 0;
  };
}

#include "Spike/Backend/Dummy/STDP/HigginsSTDP.hpp"

// STDP Parameters
struct higgins_stdp_parameters_struct : stdp_parameters_struct {
  // STDP Parameters
  float w_max = 60.0;
  float a_minus = -0.015;
  float a_plus = 0.005;
  float tau_minus = 0.025;
  float tau_plus = 0.015;
};


class HigginsSTDP : public STDP{
public:
  ~HigginsSTDP();
  ADD_BACKEND_GETTER(HigginsSTDP);

  struct higgins_stdp_parameters_struct* stdp_params = NULL;
  SpikingSynapses* syns = NULL;

  virtual void prepare_backend(Context* ctx = _global_ctx);

  // Set STDP Parameters
  virtual void Set_STDP_Parameters(SpikingSynapses* synapses, SpikingNeurons* neurons, SpikingNeurons* input_neurons, stdp_parameters_struct* stdp_parameters);
  // STDP
  virtual void Run_STDP(SpikingNeurons* neurons, float current_time_in_seconds, float timestep);

  // LTP & LTD for this model
  void apply_ltd_to_synapse_weights(float* d_last_spike_time_of_each_neuron, float current_time_in_seconds);
  void apply_ltp_to_synapse_weights(float* d_last_spike_time_of_each_neuron, float current_time_in_seconds);
};

#endif
