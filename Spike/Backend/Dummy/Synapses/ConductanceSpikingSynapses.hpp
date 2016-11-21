#pragma once

#include "SpikingSynapses.hpp"
#include "Spike/Synapses/ConductanceSpikingSynapses.hpp"

namespace Backend {
  namespace Dummy {
    class ConductanceSpikingSynapses : public virtual ::Backend::Dummy::SpikingSynapsesCommon,
                                       public ::Backend::ConductanceSpikingSynapses
    {
    public:
      virtual void prepare() {
        printf("TODO Backend::Dummy::ConductanceSpikingSynapses::prepare\n");
      }

      virtual void calculate_postsynaptic_current_injection(::SpikingNeurons * neurons, float current_time_in_seconds, float timestep) {
        // printf("TODO Dummy::ConductanceSpikingSynapses::calculate_postsynaptic_current_injection\n");
      }

      virtual void update_synaptic_conductances(float timestep, float current_time_in_seconds) {
        // printf("TODO Dummy::ConductanceSpikingSynapses::update_synaptic_conductances\n");
      }

      virtual void reset_state() {}
    };
  } // namespace Dummy
} // namespace Backend

