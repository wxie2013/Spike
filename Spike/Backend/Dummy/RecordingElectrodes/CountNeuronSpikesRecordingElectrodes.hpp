#pragma once

#include "Spike/RecordingElectrodes/CountNeuronSpikesRecordingElectrodes.hpp"
#include "RecordingElectrodes.hpp"

namespace Backend {
  namespace Dummy {
    class CountNeuronSpikesRecordingElectrodes :
      public virtual ::Backend::Dummy::RecordingElectrodesCommon,
      public ::Backend::CountNeuronSpikesRecordingElectrodes {
    public:
      MAKE_BACKEND_CONSTRUCTOR(CountNeuronSpikesRecordingElectrodes);

      virtual void reset_state() {
        // TODO
      }

      virtual void add_spikes_to_per_neuron_spike_count
      (::CountNeuronSpikesRecordingElectrodes* front, float current_time_in_seconds) {
        printf("TODO add_spikes_to_per_neuron_spike_count\n");
      }
    };
  }
}
