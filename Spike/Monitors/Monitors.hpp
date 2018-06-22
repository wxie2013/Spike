#ifndef Monitors_H
#define Monitors_H

#include <string>
using namespace std;

class Monitors; // forward definition
#include "Spike/Base.hpp"

#include "Spike/Backend/Macros.hpp"
#include "Spike/Backend/Context.hpp"
#include "Spike/Backend/Backend.hpp"
#include "Spike/Backend/Device.hpp"

#include "Spike/Models/SpikingModel.hpp"
#include "Spike/Neurons/SpikingNeurons.hpp"
#include "Spike/Synapses/SpikingSynapses.hpp"


namespace Backend {
  class Monitors : public virtual SpikeBackendBase {
  public:
    SPIKE_ADD_BACKEND_FACTORY(Monitors);
  };
}

class Monitors : public virtual SpikeBase {
public:
  Monitors();
  ~Monitors() override = default;

  SPIKE_ADD_BACKEND_GETSET(Monitors, SpikeBase);
  void init_backend(Context* ctx = _global_ctx) override;
  void reset_state() override;

  SpikingModel* model = nullptr;

private:
  std::shared_ptr<::Backend::Monitors> _backend;
};

#endif
