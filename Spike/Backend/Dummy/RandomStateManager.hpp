#pragma once

#include "Spike/Helpers/RandomStateManager.hpp"

namespace Backend {
  namespace Dummy {
    class RandomStateManager : public ::Backend::RandomStateManager {
    public:
      MAKE_BACKEND_CONSTRUCTOR(RandomStateManager);

      virtual void prepare() {
        printf("TODO Backend::Dummy::RandomStateManager::prepare\n");
      }
    };
  }
}
