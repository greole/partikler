#include <vector>
#include "SPHModels.hpp"

SPHModelFactory::map_type *SPHModelFactory::map_;

// ModelRegister<SPHModelGraph> SPHModelGraph::reg("Core::SPHModelGraph");

REGISTER_DEF_TYPE(CORE, SPHModelGraph);
