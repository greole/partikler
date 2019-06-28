#include <vector>
#include "SPHModels.hpp"

SPHModelFactory::map_type *SPHModelFactory::map_;

 ModelRegister<SPHModelGraph> SPHModelGraph::reg("SPHModelGraph");
