#include <dune/porsol/mimetic/IncompFlowSolverHybrid.hpp>
#include <dune/porsol/mortar/mortar.hpp>


namespace Dune {

  template <class GridInterface,
	    class RockInterface,
	    class BCInterface,
	    template <class GridIF, class RockIF> class InnerProduct>
  class IncompFlowSolverHybridMortar : public IncompFlowSolverHybrid<GridInterface, 
								     RockInterface,
								     BCInterface,
								     InnerProduct>
  {

  public:

    MortarHelper<typename GridInterface::GridType> mortar_;

  };
  

} // Namespace Dune
