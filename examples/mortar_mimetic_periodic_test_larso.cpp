// This code is for testing Mortar Mimetic FDM

#include "config.h"

#include <opm/core/utility/have_boost_redef.hpp>

#include <iostream>

#include <boost/array.hpp>

#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <dune/grid/CpGrid.hpp>

#include <dune/porsol/common/PeriodicHelpers.hpp>
#include <dune/porsol/common/BoundaryConditions.hpp>
#include <dune/porsol/common/GridInterfaceEuler.hpp>
#include <dune/porsol/common/ReservoirPropertyCapillary.hpp>

#include <dune/porsol/mimetic/MimeticIPEvaluator.hpp>
#include <dune/porsol/mimetic/IncompFlowSolverHybrid.hpp>

using namespace Dune;

int main(int argc, char** argv)
{

  return 0;
}
