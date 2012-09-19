// This code is for testing Mortar Mimetic FDM

#include "config.h"

#include <opm/core/utility/have_boost_redef.hpp>

#include <iostream>

#include <boost/array.hpp>

#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/eclipse/EclipseGridParser.hpp>

#include <dune/grid/CpGrid.hpp>

#include <dune/porsol/common/PeriodicHelpers.hpp>
#include <dune/porsol/common/BoundaryConditions.hpp>
#include <dune/porsol/common/GridInterfaceEuler.hpp>
#include <dune/porsol/common/ReservoirPropertyCapillary.hpp>

#include <dune/porsol/mimetic/MimeticIPEvaluator.hpp>
#include <dune/porsol/mimetic/IncompFlowSolverHybrid.hpp>

using namespace Dune;
using namespace std;

int main(int varnum, char** vararg)
{

  if (varnum != 2) { //Wrong nr of input param
    cerr << "Wrong nr of input parameters. Only eclipse grid file name is allowed." << endl;
    exit(1);
  }

  const char* ECLIPSEFILENAME(vararg[1]);
  ifstream eclipsefile(ECLIPSEFILENAME, ios::in);
  if (eclipsefile.fail()) {
    cerr << "Error: " << ECLIPSEFILENAME << " not found or readable." << endl;
    exit(1);
  }
  eclipsefile.close();
  
  cout << "Parsing Eclipse file: " << ECLIPSEFILENAME << "..." << endl;
  Opm::EclipseGridParser eclParser(ECLIPSEFILENAME, true);

  if (! (eclParser.hasField("SPECGRID") && eclParser.hasField("COORD") && eclParser.hasField("ZCORN"))) {  
    cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " << ECLIPSEFILENAME << endl;  
    exit(1);  
  }
  if (! (eclParser.hasField("PERMX") && eclParser.hasField("PORO"))) {  
    cerr << "Error: Did not find PERMX and PORO in Eclipse file " << ECLIPSEFILENAME << endl;  
    exit(1);  
  }
  
  return 0;
}
