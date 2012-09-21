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

#include <dune/porsol/mortar/mortar.hpp>

using namespace Dune;
using namespace std;

int main(int varnum, char** vararg)
{

  typedef Dune::GridInterfaceEuler<CpGrid>                       GI;
  typedef GI  ::CellIterator                                     CI;
  typedef CI  ::FaceIterator                                     FI;
  typedef Dune::BasicBoundaryConditions<true, false>             BCs;
  typedef Dune::ReservoirPropertyCapillary<3>                    RI;
  typedef Dune::IncompFlowSolverHybrid<GI, RI, BCs, Dune::MimeticIPEvaluator> FlowSolver;

  // Check input ---------------------------------------------------------------------------------------

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
  
  // Parameters
  double ztol = 1e-8;
  const int dim = 3;  
  
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

  CpGrid grid;
  grid.readEclipseFormat(ECLIPSEFILENAME, ztol, false);
  int numCells = grid.size(0);

  ReservoirPropertyCapillary<dim> rockParams;
  rockParams.init(eclParser, grid.globalCell());
  
  if (numCells < 10) {
    cout << "PORO and PERMX:" << endl;
    for (int c = 0; c < numCells; ++c) {
      cout << "Cell " << c << ":\n" << rockParams.porosity(c) << "\n" << rockParams.permeability(c) << endl;
    }
  }

  MortarHelper<CpGrid> mortar(grid);
  cout << "min = " << mortar.min()[0] << " " << mortar.min()[1] << " " << mortar.min()[2] << endl;
  cout << "max = " << mortar.max()[0] << " " << mortar.max()[1] << " " << mortar.max()[2] << endl;
  cout << "n   = " << mortar.n1() << " " << mortar.n2() << endl;

  mortar.printFace(1);
  mortar.printFace(2);
  mortar.printFace(6);

  FlowSolver solver;
  // solver.init(...)
  // Må gjøre noe med template BCs til FlowSolver!
  // Lage ny BCs class ?
  // Lage ny IncompFlowSolverHybrid class ?

  return 0;
}
