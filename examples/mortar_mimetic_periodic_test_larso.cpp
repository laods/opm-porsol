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
//#include <dune/porsol/mimetic/IncompFlowSolverHybrid.hpp>
#include <dune/porsol/mortar/IncompFlowSolverHybridMortar.hpp>

#include <dune/porsol/mortar/mortar.hpp>
#include <dune/porsol/mortar/shapefunctions.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

using namespace Dune;
using namespace std;

int main(int varnum, char** vararg)
{

  typedef Dune::GridInterfaceEuler<CpGrid>                       GI;
  typedef GI  ::CellIterator                                     CI;
  typedef CI  ::FaceIterator                                     FI;
  typedef Dune::BasicBoundaryConditions<true, false>             BCs;
  typedef Dune::ReservoirPropertyCapillary<3>                    RI;
  typedef Dune::IncompFlowSolverHybridMortar<GI, RI, BCs, Dune::MimeticIPEvaluator> FlowSolver;

  typedef double ctype; // To be used in shapefunctions

  // Check input --------------------------------

  // Parameters
  double ztol = 1e-8;
  const int dim = 3; 
  bool printSoln = false;
  bool vtk = false;
  
  CpGrid grid;
  ReservoirPropertyCapillary<dim> rockParams;

  if (varnum > 2) { //Wrong nr of input param
    cerr << "Wrong nr of input parameters. Only eclipse grid file name is allowed. If no grid file provided, a uniform cartesian grid is constructed." << endl;
    exit(1);
  }
  else if (varnum == 1) {
    array<int,3> dims = {10, 10, 5};
    array<double,3> cellsize = {0.05, 0.05, 0.02};
    grid.createCartesian(dims, cellsize);
    double uniformPORO = 0.2;
    double uniformPERM = 100.0 *Opm::prefix::milli *Opm::unit::darcy;
    rockParams.init(grid.size(0), uniformPORO, uniformPERM);
  }
  
  else {

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

    grid.readEclipseFormat(ECLIPSEFILENAME, ztol, false);
    rockParams.init(eclParser, grid.globalCell());
  }  
    
  int numCells = grid.size(0);

  if (numCells < 10) {
    cout << "PORO and PERMX:" << endl;
    for (int c = 0; c < numCells; ++c) {
      cout << "Cell " << c << ":\n" << rockParams.porosity(c) << "\n" << rockParams.permeability(c) << endl;
    }
  }


  CI::Vector gravity;
  gravity[0] = gravity[1] = gravity[2] = 0.0;
  
  grid.setUniqueBoundaryIds(true);
  GridInterfaceEuler<CpGrid> g(grid);
  
  array<FlowBC, 6> cond = {{ FlowBC(FlowBC::Periodic, 1.0*Opm::unit::barsa),
			     FlowBC(FlowBC::Periodic,-1.0*Opm::unit::barsa),
			     //FlowBC(FlowBC::Dirichlet, 1.0*Opm::unit::barsa),
			     //FlowBC(FlowBC::Dirichlet,-1.0*Opm::unit::barsa),
			     FlowBC(FlowBC::Periodic, 0.0),
			     FlowBC(FlowBC::Periodic, 0.0),
			     //FlowBC(FlowBC::Dirichlet, 0.0),
			     //FlowBC(FlowBC::Dirichlet, 0.0), 
			     FlowBC(FlowBC::Periodic, 0.0),
			     FlowBC(FlowBC::Periodic, 0.0) }};

  BCs fbc;
  createPeriodic(fbc, g, cond);

  FlowSolver solver;
  solver.init(g, rockParams, gravity, fbc);

  vector<double> src(numCells, 0.0);
  vector<double> sat(numCells, 0.0);

  solver.solve(rockParams, sat, fbc, src);

  double mod = solver.postProcessFluxes();
  cout << "Max mod: " << mod << endl;

  solver.printStats(std::cout);

  // Play with shape functions
  P1ShapeFunctionSet<ctype,ctype,2> lbasis = P1ShapeFunctionSet<ctype,ctype,2>::instance();
  P0ShapeFunctionSet<ctype,ctype,2> ubasis = P0ShapeFunctionSet<ctype,ctype,2>::instance();

  Dune::FieldVector<ctype,2> pkt(0.5);
  cout << "Size ubasis: " << ubasis.size() << endl;
  cout << "ubasis = " << ubasis[0].evaluateFunction(pkt) << endl;
  cout << "ugrad  = " << ubasis[0].evaluateGradient(pkt) << endl;
  cout << "Size lbasis: " << lbasis.size() << endl;
  cout << "lbasis = " << lbasis[0].evaluateFunction(pkt) << endl;
  cout << "lgrad  = " << lbasis[0].evaluateGradient(pkt) << endl;

  // solver.init(...) // Må gjøre noe med template BCs til FlowSolver!
  // Lage ny BCs class ?
  // Lage ny IncompFlowSolverHybrid class ?

  FlowSolver::SolutionType soln = solver.getSolution();

  // Print solution
  if (printSoln) {
    cout << "Cell Pressure:\n" << scientific << setprecision(15);
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
      cout << '\t' << soln.pressure(c) << '\n';
    }
    
    cout << "Cell (Out) Fluxes:\n";
    cout << "flux = [\n";
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
      for (FI f = c->facebegin(); f != c->faceend(); ++f) {
	cout << soln.outflux(f) << ' ';
      }
      cout << "\b\n";
    }
    cout << "]\n";
  } 
   
  vector<GI::Scalar> cellPressure;
  for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
    cellPressure.push_back(soln.pressure(c)/(1.0*Opm::unit::barsa));
  }
  
  // VTK writer
  if (vtk) {
    string vtufile = "mortar_output";
    VTKWriter<CpGrid::LeafGridView> writer(grid.leafView());
    writer.addCellData(cellPressure, "cellPressure");
    writer.write(vtufile);
  }

  MortarHelper<CpGrid> mortar(grid);
  cout << "min = " << mortar.min()[0] << " " << mortar.min()[1] << " " << mortar.min()[2] << endl;
  cout << "max = " << mortar.max()[0] << " " << mortar.max()[1] << " " << mortar.max()[2] << endl;
  cout << "n   = " << mortar.n1() << " " << mortar.n2() << endl;

  //mortar.printFace(1);
  //mortar.printFace(2);
  //mortar.printFace(6);
  mortar.periodicBCsMortar(1060);

  return 0;
}
