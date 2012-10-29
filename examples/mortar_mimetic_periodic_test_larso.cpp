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
#include <dune/porsol/mortar/IncompFlowSolverHybridMortar.hpp>
#include <dune/porsol/mortar/PeriodicHelpersMortar.hpp>


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
  typedef Dune::IncompFlowSolverHybrid<GI, RI, BCs, Dune::MimeticIPEvaluator> 
    FlowSolverOrig;
  typedef Dune::IncompFlowSolverHybridMortar<GI, RI, BCs, Dune::MimeticIPEvaluator> 
    FlowSolverMortar;

  typedef double ctype; // To be used in shapefunctions

  // Check input --------------------------------

  // Parameters
  double ztol = 1e-8;
  const int dim = 3; 
  bool printSoln = true;
  bool vtk = false;
  
  CpGrid grid;
  ReservoirPropertyCapillary<dim> rockParams;

  // Check input. If no input given, create cartesian grid and uniform rock parameters
  if (varnum > 2) { //Wrong nr of input param
    cerr << "Wrong nr of input parameters. Only eclipse grid file name is allowed. If no grid file provided, a uniform cartesian grid is constructed." << endl;
    exit(1);
  }
  else if (varnum == 1) {
    array<int,3> dims = {{3, 3, 3}};
    array<double,3> cellsize = {{1, 1, 1}};
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

    if (! (eclParser.hasField("SPECGRID") 
	   && eclParser.hasField("COORD") 
	   && eclParser.hasField("ZCORN"))) {  
      cerr << "Error: Did not find SPECGRID, COORD and ZCORN in Eclipse file " 
	   << ECLIPSEFILENAME << endl;  
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

  CI::Vector gravity;
  gravity[0] = gravity[1] = gravity[2] = 0.0;
  
  grid.setUniqueBoundaryIds(true);
  GridInterfaceEuler<CpGrid> g(grid);

  // Print gridinfo
  /*
  cout << "Grid info:\n";
  for (CI cell = g.cellbegin(); cell != g.cellend(); ++cell) {
    cout << cell->index() << "\t" << cell->centroid() << endl;
    int can_pos = 0;
    for (FI face = cell.facebegin(); face != cell.faceend(); ++face, ++can_pos) {
      cout << "  " << can_pos << "\t" << g.faceIndex(cell->index(), can_pos) << endl;
      cout << "    " << face->centroid() << endl;
    }
  }
  */
  
  // Set up Boundary Conditions
  array<FlowBC, 6> cond = {{ FlowBC(FlowBC::Periodic, 1.0*Opm::unit::barsa),
			     FlowBC(FlowBC::Periodic,-1.0*Opm::unit::barsa),
			     //FlowBC(FlowBC::Dirichlet, 1.0*Opm::unit::barsa),
			     //FlowBC(FlowBC::Dirichlet,-1.0*Opm::unit::barsa),
			     FlowBC(FlowBC::Periodic, 0.0),
			     FlowBC(FlowBC::Periodic, 0.0),
			     //FlowBC(FlowBC::Dirichlet, 0.0),
			     //FlowBC(FlowBC::Dirichlet, 0.0), 
			     //FlowBC(FlowBC::Periodic,-1.0*Opm::unit::barsa),
			     //FlowBC(FlowBC::Periodic, 1.0*Opm::unit::barsa),
			     FlowBC(FlowBC::Periodic, 0.0),
			     FlowBC(FlowBC::Periodic, 0.0) }};
  //FlowBC(FlowBC::Periodic, 1.0*Opm::unit::barsa),
  //FlowBC(FlowBC::Periodic,-1.0*Opm::unit::barsa) }};
  
  BCs fbc_orig;
  BCs fbc_mortar;

  createPeriodic(fbc_orig, g, cond);
  createPeriodicMortar(fbc_mortar, g, cond);

  // Init flow solvers
  FlowSolverOrig   solver_orig;
  FlowSolverMortar solver_mortar;
  solver_orig.init(g, rockParams, gravity, fbc_orig);
  solver_mortar.init(g, rockParams, gravity, fbc_mortar);

  // Print Mortar Matrix
  if (printSoln && numCells < 30) {
    solver_mortar.mortar_.printMortarMatrix(0);
    solver_mortar.mortar_.printMortarMatrix(1);
  }
  
  /*
  solver_mortar.mortar_.printFace(1);
  solver_mortar.mortar_.printFace(2);
  solver_mortar.mortar_.printFace(3);
  solver_mortar.mortar_.printFace(4);
  */

  vector<double> src(numCells, 0.0);
  vector<double> sat(numCells, 0.0);

  // Call solvers
  solver_orig.solve(rockParams, sat, fbc_orig, src, 1e-8, 1, 0);
  solver_orig.printStats(std::cout);
  solver_orig.printSystem("orig");

  solver_mortar.solve(rockParams, sat, fbc_mortar, src, 1e-8, 1, 0);
  //solver_mortar.printStats(std::cout);
  solver_mortar.printSystem("mortar");
  double maxIntDiff = solver_mortar.checkFluxPeriodicity();

  FlowSolverOrig::SolutionType soln_orig = solver_orig.getSolution();
  FlowSolverMortar::SolutionType soln_mortar = solver_mortar.getSolution();
  
  // Print solutions
  if (printSoln && numCells < 130) {
    cout << "\nCell pressure orig and mortar:\n" << scientific << setprecision(3);
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
      cout << c->index() 
	   << '\t' << soln_orig.pressure(c)
	   << "\t" << soln_mortar.pressure(c)
      	   << "\t" << soln_orig.pressure(c)-soln_mortar.pressure(c) << '\n';
    }

    cout << "\nOutlux orig and mortar:\n";
    for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
      for (FI f = c.facebegin(); f != c.faceend(); ++f) {
	cout << f->index()
	     << "\t" << soln_orig.outflux(f)
	     << "\t" << soln_mortar.outflux(f)
	     << "\t" << soln_orig.outflux(f)-soln_mortar.outflux(f) << "\n";
      }
    }

  }

  vector<GI::Scalar> cellPressureOrig;
  vector<GI::Scalar> cellPressureMortar;
  for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
    cellPressureOrig.push_back(soln_orig.pressure(c)/(1.0*Opm::unit::barsa));
    cellPressureMortar.push_back(soln_mortar.pressure(c)/(1.0*Opm::unit::barsa));
  }

  // VTK writer
  if (vtk) {
    string vtufile_orig   = "orig_output";
    string vtufile_mortar = "mortar_output";
    VTKWriter<CpGrid::LeafGridView> writer_orig(grid.leafView());
    VTKWriter<CpGrid::LeafGridView> writer_mortar(grid.leafView());
    writer_orig.addCellData(cellPressureOrig, "cellPressure");
    writer_mortar.addCellData(cellPressureMortar, "cellPressure");
    writer_orig.write(vtufile_orig);
    writer_mortar.write(vtufile_mortar);
  }


  return 0;
}
