// This code is for testing Mortar Mimetic FDM on uniform grid

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
  int dir_pdrop = 0;
  
  CpGrid grid;
  ReservoirPropertyCapillary<dim> rockParams;

  // Check input. If no input given, create cartesian grid and uniform rock parameters
  if (varnum > 2) { //Wrong nr of input param
    cerr << "Wrong nr of input parameters. Only grid dim is taken as input" << endl;
    exit(1);
  }
  else {
    array<int,3> dims = {{3, 3, 3}};
    if (varnum == 2) {
      int dim = atoi(vararg[1]);
      dims = {{dim, dim, dim}};
      if (dim < 3) {
	cout << "Dimension must be 3 or higher. Setting dim = 3." << endl;
	dim = 3;
      }
      else {
	cout << "Dimension set to " << dim << endl;
      }
      dims = {{dim, dim, dim}};
    }
    array<double,3> cellsize = {{1, 1, 1}};
    grid.createCartesian(dims, cellsize);
    double uniformPORO = 0.2;
    double uniformPERM = 100.0 *Opm::prefix::milli *Opm::unit::darcy;
    rockParams.init(grid.size(0), uniformPORO, uniformPERM);
  }
      
  int numCells = grid.size(0);

  // Gravity and source/sat
  CI::Vector gravity;
  gravity[0] = gravity[1] = gravity[2] = 0.0;
  vector<double> src(numCells, 0.0);
  vector<double> sat(numCells, 0.0);
  
  grid.setUniqueBoundaryIds(true);
  GridInterfaceEuler<CpGrid> g(grid);

  vector<string> pdropString;
  pdropString.push_back("-xdir");
  pdropString.push_back("-ydir");
  pdropString.push_back("-zdir");

  // Set up BCs for pdrop in x direction
  array<FlowBC,6> cond = {{ FlowBC(FlowBC::Periodic, 1.0*Opm::unit::barsa),
			    FlowBC(FlowBC::Periodic,-1.0*Opm::unit::barsa),
			    FlowBC(FlowBC::Periodic, 0.0),
			    FlowBC(FlowBC::Periodic, 0.0),
			    FlowBC(FlowBC::Periodic, 0.0),
			    FlowBC(FlowBC::Periodic, 0.0) }};

  BCs fbc_orig;
  BCs fbc_mortar;
  createPeriodic(fbc_orig, g, cond);
  createPeriodicMortar(fbc_mortar, g, cond);

  // Init flow solvers
  FlowSolverOrig   solver_orig;
  FlowSolverMortar solver_mortar;
  solver_orig.init(g, rockParams, gravity, fbc_orig);
  solver_orig.printStats(std::cout);
  solver_mortar.init(g, rockParams, gravity, fbc_mortar);  

  // Write mortar matrices to matlab
  string dimString = "uniform" + to_string(dim) + "x" + to_string(dim) + "x" + to_string(dim) + "grid";
  writeMatrixToMatlab(solver_mortar.mortar_.getMortarMatrices()[0], "L1-" + dimString + ".dat");
  writeMatrixToMatlab(solver_mortar.mortar_.getMortarMatrices()[0], "L2-" + dimString + ".dat");

  // For each direction
  for (int dir_pdrop=0; dir_pdrop<3; ++dir_pdrop) {

    if (dir_pdrop==0)  
      cout << "##### Solving for pressure drop in x direction #####" << endl << endl;
    else if (dir_pdrop==1)
      cout << "##### Solving for pressure drop in y direction #####" << endl << endl;
    else if (dir_pdrop==2)
      cout << "##### Solving for pressure drop in z direction #####" << endl << endl;
    else {
      cout << "Error! Should not be here!" << endl;
      exit(1);
    }

    // Call solvers
    solver_orig.solve(rockParams, sat, fbc_orig, src, 1e-8, 1, 0);
    double maxMod_orig = solver_orig.postProcessFluxes();
    
    solver_mortar.solve(rockParams, sat, fbc_mortar, src, 1e-8, 1, 0);
    double maxMod_mortar = solver_mortar.postProcessFluxes();
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
    }

    solver_orig.printSystem("orig-" + dimString + pdropString[dir_pdrop]);
    solver_mortar.printSystem("mortar-" + dimString + pdropString[dir_pdrop]);

    cout << endl << endl;

    // Update BCs to next iteration
    array<FlowBC,6> new_cond = cond;
    new_cond[0] = cond[4];
    new_cond[1] = cond[5];
    new_cond[2] = cond[0];
    new_cond[3] = cond[1];
    new_cond[4] = cond[2];
    new_cond[5] = cond[3];
    cond = new_cond;
     
    fbc_orig.clear();
    fbc_mortar.clear();
    createPeriodic(fbc_orig, g, cond);
    createPeriodicMortar(fbc_mortar, g, cond);
  }

  /*
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
  */

  /*
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
  */


  return 0;
}
