//===========================================================================
//
// File: mortar_mimetic_uniform_test.cpp
//
// Author(s): Lars Vingli Odsæter <lars.odsater@gmail.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  This file is part of The Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
  Tests the mortar approach implemented in dune/porsol/mortar/mortar.hpp
  on a uniform grid. Input variables:
    1) Grid size
*/

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
  bool vtk = true;
  int dir_pdrop = 0;
  int gridsize = 3;  

  CpGrid grid;
  ReservoirPropertyCapillary<dim> rockParams;

  // Check input. If no input given, create cartesian grid and uniform rock parameters
  if (varnum > 2) { //Wrong nr of input param
    cerr << "Wrong nr of input parameters. Only grid size is taken as input" << endl;
    exit(1);
  }
  else {
    array<int,3> dims = {{3, 3, 3}};
    if (varnum == 2) {
      gridsize = atoi(vararg[1]);
      if (gridsize < 3) {
	cout << "Grid size must be 3 or higher. Setting gridsize = 3." << endl;
	gridsize = 3;
      }
      else {
	cout << "Gridsize set to " << gridsize << endl;
      }
      dims[0] = gridsize;
      dims[1] = gridsize;
      dims[2] = gridsize;
    }
    array<double,3> cellsize = {{1.0/gridsize, 1.0/gridsize, 1.0/gridsize}};
    grid.createCartesian(dims, cellsize);
    double uniformPORO = 0.2;
    //double uniformPERM = 100.0 *Opm::prefix::milli *Opm::unit::darcy;
    double uniformPERM = 1.0; 
    rockParams.init(grid.size(0), uniformPORO, uniformPERM);
    
    for (int zlayer=0; zlayer<dims[2]; ++zlayer) {
      cout << double(zlayer+1)/gridsize << endl;
      if (double(zlayer+1)/gridsize<0.26 || double(zlayer+1)/gridsize>0.76)
	continue;
      int cell_idx = zlayer*dims[0]*dims[1];
      for (int i=0; i<dims[0]*dims[1]; ++i) {
        rockParams.permeabilityModifiable(cell_idx) *= 10;
	++cell_idx;
      }
    }
   
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
  if (numCells<126) 
    solver_orig.printStats(std::cout);
  cout << endl;
  solver_mortar.init(g, rockParams, gravity, fbc_mortar);  

  // Write mortar matrices to matlab
  string dimString = "uniform" + to_string(static_cast<long long>(gridsize)) + "x" + 
    to_string(static_cast<long long>(gridsize)) + "x" + 
    to_string(static_cast<long long>(gridsize)) + "grid";
  writeMatrixToMatlab(solver_mortar.mortar_.getMortarMatrices()[0], "L1-" + dimString + ".dat");
  writeMatrixToMatlab(solver_mortar.mortar_.getMortarMatrices()[0], "L2-" + dimString + ".dat");

  vector<vector<GI::Scalar> > cellPressureOrig;
  vector<vector<GI::Scalar> > cellPressureMortar;
  vector<vector<GI::Scalar> > cellVelocity;

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
    solver_orig.solve(rockParams, sat, fbc_orig, src, 1e-14, 1, 0);
    solver_mortar.solve(rockParams, sat, fbc_mortar, src, 1e-14, 1, 0);
    cout << "Flux integrals mortar:\n";
    double maxIntDiff_mortar = solver_mortar.checkFluxPeriodicity();
    double maxMod_mortar = solver_mortar.postProcessFluxes();
    cout << "Flux integrals orig:\n";
    double maxIntDiff_orig = solver_orig.checkFluxPeriodicity();
    double maxMod_orig = solver_orig.postProcessFluxes();

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
      
      vector<double> pi_orig = solver_orig.getContactPressureSoln();
      vector<double> pi_mortar = solver_mortar.getContactPressureSoln();
      cout << "\nContact pressures orig and mortar:\n";
      cout << pi_orig.size() << " " << pi_mortar.size() << endl;
      for (int i=0; i<pi_orig.size(); ++i) {
	cout << i
	     << "\t" << pi_orig[i]
	     << "\t" << pi_mortar[i]
	     << "\t" << pi_orig[i]-pi_mortar[i] << endl;
      }

      
    }

    solver_orig.printSystem("orig-" + dimString + pdropString[dir_pdrop]);
    solver_mortar.printSystem("mortar-" + dimString + pdropString[dir_pdrop]);

    if (vtk) {
      vector<GI::Scalar> temp;
      cellPressureOrig.push_back(temp);
      cellPressureMortar.push_back(temp);
      cellVelocity.push_back(temp);
      for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
	cellPressureOrig.back().push_back(soln_orig.pressure(c)/(1.0*Opm::unit::barsa) 
					  + 1.0 - 1.0/gridsize);
	cellPressureMortar.back().push_back(soln_mortar.pressure(c)/(1.0*Opm::unit::barsa)
					    + 1.0 - 1.0/gridsize);
	GI::Vector velocity(0.0);
	
	double totalWeight = 0.0;
	for (FI f = c.facebegin(); f != c.faceend(); ++f) {
	  GI::Vector diff = c.centroid() - f.centroid();
	  totalWeight += diff.two_norm();
	}

	for (FI f = c.facebegin(); f != c.faceend(); ++f) {
	  GI::Vector diff = c.centroid() - f.centroid();
	  //double weight = 1/(diff.two_norm());
	  // The original weights don't work as intended. 
	  double weight = 0.5; 
	  GI::Vector faceContrib = f.normal();
	  faceContrib *= weight*soln_mortar.outflux(f);
	  velocity += faceContrib;
	}
	cellVelocity.back().push_back(velocity[0]);
	cellVelocity.back().push_back(velocity[1]);
	cellVelocity.back().push_back(velocity[2]);

      }

    }

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

  vector<GI::Scalar> permField(numCells, 0.0);
  for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
    int c_idx = c.index();
    permField[c_idx] = rockParams.permeability(c_idx)(0,0);
  }

  if (vtk) {
    string vtufile = "solution_" + dimString;
    VTKWriter<CpGrid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.addCellData(permField, "permeability");
    vtkwriter.addCellData(cellPressureOrig[0], "cellPressureOrig-pdd0");
    vtkwriter.addCellData(cellPressureOrig[1], "cellPressureOrig-pdd1");
    vtkwriter.addCellData(cellPressureOrig[2], "cellPressureOrig-pdd2");
    vtkwriter.addCellData(cellPressureMortar[0], "cellPressureMortar-pdd0");
    vtkwriter.addCellData(cellPressureMortar[1], "cellPressureMortar-pdd1");
    vtkwriter.addCellData(cellPressureMortar[2], "cellPressureMortar-pdd2");
    vtkwriter.addCellData(cellVelocity[0], "cellVelocityMortar-pdd0", 3);
    vtkwriter.addCellData(cellVelocity[1], "cellVelocityMortar-pdd1", 3);
    vtkwriter.addCellData(cellVelocity[2], "cellVelocityMortar-pdd2", 3);
    vtkwriter.write(vtufile);
  }


  return 0;
}
