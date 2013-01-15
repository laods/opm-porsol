//===========================================================================
//
// File: mortar.cpp
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


#ifndef MORTAR_HPP
#define MORTAR_HPP

#include <vector>

#include <dune/grid/CpGrid.hpp>
#include <dune/porsol/mortar/boundarygrid.cc>
#include <dune/porsol/mortar/matrixops.hpp>
#include <dune/porsol/mortar/shapefunctions.hh>

#include <dune/geometry/quadraturerules.hh>

//! \brief A sparse matrix holding our operator
typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > Matrix;

//! \brief A vector holding our RHS
typedef Dune::BlockVector<Dune::FieldVector<double,1> > Vector;

//! \brief For storing matrix adjacency/sparsity patterns
typedef std::vector< std::set<int> > AdjacencyPattern;

// Helper class for mortar methods
template<class GridInterface>
class MortarHelper 
{
public:

  typedef typename GridInterface::GridType GridType;
  static const int dim = GridType::dimension;

  typedef typename GridType::LeafGridView LeafGridView;
  typedef typename LeafGridView::ctype ctype;
  typedef typename LeafGridView::IndexSet LeafIndexSet;
  typedef typename LeafGridView::template Codim<dim>::Iterator LeafVertexIterator;
  typedef typename LeafGridView::template Codim<dim-1>::Iterator LeafFaceIterator;
  typedef typename LeafGridView::template Codim<0>::Iterator LeafIterator;
  typedef typename LeafGridView::template Codim<1>::Geometry::GlobalCoordinate GlobalCoordinate;
  typedef typename LeafGridView::template Codim<1>::Geometry::LocalCoordinate LocalCoordinate;

  // Default constructor
  MortarHelper() {};

  // Constructor 1
  MortarHelper(const GridInterface& grid)
    : pgrid_(&grid) 
  {
    findMinMax();
    find_n();
    tol_ = 1e-8;
    nEqns_ = 0;
    dco_ = true;
  };

  // Constructor 2
  MortarHelper(const GridInterface& grid, int nEqns, Opm::SparseTable<int> cellFaces)
    : pgrid_(&grid), nEqns_(nEqns), cellFaces_(cellFaces)
  {
    findMinMax();
    find_n();
    tol_ = 1e-8;
    dco_ = true;
  };

  // Constructor 3
  MortarHelper(const GridInterface& grid, std::vector<ctype> min,
	       std::vector<ctype> max, int n1, int n2_, int nEqns, 
	       Opm::SparseTable<int> cellFaces,
	       double tol = 1e-8, bool dco = true)
    : pgrid_(&grid), min_(min), max_(max), n1_(n1), n2_(n2), 
      nEqns_(nEqns), cellFaces_(cellFaces), tol_(tol), dco_(dco) {};

  // Initializer 1. Used in addition to default contructor.
  void init(const GridInterface& grid) 
  {
    pgrid_ = &grid;
    findMinMax();
    find_n();
    tol_ = 1e-8;
    nEqns_ = 0;
    dco_ = true;
  }

  // Initializer 2. Used in addition to default contructor.
  void init(const GridInterface& grid, int nEqns) 
  {
    pgrid_ = &grid;
    findMinMax();
    find_n();
    tol_ = 1e-8;
    nEqns_ = nEqns;
    dco_ = true;
  }

  // Initializer 3. Used in addition to default contructor.
  void init(const GridInterface& grid, int nEqns, Opm::SparseTable<int> cellFaces) 
  {
    pgrid_ = &grid;
    findMinMax();
    find_n();
    tol_ = 1e-8;
    nEqns_ = nEqns;
    cellFaces_ = cellFaces;
    dco_ = true;
  }
  
  // Initializer 4. Used in addition to default contructor.
  void init(const GridInterface& grid, std::vector<ctype> min,
	    std::vector<ctype> max, int n1, int n2_, int nEqns, 
	    Opm::SparseTable<int> cellFaces,
	    double tol = 1e-8, bool dco = true)
  {
    pgrid_ = &grid;
    min_ = min;  max_ = max;
    n1_  = n1;   n2 = n2_;
    tol_ = tol;
    nEqns_ = nEqns;
    cellFaces_ = cellFaces;
    dco_ = dco;
  }

  // Clear member variables
  void clear()
  {
    n1_ = n2_ = 0;
    min_ = max_ = std::vector<double>(3,0.0);
    pgrid_ = 0;
    nEqns_ = 0;
    tol_ = 1e-8;
    cellFaces_.clear();
    dco_ = true;
    master.clear();
    slave.clear();
    L.clear();
    rhs.clear();
  }

  // Calculate minimum coordinate of grid
  std::vector<double> min() {
    return min_;
  }

  // Calculate maximum coordinate of grid
  std::vector<double> max() {
    return max_;
  }

  // Set and get functions
  int n1() {
    return n1_;
  } 
  int n2() {
    return n2_;
  }
  void setMinMax(std::vector<double> min, std::vector<double> max) {
    min_ = min;
    max_ = max;
  }
  void set_n(int n1, int n2) {
    n1_ = n1;
    n2_ = n2;
  }
  const std::vector<Matrix> getMortarMatrices() {
    return L;
  }
  const Vector getRhs(int dir) {
    return rhs[dir];
  }

  void findMinMax();
  void find_n();
  
  // Test function for debugging
  // Print vertices on face quads with global index and coord
  void printFace(int face);

  // Print mortar matrix in direction dir
  void printMortarMatrix(int dir);

  // Calculate mortar matrices
  void periodicBCsMortar();

private:
  std::vector<double> min_; // minimum coordinate of grid
  std::vector<double> max_; // maximum coordinate of grid
  int n1_; // Number of discretization cells in x-dir
  int n2_; // Number of discretization cells in y-dir
  const GridInterface* pgrid_;
  int nEqns_; // Number of equations in system
  bool dco_; // Dune convention ordering of vertices in quad
  Opm::SparseTable<int> cellFaces_;
  double tol_; // Tolerance to decide if point is on line for instance

  enum SIDE {
    LEFT,
    RIGHT
  };

  enum Direction { NONE = 0, X = 1, Y = 2, Z = 4,
                     XY = 1+2, XZ = 1+4, YZ = 2+4, 
                    XYZ = 1+2+4 };

  std::vector<BoundaryGrid> master; // Master grids
  std::vector< std::vector<BoundaryGrid::Vertex> > slave; // Slave grids
  std::vector<Matrix> L; // Mortar matrices
  std::vector<Vector> rhs; // Right hand sides

  // Extract vertices on boundary face 
  std::vector<BoundaryGrid::Vertex> extractFace(Direction dir, ctype coord);

  // Extract boundary grid on master face in direction dir
  BoundaryGrid extractMasterFace(Direction dir, ctype coord, SIDE side=LEFT, bool dc=false);
  
  // Check if coord is on plane
  bool isOnPlane(Direction plane, GlobalCoordinate coord, ctype value);

  // Check if coord is on line
  bool isOnLine(Direction dir, GlobalCoordinate coord, ctype x, ctype y);

  // Check is coord is on point
  bool isOnPoint(GlobalCoordinate coord, GlobalCoordinate point);

  // Calculate mortar matrix in direction dir
  Matrix findLMatrixMortar(const BoundaryGrid& b1, const BoundaryGrid& interface, int dir);

  // Returns equation number for a given quad
  int getEquationForDof(const BoundaryGrid::Quad& quad);
  int getEquationForDof(const BoundaryGrid::Quad& quad, int dir);

  // Calculates centroid of quad defined by four vertices
  GlobalCoordinate centroid(std::vector<GlobalCoordinate> vertices, int dir);

}; // MortarHelper


template<class GridInterface>
void MortarHelper<GridInterface>::findMinMax() {
  // Based on ElasticityUpscale::findBoundaries()
  double big = 1e5;
  if (min_.empty()) min_.resize(3,big);
  else min_[0] = min_[1] = min_[2] = big;
  if (max_.empty()) max_.resize(3,-big);
  else max_[0] = max_[1] = max_[2] = -big;

  const LeafVertexIterator itend = pgrid_->grid().leafView().template end<dim>();
  // iterate over vertices and find slaves
  LeafVertexIterator start = pgrid_->grid().leafView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    for (int i=0;i<3;++i) {
      min_[i] = std::min(min_[i],it->geometry().corner(0)[i]);
      max_[i] = std::max(max_[i],it->geometry().corner(0)[i]);
    }
  }
}

template<class GridInterface>
void MortarHelper<GridInterface>::find_n() {
  n1_ = pgrid_->grid().logicalCartesianSize()[0];
  n2_ = pgrid_->grid().logicalCartesianSize()[1];
}

template<class GridInterface>
std::vector<BoundaryGrid::Vertex> MortarHelper<GridInterface>::extractFace(Direction dir, ctype coord) {
  // Based on ElasticityUpscale::extractFaces()
  GridType gv =  pgrid_->grid();
  std::vector<BoundaryGrid::Vertex> result;
  const LeafIndexSet& set = gv.leafView().indexSet();
  const LeafVertexIterator itend = gv.leafView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(gv);
  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    if (isOnPlane(dir,it->geometry().corner(0),coord)) {
      BoundaryGrid::Vertex v;
      v.i = mapper.map(*it);
      BoundaryGrid::extract(v.c,it->geometry().corner(0),log2(dir));
      result.push_back(v);
    }
  }

  return result;
}

template<class GridInterface>
BoundaryGrid MortarHelper<GridInterface>::extractMasterFace(Direction dir, 
						       ctype coord,
						       SIDE side,
						       bool dc) 
{
  // Based on ElasticityUpscale::extractMasterFace()
  GridType gv =  pgrid_->grid();

  static const int V1[3][4] = {{0,2,4,6},
                               {0,1,4,5},
                               {0,1,2,3}};
  static const int V2[3][4] = {{1,3,5,7},
                               {2,3,6,7},
                               {4,5,6,7}};
  const LeafIndexSet& set = gv.leafView().indexSet();
  const LeafVertexIterator itend = gv.leafView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(gv);
  LeafVertexIterator start = gv.leafView().template begin<dim>();
  LeafIterator cellend = gv.leafView().template end<0>();
  int c = 0;
  int i = log2(dir);
  BoundaryGrid result;
  for (LeafIterator cell  = gv.leafView().template begin<0>(); 
                    cell != cellend; ++cell, ++c) {
    std::vector<BoundaryGrid::Vertex> verts;
    int idx; 
    if (side == LEFT)
     idx = set.subIndex(*cell,V1[i][0],dim);
    else if (side == RIGHT)
     idx = set.subIndex(*cell,V2[i][0],dim);
    LeafVertexIterator it=start;
    for (it ; it != itend; ++it) {
      if (mapper.map(*it) == idx)
        break;
    }
    if (isOnPlane(dir,it->geometry().corner(0),coord)) {
      for (int j=0;j<4;++j) {
        if (side == LEFT)
          idx = set.subIndex(*cell,V1[i][j],dim);
        if (side == RIGHT)
          idx = set.subIndex(*cell,V2[i][j],dim);
        LeafVertexIterator it=start;
        for (it ; it != itend; ++it) {
          if (mapper.map(*it) == idx)
            break;
        }
        if (!isOnPlane(dir,it->geometry().corner(0),coord))
          continue;
        BoundaryGrid::Vertex v;
        BoundaryGrid::extract(v,it->geometry().corner(0),i);
        v.i = idx;
        verts.push_back(v);
      }
    }
    if (verts.size() == 4) {
      BoundaryGrid::Quad q;
      q.v[0] = minXminY(verts);
      q.v[1] = maxXminY(verts);
      if (dc) {
        q.v[2] = minXmaxY(verts);
        q.v[3] = maxXmaxY(verts);
      } else {
        q.v[2] = maxXmaxY(verts);
        q.v[3] = minXmaxY(verts);
      }

      // Find global index of quad(face) if mapping is present
      if (!cellFaces_.empty()) {
	if (side == LEFT) 
	  q.globalFaceIndex = cellFaces_[cell->index()][2*i];
	else
	  q.globalFaceIndex = cellFaces_[cell->index()][2*i+1];
      }
      result.add(q);
    }
  }

  return result;
}

template<class GridInterface>
bool MortarHelper<GridInterface>::isOnPlane(Direction plane,
                                            GlobalCoordinate coord,
                                            ctype value)
{
  if (plane < X || plane > Z)
    return false;
  int p = log2(plane);
  ctype delta = fabs(value-coord[p]);
  return delta < tol_;
}

template<class GridInterface>
bool MortarHelper<GridInterface>::isOnLine(Direction dir,
                                           GlobalCoordinate coord,
                                           ctype x, ctype y)
{
  if (dir < X || dir > Z)
    return false;
  int ix = int(log2(dir)+1) % 3;
  int iy = int(log2(dir)+2) % 3;
  ctype delta = x-coord[ix];
  if (delta > tol_ || delta < -tol_)
    return false;
  delta = y-coord[iy];
  if (delta > tol_ || delta < -tol_)
    return false;

  return true;
}

template<class GridInterface>
bool MortarHelper<GridInterface>::isOnPoint(GlobalCoordinate coord,
                                            GlobalCoordinate point)
{
  GlobalCoordinate delta = point-coord;
  return delta.one_norm() < tol_;
}

template<class GridInterface>
void MortarHelper<GridInterface>::printFace(int face) {
  BoundaryGrid bg;
  switch (face) {
  case 1: // xmin
    bg = extractMasterFace(X, min_[0], LEFT);
    break;
  case 2: // xmax
    bg = extractMasterFace(X, max_[0], RIGHT);
    break;
  case 3: // ymin
    bg = extractMasterFace(Y, min_[1], LEFT);
    break;
  case 4: // ymax
    bg = extractMasterFace(Y, max_[1], RIGHT);
    break;
  case 5: // zmin
    bg = extractMasterFace(Z, min_[2], LEFT);
    break;
  case 6: // zmax
    bg = extractMasterFace(Z, max_[2], RIGHT);
    break;
  default:
    std::cerr << "Face nr must be >0 and <7" << std::endl;
    return;
  }
  std::cout << bg << std::endl;
  for (int i=0; i<bg.size(); ++i) {
    BoundaryGrid::Quad quad = bg[i];
    std::cout << "  " << i+1 << ", Quad centroid: " << quad.pos(0.5,0.5) 
	      << ", global face index: " << quad.globalFaceIndex << std::endl;
  }
}

template<class GridInterface>
void MortarHelper<GridInterface>::printMortarMatrix(int dir) 
{
  ASSERT(dir == 0 || dir == 1);
  ASSERT(!L.empty());
 
  if (dir == 0) std::cout << "\nMortar matrix X direction:" << std::endl;
  else          std::cout << "\nMortar matrix Y direction:" << std::endl;

  std::cout << std::setprecision(3) << std::fixed;

  Matrix MM(L[dir]);
  int n = MM.N();
  int m = MM.M();
  ASSERT(n == nEqns_);
  for (int i = 0; i < n; ++i) {
    std::cout << i << ": ";
    for (int j = 0; j < m; ++j) {
      if (MM.exists(i,j)) std::cout << MM[i][j];
      else std::cout << 0.0;
      std::cout << " ";
    }
    std::cout << std::endl;
  }
}


template<class GridInterface>
void MortarHelper<GridInterface>::periodicBCsMortar() {

  // Based on ElasticityUpscale::periodicBCsMortar()
  // But not MPC part (that is only step 3-6)

  // Step 3: extracts and establishes a quad grid 
  //         for the left/right/front/back sides
  master.push_back(extractMasterFace(X, min_[0], LEFT, dco_));
  master.push_back(extractMasterFace(X, max_[0], RIGHT, dco_));
  master.push_back(extractMasterFace(Y, min_[1], LEFT, dco_));
  master.push_back(extractMasterFace(Y, max_[1], RIGHT, dco_));

  // Step 4: Establishes grids for the dual dofs
  BoundaryGrid::FaceCoord fmin, fmax;
  fmin[0] = min_[1]; fmin[1] = min_[2];
  fmax[0] = max_[1]; fmax[1] = max_[2];
  BoundaryGrid lambdax = BoundaryGrid::uniform(fmin, fmax, n2_, 1, dco_);
  
  fmin[0] = min_[0]; fmin[1] = min_[2];
  fmax[0] = max_[0]; fmax[1] = max_[2];
  BoundaryGrid lambday = BoundaryGrid::uniform(fmin, fmax, n1_, 1, dco_);

  // Step 5: Calculates the coupling matrix L1 between left/right sides
  Matrix L1_left  = findLMatrixMortar(master[0], lambdax, 0);
  Matrix L1_right = findLMatrixMortar(master[1], lambdax, 0);
  L.push_back(MatrixOps::Axpy(L1_left, L1_right, -1));

  // Step 6: Calculates the coupling matrix L2 between front/back sides
  Matrix L2_left  = findLMatrixMortar(master[2], lambday, 1);
  Matrix L2_right = findLMatrixMortar(master[3], lambday, 1);
  L.push_back(MatrixOps::Axpy(L2_left, L2_right, -1));
  
}


template<class GridInterface>
Matrix MortarHelper<GridInterface>::findLMatrixMortar(const BoundaryGrid& b1,
						 const BoundaryGrid& interface,
						 int dir)
{

  std::vector< std::set<int> > adj;
  adj.resize(nEqns_);

  // process pillar by pillar
  size_t per_pillar = b1.size()/interface.size();
  for (size_t p=0;p<interface.size();++p) {
    for (size_t q=0;q<per_pillar;++q) {
      for (size_t i=0;i<1;++i) {
	const BoundaryGrid::Quad& qu(b1[p+per_pillar*q]);
	int dof;
	if (cellFaces_.empty()) 
	  dof = getEquationForDof(qu, dir);
	else
	  dof = getEquationForDof(qu);
	for (int j=0;j<4;++j) {
	  adj[dof].insert(interface[p].v[j].i);
	}
      }
    }
  }

  Matrix B;
  MatrixOps::fromAdjacency(B,adj,nEqns_,interface.totalNodes());

  // get a set of P0 shape functions for the face pressures
  P0ShapeFunctionSet<ctype,ctype,2> pbasis = P0ShapeFunctionSet<ctype,ctype,2>::instance();
  // get a set of PN shape functions for multipliers
  P1ShapeFunctionSet<ctype,ctype,2> lbasis = P1ShapeFunctionSet<ctype,ctype,2>::instance();
  // get a reference element
  Dune::GeometryType gt;
  gt.makeCube(2);
  const Dune::template GenericReferenceElement<ctype,2> &ref =
    Dune::GenericReferenceElements<ctype,2>::general(gt);
  // get a quadrature rule
  const Dune::QuadratureRule<ctype,2>& rule = 
    Dune::QuadratureRules<ctype,2>::rule(gt,1);
  
  // do the assembly loop
  typename Dune::QuadratureRule<ctype,2>::const_iterator r;
  for (size_t p=0;p<interface.size();++p) {
    const BoundaryGrid::Quad& qi(interface[p]);
    HexGeometry<2,2,GridInterface> lg(qi);
    for (size_t q=0;q<per_pillar;++q) {
      const BoundaryGrid::Quad& qu(b1[p+per_pillar*q]);
      HexGeometry<2,2,GridType> hex(qu,pgrid_->grid(),dir);
      Dune::FieldMatrix<ctype,1,4> E; // One row
      E = 0;
      for (r = rule.begin(); r != rule.end();++r) {
        ctype detJ = hex.integrationElement(r->position());
        if (detJ < 0)
          continue;
        typename HexGeometry<2,2,GridInterface>::LocalCoordinate loc = 
                                        lg.local(hex.global(r->position()));
        for (int i=0;i<pbasis.size();++i) {
          for (int j=0;j<lbasis.size();++j)
            E[i][j] += pbasis[i].evaluateFunction(r->position())*
                       lbasis[j].evaluateFunction(loc)*detJ*r->weight();
        }
      }
      // and assemble element contributions 
      for (int i=0;i<1;++i) {
	int indexi;
	if (cellFaces_.empty())
	  indexi = getEquationForDof(qu,dir);
	else
	  indexi = getEquationForDof(qu);
	if (indexi > -1) {
	  for (int j=0;j<4;++j) {
	    int indexj = qi.v[j].i;
	    B[indexi][indexj] += E[i][j];
	  }
	}
      }
    }
  }

  return B;
}

template<class GridInterface>
int MortarHelper<GridInterface>::getEquationForDof(const BoundaryGrid::Quad& quad) {
  ASSERT(!cellFaces_.empty());
  return quad.globalFaceIndex;
}

template<class GridInterface>
int MortarHelper<GridInterface>::getEquationForDof(const BoundaryGrid::Quad& quad, int dir)
{
  if (!cellFaces_.empty()) {
    std::cout << "Warning! Mapping from quad to global face index exists,\n"
	      << "so you should call getEquationForDof(quad) instead.\n"
	      << "Calling getEquationForDof(quad)...\n";
    return getEquationForDof(quad);
  }
 
  GridType gv = pgrid_->grid();
  LeafFaceIterator itFaceStart = gv.leafView().template begin<dim-1>();
  LeafFaceIterator itFaceEnd   = gv.leafView().template end<dim-1>();

  int vertexIdx[4] = {quad.v[0].i, 
		      quad.v[1].i,
		      quad.v[2].i,
		      quad.v[3].i};

  // Find global vertex coordinates
  std::vector<GlobalCoordinate> vertices;
  vertices.push_back(gv.vertexPosition(vertexIdx[0]));
  vertices.push_back(gv.vertexPosition(vertexIdx[1]));
  if (dco_) {
    vertices.push_back(gv.vertexPosition(vertexIdx[3]));
    vertices.push_back(gv.vertexPosition(vertexIdx[2]));
  }
  else {
    vertices.push_back(gv.vertexPosition(vertexIdx[2]));
    vertices.push_back(gv.vertexPosition(vertexIdx[3]));
  }

  GlobalCoordinate faceCentroid = centroid(vertices,dir);

  typedef typename GridInterface::CellIterator CI;
  typedef typename CI::FaceIterator FI;
  
  for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
      if ( (f->boundary()) &&  (isOnPoint(f->centroid(), faceCentroid)) ) {
	return f->index();
      }
    }
  }

  std::cerr << "Global index of current quad not found" << std::endl;
  return -1;
}

template<class GridInterface>
typename MortarHelper<GridInterface>::GlobalCoordinate 
MortarHelper<GridInterface>::centroid(std::vector<MortarHelper<GridInterface>::GlobalCoordinate> 
				      vertices, int dir)
{
  // Calculates 3D centroid of a quad defined by vertices. 
  // Assumes that all vertices lie in the same plane (XY, XZ or YZ)

  int coord1, coord2;
  GlobalCoordinate result(0.0);
  switch (dir) {
  case 0: // x fixed
    coord1 = 1;
    coord2 = 2;
    result[0] = vertices[0][0];
    break;
  case 1: // y fixed
    coord1 = 0;
    coord2 = 2;
    result[1] = vertices[0][1];
    break;
  case 2: // z fixed
    coord1 = 0;
    coord2 = 1;
    result[2] = vertices[0][2];
    break;
  default:
    std::cerr << "dir must be either 0(x), 1(y) or 2(z)" << std::endl;
    break;
  }

  double signedArea = 0.0;
  double x0 = 0.0; // Current vertex X (coord1)
  double y0 = 0.0; // Current vertex Y (coord2)
  double x1 = 0.0; // Next vertex X
  double y1 = 0.0; // Next vertex Y  
  double a = 0.0;  // Partial signed area
  int vertexCount = vertices.size();

  // For all vertices except last
  int i=0;
  for (i=0; i<vertexCount-1; ++i)
    {
      x0 = vertices[i][coord1];
      y0 = vertices[i][coord2];
      x1 = vertices[i+1][coord1];
      y1 = vertices[i+1][coord2];
      a = x0*y1 - x1*y0;
      signedArea += a;
      result[coord1] += (x0 + x1)*a;
      result[coord2] += (y0 + y1)*a;
    }

  // Do last vertex
  x0 = vertices[i][coord1];
  y0 = vertices[i][coord2];
  x1 = vertices[0][coord1];
  y1 = vertices[0][coord2];
  a = x0*y1 - x1*y0;
  signedArea += a;
  result[coord1] += (x0 + x1)*a;
  result[coord2] += (y0 + y1)*a;

  signedArea *= 0.5;
  result[coord1] /= (6*signedArea);
  result[coord2] /= (6*signedArea);

  return result;
}


#endif // MORTAR_HPP
