
#ifndef MORTAR_HPP
#define MORTAR_HPP


#include <vector>

#include <dune/grid/CpGrid.hpp>
#include <dune/porsol/mortar/boundarygrid.cc>
#include <dune/porsol/mortar/matrixops.hpp>
#include <dune/porsol/mortar/shapefunctions.hh>

//! \brief A sparse matrix holding our operator
typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > Matrix;

//! \brief A vector holding our RHS
typedef Dune::BlockVector<Dune::FieldVector<double,1> > Vector;

//! \brief For storing matrix adjacency/sparsity patterns
typedef std::vector< std::set<int> > AdjacencyPattern;

// Helper class for mortar methods
template<class GridType>
class MortarHelper 
{
public:

  static const int dim = GridType::dimension;
  typedef typename GridType::LeafGridView::ctype ctype;
  typedef typename GridType::LeafGridView::IndexSet LeafIndexSet;
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator LeafVertexIterator;
  typedef typename GridType::LeafGridView::template Codim<0>::Iterator LeafIterator;
  typedef typename GridType::LeafGridView::template Codim<1>::Geometry::GlobalCoordinate 
  GlobalCoordinate;
  typedef typename GridType::LeafGridView::template Codim<1>::Geometry::LocalCoordinate 
  LocalCoordinate;

  MortarHelper() {};

  MortarHelper(const GridType& gv)
    : pgv_(&gv) 
  {
    findMinMax();
    find_n();
    tol_ = 1e-8;
  };

  MortarHelper(const GridType& gv, std::vector<ctype> min,
	       std::vector<ctype> max, int n1, int n2_, double tol = 1e-8)
    : pgv_(&gv), min_(min), max_(max), n1_(n1), n2_(n2), tol_(tol) {};

  std::vector<double> min();
  std::vector<double> max();
  int n1();
  int n2();
  double tol_;
  void setMinMax(std::vector<double> min, std::vector<double> max);
  void findMinMax();
  void set_n(int n1, int n2);
  void find_n();
  
  // Test function for debugging
  // Print vertices on face quads with global index and coord
  void printFace(int face);

  void periodicBCsMortar();

private:
  std::vector<double> min_;
  std::vector<double> max_;
  int n1_;
  int n2_;
  const GridType* pgv_;
  
  enum SIDE {
    LEFT,
    RIGHT
  };
  enum Direction { NONE = 0, X = 1, Y = 2, Z = 4,
                     XY = 1+2, XZ = 1+4, YZ = 2+4, 
                    XYZ = 1+2+4 };

  std::vector<BoundaryGrid> master;
  std::vector< std::vector<BoundaryGrid::Vertex> > slave;
  std::vector<Matrix> L;

  std::vector<BoundaryGrid::Vertex> extractFace(Direction dir, ctype coord);
  BoundaryGrid extractMasterFace(Direction dir, ctype coord, SIDE side=LEFT, bool dc=false);
  bool isOnPlane(Direction plane, GlobalCoordinate coord, ctype value);
  bool isOnLine(Direction dir, GlobalCoordinate coord, ctype x, ctype y);
  bool isOnPoint(GlobalCoordinate coord, GlobalCoordinate point);

  Matrix findLMatrixMortar(const BoundaryGrid& b1, const BoundaryGrid& interface, int dir);

};


template<class GridType>
std::vector<double> MortarHelper<GridType>::min() {
  return min_;
}

template<class GridType>
std::vector<double> MortarHelper<GridType>::max() {
  return max_;
}

template<class GridType>
int MortarHelper<GridType>::n1() {
  return n1_;
}

template<class GridType>
int MortarHelper<GridType>::n2() {
  return n2_;
}

template<class GridType>
void MortarHelper<GridType>::setMinMax(std::vector<double> min, std::vector<double> max) {
  min_ = min;
  max_ = max;
}

template<class GridType>
void MortarHelper<GridType>::set_n(int n1, int n2) {
  n1_ = n1;
  n2_ = n2;
}

template<class GridType>
void MortarHelper<GridType>::findMinMax() {
  // Based on ElasticityUpscale::findBoundaries()
  double big = 1e5;
  if (min_.empty()) min_.resize(3,big);
  else min_[0] = min_[1] = min_[2] = big;
  if (max_.empty()) max_.resize(3,-big);
  else max_[0] = max_[1] = max_[2] = -big;

  const LeafVertexIterator itend = pgv_->leafView().template end<dim>();
  // iterate over vertices and find slaves
  LeafVertexIterator start = pgv_->leafView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    for (int i=0;i<3;++i) {
      min_[i] = std::min(min_[i],it->geometry().corner(0)[i]);
      max_[i] = std::max(max_[i],it->geometry().corner(0)[i]);
    }
  }
}

template<class GridType>
void MortarHelper<GridType>::find_n() {
  n1_ = pgv_->logicalCartesianSize()[0];
  n2_ = pgv_->logicalCartesianSize()[1];
}

template<class GridType>
std::vector<BoundaryGrid::Vertex> MortarHelper<GridType>::extractFace(Direction dir, ctype coord) {
  // Based on ElasticityUpscale::extractFaces
  std::vector<BoundaryGrid::Vertex> result;
  const LeafIndexSet& set = pgv_->leafView().indexSet();
  const LeafVertexIterator itend = pgv_->leafView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(*pgv_);
  // iterate over vertices and find slaves
  LeafVertexIterator start = pgv_->leafView().template begin<dim>();
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

template<class GridType>
BoundaryGrid MortarHelper<GridType>::extractMasterFace(Direction dir, 
						       ctype coord,
						       SIDE side,
						       bool dc) 
{
  // Based on ElasticityUpscale::extractMasterFace
  static const int V1[3][4] = {{0,2,4,6},
                               {0,1,4,5},
                               {0,1,2,3}};
  static const int V2[3][4] = {{1,3,5,7},
                               {2,3,6,7},
                               {4,5,6,7}};
  const LeafIndexSet& set = pgv_->leafView().indexSet();
  const LeafVertexIterator itend = pgv_->leafView().template end<dim>();

  // make a mapper for codim dim entities in the leaf grid 
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,
                                            Dune::MCMGVertexLayout> mapper(*pgv_);
  LeafVertexIterator start = pgv_->leafView().template begin<dim>();
  LeafIterator cellend = pgv_->leafView().template end<0>();
  int c = 0;
  int i = log2(dir);
  BoundaryGrid result;
  for (LeafIterator cell  = pgv_->leafView().template begin<0>(); 
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
      result.add(q);
    }
  }

  return result;
}

template<class GridType>
bool MortarHelper<GridType>::isOnPlane(Direction plane,
                                            GlobalCoordinate coord,
                                            ctype value)
{
  if (plane < X || plane > Z)
    return false;
  int p = log2(plane);
  ctype delta = fabs(value-coord[p]);
  return delta < tol_;
}

template<class GridType>
bool MortarHelper<GridType>::isOnLine(Direction dir,
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

template<class GridType>
bool MortarHelper<GridType>::isOnPoint(GlobalCoordinate coord,
                                            GlobalCoordinate point)
{
  GlobalCoordinate delta = point-coord;
  return delta.one_norm() < tol_;
}

template<class GridType>
void MortarHelper<GridType>::printFace(int face) {
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
}

template<class GridType>
void MortarHelper<GridType>::periodicBCsMortar() {

  // Based on ElasticityUpscale::periodicBCsMortar()
  // But not MPC part (that is only step 3-6)

  // Step 3: extracts and establishes a quad grid 
  //         for the left/right/front/back sides
  master.push_back(extractMasterFace(X, min_[0], LEFT, true));
  master.push_back(extractMasterFace(X, max_[0], RIGHT, true));
  master.push_back(extractMasterFace(Y, min_[1], LEFT, true));
  master.push_back(extractMasterFace(Y, max_[1], RIGHT, true));

  // Step 4: Establishes grids for the dual dofs
  BoundaryGrid::FaceCoord fmin, fmax;
  fmin[0] = min_[1]; fmin[1] = min_[2];
  fmax[0] = max_[1]; fmax[1] = max_[2];
  BoundaryGrid lambdax = BoundaryGrid::uniform(fmin, fmax, n2_, 1, true);
  
  fmin[0] = min_[0]; fmin[1] = min_[2];
  fmax[0] = max_[0]; fmax[1] = max_[2];
  BoundaryGrid lambday = BoundaryGrid::uniform(fmin, fmax, n1_, 1, true);

  /*
  // Step 5: Calculates the coupling matrix L1 between left/right sides
  Matrix L1_left  = findLMatrixMortar(master[0], lambdax, 0);
  Matrix L1_right = findLMatrixMortar(master[1], lambdax, 0);
  L.push_back(MatrixOps::Axpy(L1_left, L1_right, -1));

  // Step 6: Calculates the coupling matrix L1 between front/back sides
  Matrix L2_left  = findLMatrixMortar(master[2], lambday, 1);
  Matrix L2_right = findLMatrixMortar(master[3], lambday, 1);
  L.push_back(MatrixOps::Axpy(L2_left, L2_right, -1));
  */
}


template<class GridType>
Matrix MortarHelper<GridType>::findLMatrixMortar(const BoundaryGrid& b1,
						 const BoundaryGrid& interface,
						 int dir)
{
  std::vector< std::set<int> > adj;
  //adj.resize(A.getEqns());

  // process pillar by pillar
  size_t per_pillar = b1.size()/interface.size();
  for (size_t p=0;p<interface.size();++p) {
    for (size_t q=0;q<per_pillar;++q) {
      for (size_t i=0;i<1;++i) { // ?? i<4 or i<1 ??
        //for (size_t d=0;d<3;++d) { // We only have one DOF per vertex
	//MPC* mpc = A.getMPC(b1[p*per_pillar+q].v[i].i,d);

	// A given DOF here is never a MPC since the DOF is at face centers, not vertices!
	// Check for MPC not nescessary (as it is in elasticity upscale)
	
	// Own code:
	/*
	int dof = getEquationForDof(face);
	for (int j=0;j<4;++j) {
	  adj[dof].insert(interface[p].v[j].i); 
	  // No need for multiplying with 3 (only one DOF per vertex)
	}
	*/
	  
	// Original implementation from elasticity:
	/*
          if (mpc) {
            for (int n=0;n<mpc->getNoMaster();++n) {
              int dof = A.getEquationForDof(mpc->getMaster(n).node,d);
              if (dof > -1) {
                for (int j=0;j<4;++j)
                  adj[dof].insert(3*interface[p].v[j].i+d);
              }
            }
          } else {
            int dof = A.getEquationForDof(b1[p*per_pillar+q].v[i].i,d);
            if (dof > -1) {
              for (int j=0;j<4;++j)
                adj[dof].insert(3*interface[p].v[j].i+d);
            }
          }
	*/
	//}
      }
    }
  }

  /*
  Matrix B;
  MatrixOps::fromAdjacency(B,adj,A.getEqns(),interface.totalNodes()); // Fix getEqns()

  // get a set of P0 shape functions for the face pressures
  P0ShapeFunctionSet<ctype,ctype,2> pbasis = P0ShapeFunctionSet<ctype,ctype,2>::instance();
  // get a set of P1 shape functions for multipliers
  P1ShapeFunctionSet<ctype,ctype,2> lbasis = P1ShapeFunctionSet<ctype,ctype,2>::instance();
  // get a reference element
  Dune::GeometryType gt;
  gt.makeCube(2);
  const Dune::template GenericReferenceElement<ctype,2> &ref =
    Dune::GenericReferenceElements<ctype,2>::general(gt);
  // get a quadrature rule
  const Dune::QuadratureRule<ctype,2>& rule = 
                        Dune::QuadratureRules<ctype,2>::rule(gt,2);

  // do the assembly loop
  typename Dune::QuadratureRule<ctype,2>::const_iterator r;
  for (size_t p=0;p<interface.size();++p) {
    const BoundaryGrid::Quad& qi(interface[p]);
    HexGeometry<2,2,GridType> lg(qi);
    for (size_t q=0;q<per_pillar;++q) {
      const BoundaryGrid::Quad& qu(b1[p*per_pillar+q]);
      HexGeometry<2,2,GridType> hex(qu,&pgv_,dir);
      Dune::FieldMatrix<ctype,1,4> E; // One row
      E = 0;
      for (r = rule.begin(); r != rule.end();++r) {
        ctype detJ = hex.integrationElement(r->position());
        if (detJ < 0)
          continue;
        typename HexGeometry<2,2,GridType>::LocalCoordinate loc = 
                                        lg.local(hex.global(r->position()));
        for (int i=0;i<pbasis.size();++i) {
          for (int j=0;j<lbasis.size();++j)
            E[i][j] += pbasis[i].evaluateFunction(r->position())*
                       lbasis[j].evaluateFunction(loc)*detJ*r->weight();
        }
      }
      // and assemble element contributions
      
      // Own code:
      for (int i=0;i<1;++i) {
	// No need to check for MPC
	int indexi = A.getEquationForDof(qu.v[i].i,d); // Fix
	if (indexi > -1) {
	  for (int j=0;j<4;++j) {
	    int indexj = qi.v[j].i;
	    B[indexi][indexj] += E[i][j];
	  }
	}
      }

      // Old:
      /*
      for (int d=0;d<3;++d) {
        for (int i=0;i<4;++i) {
          MPC* mpc = A.getMPC(qu.v[i].i,d);
          if (mpc) {
            for (int n=0;n<mpc->getNoMaster();++n) {
              int indexi = A.getEquationForDof(mpc->getMaster(n).node,d);
              if (indexi > -1) {
                for (int j=0;j<4;++j) {
                  int indexj = qi.v[j].i*3+d;
                  B[indexi][indexj] += E[i][j];
                }
              }
            }
          } else {
            int indexi = A.getEquationForDof(qu.v[i].i,d);
            if (indexi > -1) {
              for (int j=0;j<4;++j) {
                int indexj = qi.v[j].i*3+d;
                B[indexi][indexj] += E[i][j];
              }
            }
          }
        }
      }
      
    }
  }

  return B;
  */
}

#endif // MORTAR_HPP
