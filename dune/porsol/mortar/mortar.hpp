#include <vector>

#include <dune/grid/CpGrid.hpp>

using namespace std;

// Helper class for mortar methods
template<class GridType, int dim>
class MortarHelper 
{
public:

  typedef typename GridType::LeafGridView::IndexSet LeafIndexSet;
  typedef typename GridType::LeafGridView::template Codim<dim>::Iterator LeafVertexIterator;

  MortarHelper() {};
  MortarHelper(vector<double> min, vector<double> max, int n1, int n2)
    : min_(min), max_(max), n1_(n1), n2_(n2) {};
  vector<double> min();
  vector<double> max();
  int n1();
  int n2();
  void setMinMax(vector<double> min, vector<double> max);
  void findMinMax(GridType gv);
  void set_n(int n1, int n2);
  void find_n(GridType gv);
private:
  vector<double> min_;
  vector<double> max_;
  int n1_;
  int n2_;
};

/*
template<class GridType, int dim>
MortarHelper<GridType,dim>::MortarHelper() {
  n1_ = -1;
  n2_ = -1;
  vector<double> empty;
  min_ = empty;
  max_ = empty;
}
*/

template<class GridType, int dim>
vector<double> MortarHelper<GridType,dim>::min() {
  return min_;
}

template<class GridType, int dim>
vector<double> MortarHelper<GridType,dim>::max() {
  return max_;
}

template<class GridType, int dim>
int MortarHelper<GridType,dim>::n1() {
  return n1_;
}

template<class GridType, int dim>
int MortarHelper<GridType,dim>::n2() {
  return n2_;
}

template<class GridType, int dim>
void MortarHelper<GridType,dim>::setMinMax(vector<double> min, vector<double> max) {
  min_ = min;
  max_ = max;
}

template<class GridType, int dim>
void MortarHelper<GridType,dim>::set_n(int n1, int n2) {
  n1_ = n1;
  n2_ = n2;
}

template<class GridType, int dim>
void MortarHelper<GridType,dim>::findMinMax(GridType gv) {
  // Based on ElasticityUpscale::findBoundaries()
  double big = 1e5;
  if (min_.empty()) min_.resize(3,big);
  else min_[0] = min_[1] = min_[2] = big;
  if (max_.empty()) max_.resize(3,-big);
  else max_[0] = max_[1] = max_[2] = -big;

  const LeafVertexIterator itend = gv.leafView().template end<dim>();
  // iterate over vertices and find slaves
  LeafVertexIterator start = gv.leafView().template begin<dim>();
  for (LeafVertexIterator it = start; it != itend; ++it) {
    for (int i=0;i<3;++i) {
      min_[i] = std::min(min_[i],it->geometry().corner(0)[i]);
      max_[i] = std::max(max_[i],it->geometry().corner(0)[i]);
    }
  }
}

template<class GridType, int dim>
void MortarHelper<GridType,dim>::find_n(GridType gv) {
  n1_ = gv.logicalCartesianSize()[0];
  n2_ = gv.logicalCartesianSize()[1];
}
