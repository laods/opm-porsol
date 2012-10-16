// This header includes createPeriodic functions analogous to those in 
// dune/porsol/common/PeriodicHelpers.hpp, but for the Mortar case.
// The only change is that it does not search for periodic partners
// in the X og Y directions.

#include <dune/porsol/common/PeriodicHelpers.hpp>

namespace Dune 
{
  
  template <class BCs, class GridInterface>
  void createPeriodicMortar(BCs& fbcs,
			    const GridInterface& g,
			    const array<FlowBC, 2*GridInterface::Dimension>& fconditions,
			    double spatial_tolerance = 1e-6)
  {
    BOOST_STATIC_ASSERT(BCs::HasFlowConds);
    BOOST_STATIC_ASSERT(!BCs::HasSatConds);
    // Check the conditions given.
    for (int i = 0; i < GridInterface::Dimension; ++i) {
      if (fconditions[2*i].isPeriodic()) {
	ASSERT(fconditions[2*i + 1].isPeriodic());
	ASSERT(fconditions[2*i].pressureDifference() == -fconditions[2*i + 1].pressureDifference());
      }
    }
    std::vector<BoundaryFaceInfo> bfinfo;
    array<double, 6> side_areas;
    createPeriodicMortarImpl(fbcs, bfinfo, side_areas, g, 
		       extractPeriodic(fconditions), spatial_tolerance);
    storeFlowCond(fbcs, bfinfo, fconditions, side_areas);
  }


  template <class BCs, class GridInterface>
  void createPeriodicMortarImpl(BCs& fbcs,
			  std::vector<BoundaryFaceInfo>& bfinfo,
			  array<double, 6>& side_areas,
			  const GridInterface& g,
			  const array<bool, 2*GridInterface::Dimension>& is_periodic,
			  double spatial_tolerance = 1e-6)
  {
    // Fake a variable s.t. we don't search for periodic partners in X/Y direction
    array<bool, 2*GridInterface::Dimension> is_periodic_fake = is_periodic;
    for (int can_pos = 0; can_pos < 4; ++can_pos) {
      is_periodic_fake[can_pos] = false;
    }

    findPeriodicPartners(bfinfo, side_areas, g.grid().leafView(), 
			 is_periodic_fake, spatial_tolerance);

    int num_bdy = bfinfo.size();
    // This will likely change with boundarySegmentIndex() instead of boundaryId():
    int max_bid = num_bdy;

    // Resize the conditions object. We clear it first to make sure it's all defaults.
    fbcs.clear();
    fbcs.resize(max_bid + 1);

    // Now we have all the info, and the object to store it in. So we store it...
    for (int i = 0; i < num_bdy; ++i) {
      int bid1 = bfinfo[i].bid;
      int bid2 = bfinfo[i].partner_bid;
      if (bid1 < bid2) {
	// If there is no periodic partner, bid2 will be zero, so we will not end up here.
	fbcs.setPeriodicPartners(bid1, bid2);
      }
      fbcs.setCanonicalBoundaryId(bid1, bfinfo[i].canon_pos + 1);
    }

  }


} // Namespace Dune
