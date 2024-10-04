#ifndef VORONOTALT_PERIODIC_BOX_H_
#define VORONOTALT_PERIODIC_BOX_H_

#include <vector>

#include "basic_types_and_functions.h"

namespace voronotalt
{

class PeriodicBox
{
public:
	PeriodicBox() noexcept : enabled_(false)
	{
	}

	bool enabled() const noexcept
	{
		return enabled_;
	}

	bool equals(const PeriodicBox& pb) const noexcept
	{
		return (enabled_==pb.enabled_ && point_equals_point(shift_direction_a_, pb.shift_direction_a_) && point_equals_point(shift_direction_b_, pb.shift_direction_b_) && point_equals_point(shift_direction_c_, pb.shift_direction_c_));
	}

	const SimpleSphere shift_by_weighted_directions(const SimpleSphere& o, const Float wa, const Float wb, const Float wc) const
	{
		return SimpleSphere(
				o.p.x+shift_direction_a_.x*wa+shift_direction_b_.x*wb+shift_direction_c_.x*wc,
				o.p.y+shift_direction_a_.y*wa+shift_direction_b_.y*wb+shift_direction_c_.y*wc,
				o.p.z+shift_direction_a_.z*wa+shift_direction_b_.z*wb+shift_direction_c_.z*wc,
				o.r);
	}

	static PeriodicBox create_periodic_box_from_corners(const std::vector<SimplePoint>& periodic_box_corners) noexcept
	{
		PeriodicBox pb;
		if(periodic_box_corners.size()>=2)
		{
			pb.enabled_=true;
			SimplePoint corner_min=periodic_box_corners[0];
			SimplePoint corner_max=periodic_box_corners[0];
			for(UnsignedInt i=1;i<periodic_box_corners.size();i++)
			{
				const SimplePoint& corner=periodic_box_corners[i];
				corner_min.x=std::min(corner_min.x, corner.x);
				corner_min.y=std::min(corner_min.y, corner.y);
				corner_min.z=std::min(corner_min.z, corner.z);
				corner_max.x=std::max(corner_max.x, corner.x);
				corner_max.y=std::max(corner_max.y, corner.y);
				corner_max.z=std::max(corner_max.z, corner.z);
			}
			pb.shift_direction_a_.x=corner_max.x-corner_min.x;
			pb.shift_direction_b_.y=corner_max.y-corner_min.y;
			pb.shift_direction_c_.z=corner_max.z-corner_min.z;
		}
		return pb;
	}

	static PeriodicBox create_periodic_box_from_shift_directions(const std::vector<SimplePoint>& shift_directions) noexcept
	{
		PeriodicBox pb;
		if(shift_directions.size()==3)
		{
			pb.enabled_=true;
			pb.shift_direction_a_=shift_directions[0];
			pb.shift_direction_b_=shift_directions[1];
			pb.shift_direction_c_=shift_directions[2];
		}
		return pb;
	}

	static PeriodicBox create_periodic_box_from_shift_directions_or_from_corners(const std::vector<SimplePoint>& shift_directions, const std::vector<SimplePoint>& periodic_box_corners) noexcept
	{
		PeriodicBox pb=create_periodic_box_from_shift_directions(shift_directions);
		if(!pb.enabled())
		{
			pb=create_periodic_box_from_corners(periodic_box_corners);
		}
		return pb;
	}

private:
	bool enabled_;
	SimplePoint shift_direction_a_;
	SimplePoint shift_direction_b_;
	SimplePoint shift_direction_c_;
};

}

#endif /* VORONOTALT_PERIODIC_BOX_H_ */
