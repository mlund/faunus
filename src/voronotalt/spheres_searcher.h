#ifndef VORONOTALT_SPHERES_SEARCHER_H_
#define VORONOTALT_SPHERES_SEARCHER_H_

#include <vector>

#include "basic_types_and_functions.h"

namespace voronotalt
{

class SpheresSearcher
{
public:
	explicit SpheresSearcher(const std::vector<SimpleSphere>& spheres) : spheres_(spheres), box_size_(FLOATCONST(1.0))
	{
		for(UnsignedInt i=0;i<spheres_.size();i++)
		{
			const SimpleSphere& s=spheres_[i];
			box_size_=std::max(box_size_, s.r*FLOATCONST(2.0)+FLOATCONST(0.25));
		}

		for(UnsignedInt i=0;i<spheres_.size();i++)
		{
			const GridPoint gp(spheres_[i], box_size_);
			if(i==0)
			{
				grid_offset_=gp;
				grid_size_=gp;
			}
			else
			{
				grid_offset_.x=std::min(grid_offset_.x, gp.x);
				grid_offset_.y=std::min(grid_offset_.y, gp.y);
				grid_offset_.z=std::min(grid_offset_.z, gp.z);
				grid_size_.x=std::max(grid_size_.x, gp.x);
				grid_size_.y=std::max(grid_size_.y, gp.y);
				grid_size_.z=std::max(grid_size_.z, gp.z);
			}
		}

		grid_size_.x=grid_size_.x-grid_offset_.x+1;
		grid_size_.y=grid_size_.y-grid_offset_.y+1;
		grid_size_.z=grid_size_.z-grid_offset_.z+1;

		map_of_boxes_.resize(grid_size_.x*grid_size_.y*grid_size_.z, -1);

		for(UnsignedInt i=0;i<spheres_.size();i++)
		{
			const GridPoint gp(spheres_[i], box_size_, grid_offset_);
			const int index=gp.index(grid_size_);
			const int box_id=map_of_boxes_[index];
			if(box_id<0)
			{
				map_of_boxes_[index]=static_cast<int>(boxes_.size());
				boxes_.push_back(std::vector<UnsignedInt>(1, i));
			}
			else
			{
				boxes_[box_id].push_back(i);
			}
		}
	}

	bool find_colliding_ids(const UnsignedInt& central_id, std::vector<UnsignedInt>& colliding_ids, const bool discard_hidden, int& exclusion_status) const
	{
		colliding_ids.clear();
		exclusion_status=0;
		if(central_id<spheres_.size())
		{
			const SimpleSphere& central_sphere=spheres_[central_id];
			const GridPoint gp(central_sphere, box_size_, grid_offset_);
			GridPoint dgp=gp;
			for(int dx=-1;dx<=1;dx++)
			{
				dgp.x=gp.x+dx;
				for(int dy=-1;dy<=1;dy++)
				{
					dgp.y=gp.y+dy;
					for(int dz=-1;dz<=1;dz++)
					{
						dgp.z=gp.z+dz;
						const int index=dgp.index(grid_size_);
						if(index>=0)
						{
							const int box_id=map_of_boxes_[index];
							if(box_id>=0)
							{
								const std::vector<UnsignedInt>& ids=boxes_[box_id];
								for(UnsignedInt i=0;i<ids.size();i++)
								{
									const UnsignedInt id=ids[i];
									const SimpleSphere& candidate_sphere=spheres_[id];
									if(id!=central_id && sphere_intersects_sphere(central_sphere, candidate_sphere))
									{
										if(discard_hidden
												&& sphere_contains_sphere(candidate_sphere, central_sphere)
												&& (!sphere_equals_sphere(candidate_sphere, central_sphere) || central_id>id))
										{
											colliding_ids.clear();
											exclusion_status=1;
											return false;
										}
										else if(!discard_hidden
												|| !sphere_contains_sphere(central_sphere, candidate_sphere))
										{
											colliding_ids.push_back(id);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return (!colliding_ids.empty());
	}

private:
	struct GridPoint
	{
		int x;
		int y;
		int z;

		GridPoint() : x(0), y(0), z(0)
		{
		}

		GridPoint(const SimpleSphere& s, const Float grid_step)
		{
			init(s, grid_step);
		}

		GridPoint(const SimpleSphere& s, const Float grid_step, const GridPoint& grid_offset)
		{
			init(s, grid_step, grid_offset);
		}

		void init(const SimpleSphere& s, const Float grid_step)
		{
			x=static_cast<int>(s.p.x/grid_step);
			y=static_cast<int>(s.p.y/grid_step);
			z=static_cast<int>(s.p.z/grid_step);
		}

		void init(const SimpleSphere& s, const Float grid_step, const GridPoint& grid_offset)
		{
			x=static_cast<int>(s.p.x/grid_step)-grid_offset.x;
			y=static_cast<int>(s.p.y/grid_step)-grid_offset.y;
			z=static_cast<int>(s.p.z/grid_step)-grid_offset.z;
		}

		int index(const GridPoint& grid_size) const
		{
			return ((x>=0 && y>=0 && z>=0 && x<grid_size.x && y<grid_size.y &&z<grid_size.z) ? (z*grid_size.x*grid_size.y+y*grid_size.x+x) : (-1));
		}

		bool operator<(const GridPoint& gp) const
		{
			return (x<gp.x || (x==gp.x && y<gp.y) || (x==gp.x && y==gp.y && z<gp.z));
		}
	};

	const std::vector<SimpleSphere>& spheres_;
	GridPoint grid_offset_;
	GridPoint grid_size_;
	std::vector<int> map_of_boxes_;
	std::vector< std::vector<UnsignedInt> > boxes_;
	Float box_size_;
};

}

#endif /* VORONOTALT_SPHERES_SEARCHER_H_ */
