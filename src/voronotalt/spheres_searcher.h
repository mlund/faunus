#ifndef VORONOTALT_SPHERES_SEARCHER_H_
#define VORONOTALT_SPHERES_SEARCHER_H_

#include <vector>
#include <algorithm>

#include "basic_types_and_functions.h"

namespace voronotalt
{

class SpheresSearcher
{
public:
	SpheresSearcher() noexcept
	{
	}

	explicit SpheresSearcher(const std::vector<SimpleSphere>& spheres) noexcept
	{
		init(spheres);
	}

	void init(const std::vector<SimpleSphere>& spheres) noexcept
	{
		spheres_=spheres;
		grid_parameters_=GridParameters(spheres_);
		init_boxes();
	}

	bool update(const std::vector<SimpleSphere>& spheres, const std::vector<UnsignedInt>& ids_to_update) noexcept
	{
		if(ids_to_update.empty())
		{
			return false;
		}

		if(ids_to_update.size()>size_threshold_for_full_rebuild())
		{
			init(spheres);
			return true;
		}

		for(UnsignedInt i=0;i<ids_to_update.size();i++)
		{
			if(ids_to_update[i]>=spheres.size() || !update_sphere(ids_to_update[i], spheres[ids_to_update[i]]))
			{
				init(spheres);
				return true;
			}
		}

		return true;
	}

	void assign(const SpheresSearcher& obj) noexcept
	{
		grid_parameters_=obj.grid_parameters_;

		spheres_.resize(obj.spheres_.size());
		map_of_boxes_.resize(obj.map_of_boxes_.size());
		boxes_.resize(obj.boxes_.size());

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<obj.spheres_.size();i++)
			{
				spheres_[i]=obj.spheres_[i];
			}
		}

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<obj.map_of_boxes_.size();i++)
			{
				map_of_boxes_[i]=obj.map_of_boxes_[i];
			}
		}

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<obj.boxes_.size();i++)
			{
				boxes_[i]=obj.boxes_[i];
			}
		}
	}

	const std::vector<SimpleSphere>& searchable_spheres() const noexcept
	{
		return spheres_;
	}

	bool find_colliding_ids(const UnsignedInt& central_id, std::vector<ValuedID>& colliding_ids, const bool discard_hidden, int& exclusion_status) const noexcept
	{
		colliding_ids.clear();
		exclusion_status=0;
		if(central_id<spheres_.size())
		{
			const SimpleSphere& central_sphere=spheres_[central_id];
			const GridPoint gp(central_sphere, grid_parameters_.box_size, grid_parameters_.grid_offset);
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
						const int index=dgp.index(grid_parameters_.grid_size);
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
											colliding_ids.push_back(ValuedID(distance_to_center_of_intersection_circle_of_two_spheres(central_sphere, candidate_sphere), id));
										}
									}
								}
							}
						}
					}
				}
			}
		}
		if(!colliding_ids.empty())
		{
			std::sort(colliding_ids.begin(), colliding_ids.end());
		}
		return (!colliding_ids.empty());
	}

private:
	struct GridPoint
	{
		int x;
		int y;
		int z;

		GridPoint()  noexcept: x(0), y(0), z(0)
		{
		}

		GridPoint(const SimpleSphere& s, const Float grid_step) noexcept
		{
			init(s, grid_step);
		}

		GridPoint(const SimpleSphere& s, const Float grid_step, const GridPoint& grid_offset) noexcept
		{
			init(s, grid_step, grid_offset);
		}

		void init(const SimpleSphere& s, const Float grid_step) noexcept
		{
			x=static_cast<int>(s.p.x/grid_step);
			y=static_cast<int>(s.p.y/grid_step);
			z=static_cast<int>(s.p.z/grid_step);
		}

		void init(const SimpleSphere& s, const Float grid_step, const GridPoint& grid_offset) noexcept
		{
			x=static_cast<int>(s.p.x/grid_step)-grid_offset.x;
			y=static_cast<int>(s.p.y/grid_step)-grid_offset.y;
			z=static_cast<int>(s.p.z/grid_step)-grid_offset.z;
		}

		int index(const GridPoint& grid_size) const noexcept
		{
			return ((x>=0 && y>=0 && z>=0 && x<grid_size.x && y<grid_size.y && z<grid_size.z) ? (z*grid_size.x*grid_size.y+y*grid_size.x+x) : (-1));
		}

		bool operator<(const GridPoint& gp) const noexcept
		{
			return (x<gp.x || (x==gp.x && y<gp.y) || (x==gp.x && y==gp.y && z<gp.z));
		}

		bool operator==(const GridPoint& gp) const noexcept
		{
			return (x==gp.x && y==gp.y && z==gp.z);
		}
	};

	struct GridParameters
	{
		GridPoint grid_offset;
		GridPoint grid_size;
		Float box_size;
		int padding_;

		GridParameters() noexcept : box_size(FLOATCONST(1.0)), padding_(1)
		{
			grid_size.x=1;
			grid_size.y=1;
			grid_size.z=1;
		}

		explicit GridParameters(const std::vector<SimpleSphere>& spheres) noexcept : box_size(FLOATCONST(1.0)), padding_(1)
		{
			init(spheres);
		}

		GridParameters(const std::vector<SimpleSphere>& spheres, const int padding) noexcept : box_size(FLOATCONST(1.0)), padding_(padding)
		{
			init(spheres);
		}

		void init(const std::vector<SimpleSphere>& spheres) noexcept
		{
			box_size=FLOATCONST(1.0);
			padding_=std::max(0, padding_);

			for(UnsignedInt i=0;i<spheres.size();i++)
			{
				const SimpleSphere& s=spheres[i];
				box_size=std::max(box_size, s.r*FLOATCONST(2.0)+FLOATCONST(0.25));
			}

			for(UnsignedInt i=0;i<spheres.size();i++)
			{
				const GridPoint gp(spheres[i], box_size);
				if(i==0)
				{
					grid_offset=gp;
					grid_size=gp;
				}
				else
				{
					grid_offset.x=std::min(grid_offset.x, gp.x-padding_);
					grid_offset.y=std::min(grid_offset.y, gp.y-padding_);
					grid_offset.z=std::min(grid_offset.z, gp.z-padding_);
					grid_size.x=std::max(grid_size.x, gp.x+padding_);
					grid_size.y=std::max(grid_size.y, gp.y+padding_);
					grid_size.z=std::max(grid_size.z, gp.z+padding_);
				}
			}

			grid_size.x=grid_size.x-grid_offset.x+1;
			grid_size.y=grid_size.y-grid_offset.y+1;
			grid_size.z=grid_size.z-grid_offset.z+1;
		}

		bool operator==(const GridParameters& gp) const noexcept
		{
			return (box_size==gp.box_size && grid_offset==gp.grid_offset && grid_size==gp.grid_size);
		}
	};

	UnsignedInt size_threshold_for_full_rebuild() const noexcept
	{
		return static_cast<UnsignedInt>(spheres_.size()/2);
	}

	void init_boxes() noexcept
	{
		map_of_boxes_.clear();
		boxes_.clear();

		map_of_boxes_.resize(grid_parameters_.grid_size.x*grid_parameters_.grid_size.y*grid_parameters_.grid_size.z, -1);

		for(UnsignedInt i=0;i<spheres_.size();i++)
		{
			const GridPoint gp(spheres_[i], grid_parameters_.box_size, grid_parameters_.grid_offset);
			const int index=gp.index(grid_parameters_.grid_size);
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

	bool update_sphere(const UnsignedInt& sphere_id, const SimpleSphere& moved_sphere) noexcept
	{
		if(sphere_id<spheres_.size() && moved_sphere.r<=spheres_[sphere_id].r)
		{
			const GridPoint gp1(moved_sphere, grid_parameters_.box_size, grid_parameters_.grid_offset);
			const int index1=gp1.index(grid_parameters_.grid_size);
			if(index1>=0)
			{
				const GridPoint gp0(spheres_[sphere_id], grid_parameters_.box_size, grid_parameters_.grid_offset);
				const int index0=gp0.index(grid_parameters_.grid_size);
				if(index0>=0)
				{
					if(index1!=index0)
					{
						{
							const int box_id0=map_of_boxes_[index0];
							if(box_id0<0)
							{
								return false;
							}
							else
							{
								std::vector<UnsignedInt>& v0=boxes_[box_id0];
								std::vector<UnsignedInt>::iterator it=std::lower_bound(v0.begin(), v0.end(), sphere_id);
								if(it!=v0.end() && (*it)==sphere_id)
								{
									v0.erase(it);
								}
							}
						}

						{
							const int box_id1=map_of_boxes_[index1];
							if(box_id1<0)
							{
								map_of_boxes_[index1]=static_cast<int>(boxes_.size());
								boxes_.push_back(std::vector<UnsignedInt>(1, sphere_id));
							}
							else
							{
								std::vector<UnsignedInt>& v1=boxes_[box_id1];
								std::vector<UnsignedInt>::iterator it=std::lower_bound(v1.begin(), v1.end(), sphere_id);
								if(it==v1.end() || (*it)!=sphere_id)
								{
									v1.insert(it, sphere_id);
								}
							}
						}
					}
					spheres_[sphere_id]=moved_sphere;
					return true;
				}
			}
		}
		return false;
	}

	std::vector<SimpleSphere> spheres_;
	GridParameters grid_parameters_;
	std::vector<int> map_of_boxes_;
	std::vector< std::vector<UnsignedInt> > boxes_;
};

}

#endif /* VORONOTALT_SPHERES_SEARCHER_H_ */
