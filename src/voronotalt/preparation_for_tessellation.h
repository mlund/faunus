#ifndef VORONOTALT_PREPARATION_FOR_TESSELLATION_H_
#define VORONOTALT_PREPARATION_FOR_TESSELLATION_H_

#include <algorithm>

#include "spheres_searcher.h"
#include "time_recorder.h"

namespace voronotalt
{

class PreparationForTessellation
{
public:
	struct Result
	{
		std::vector<SimpleSphere> spheres_with_periodic_instances;
		std::vector<UnsignedInt> periodic_links_of_spheres;
		std::vector<int> all_exclusion_statuses;
		std::vector< std::vector<UnsignedInt> > all_colliding_ids;
		std::vector< std::pair<UnsignedInt, UnsignedInt> > relevant_collision_ids;
		UnsignedInt total_input_spheres;
		UnsignedInt total_spheres;
		UnsignedInt total_collisions;

		Result() : total_input_spheres(0), total_spheres(0), total_collisions(0)
		{
		}
	};

	static void prepare_for_tessellation(
			const std::vector<SimpleSphere>& spheres,
			const std::vector<int>& grouping_of_spheres,
			Result& result,
			TimeRecorder& time_recorder)
	{
		prepare_for_tessellation(spheres, grouping_of_spheres, std::vector<SimplePoint>(), result, time_recorder);
	}

	static void prepare_for_tessellation(
			const std::vector<SimpleSphere>& input_spheres,
			const std::vector<int>& grouping_of_spheres,
			const std::vector<SimplePoint>& periodic_box_corners,
			Result& result,
			TimeRecorder& time_recorder)
	{
		time_recorder.reset();

		result=Result();

		result.total_input_spheres=input_spheres.size();

		if(periodic_box_corners.size()>=2)
		{
			SimplePoint corner_a=periodic_box_corners[0];
			SimplePoint corner_b=periodic_box_corners[0];
			for(UnsignedInt i=1;i<periodic_box_corners.size();i++)
			{
				const SimplePoint& corner=periodic_box_corners[i];
				corner_a.x=std::min(corner_a.x, corner.x);
				corner_a.y=std::min(corner_a.y, corner.y);
				corner_a.z=std::min(corner_a.z, corner.z);
				corner_b.x=std::max(corner_b.x, corner.x);
				corner_b.y=std::max(corner_b.y, corner.y);
				corner_b.z=std::max(corner_b.z, corner.z);
			}
			const SimplePoint shift=sub_of_points(corner_b, corner_a);

			Float buffer_distance=SpheresSearcher::calculate_grid_box_size(input_spheres, FLOATCONST(0.1));

			result.spheres_with_periodic_instances=input_spheres;
			result.periodic_links_of_spheres.resize(input_spheres.size());

			for(UnsignedInt i=0;i<input_spheres.size();i++)
			{
				result.periodic_links_of_spheres[i]=i;
				const SimpleSphere& o=input_spheres[i];
				SimpleSphere m=o;
				for(int sx=-1;sx<=1;sx++)
				{
					m.p.x=o.p.x+(shift.x*static_cast<Float>(sx));
					const bool useful_x=((sx<0 && (corner_a.x-m.p.x)<=buffer_distance) || (sx>0 && (m.p.x-corner_b.x)<=buffer_distance));
					for(int sy=-1;sy<=1;sy++)
					{
						m.p.y=o.p.y+(shift.y*static_cast<Float>(sy));
						const bool useful_y=((sy<0 && (corner_a.y-m.p.y)<=buffer_distance) || (sy>0 && (m.p.y-corner_b.y)<=buffer_distance));
						for(int sz=-1;sz<=1;sz++)
						{
							if(sx!=0 || sy!=0 || sz!=0)
							{
								m.p.z=o.p.z+(shift.z*static_cast<Float>(sz));
								const bool useful_z=((sz<0 && (corner_a.z-m.p.z)<=buffer_distance) || (sz>0 && (m.p.z-corner_b.z)<=buffer_distance));
								if(useful_x || useful_y || useful_z)
								{
									result.spheres_with_periodic_instances.push_back(m);
									result.periodic_links_of_spheres.push_back(i);
								}
							}
						}
					}
				}
			}
		}

		const std::vector<SimpleSphere>& spheres=(result.spheres_with_periodic_instances.empty() ? input_spheres : result.spheres_with_periodic_instances);

		const UnsignedInt N=spheres.size();
		result.total_spheres=N;

		SpheresSearcher spheres_searcher(spheres);

		time_recorder.record_elapsed_miliseconds_and_reset("init spheres searcher");

		result.all_exclusion_statuses.resize(N, 0);

		result.all_colliding_ids.resize(N);
		for(UnsignedInt i=0;i<N;i++)
		{
			result.all_colliding_ids[i].reserve(100);
		}

		time_recorder.record_elapsed_miliseconds_and_reset("pre-allocate colliding IDs");

		#pragma omp parallel
		{
			std::vector<UnsignedInt> colliding_ids;
			colliding_ids.reserve(100);

			std::vector< std::pair<Float, UnsignedInt> > distances_of_colliding_ids;
			distances_of_colliding_ids.reserve(100);

			#pragma omp for
			for(UnsignedInt i=0;i<N;i++)
			{
				spheres_searcher.find_colliding_ids(i, colliding_ids, true, result.all_exclusion_statuses[i]);
				if(!colliding_ids.empty())
				{
					distances_of_colliding_ids.resize(colliding_ids.size());
					for(std::size_t j=0;j<colliding_ids.size();j++)
					{
						distances_of_colliding_ids[j].first=distance_to_center_of_intersection_circle_of_two_spheres(spheres[i], spheres[colliding_ids[j]]);
						distances_of_colliding_ids[j].second=colliding_ids[j];
					}
					std::sort(distances_of_colliding_ids.begin(), distances_of_colliding_ids.end());
					for(std::size_t j=0;j<colliding_ids.size();j++)
					{
						colliding_ids[j]=distances_of_colliding_ids[j].second;
					}
					result.all_colliding_ids[i]=colliding_ids;
				}
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("detect all collisions");

		for(UnsignedInt i=0;i<result.all_colliding_ids.size();i++)
		{
			result.total_collisions+=result.all_colliding_ids[i].size();
		}

		result.total_collisions=result.total_collisions/2;

		time_recorder.record_elapsed_miliseconds_and_reset("count all collisions");

		result.relevant_collision_ids.reserve(result.total_collisions);
		for(UnsignedInt id_a=0;id_a<N;id_a++)
		{
			if(result.all_exclusion_statuses[id_a]==0)
			{
				for(UnsignedInt j=0;j<result.all_colliding_ids[id_a].size();j++)
				{
					const UnsignedInt id_b=result.all_colliding_ids[id_a][j];
					if((result.all_exclusion_statuses[id_b]==0 && result.all_colliding_ids[id_a].size()<result.all_colliding_ids[id_b].size()) || (id_a<id_b && result.all_colliding_ids[id_a].size()==result.all_colliding_ids[id_b].size()))
					{
						if(grouping_of_spheres.empty() || id_a>=grouping_of_spheres.size() || id_b>=grouping_of_spheres.size() || grouping_of_spheres[id_a]!=grouping_of_spheres[id_b])
						{
							if(result.periodic_links_of_spheres.empty() || id_a>=result.periodic_links_of_spheres.size() || id_b>=result.periodic_links_of_spheres.size() || result.periodic_links_of_spheres[id_a]==id_a || result.periodic_links_of_spheres[id_b]==id_b)
							{
								result.relevant_collision_ids.push_back(std::pair<UnsignedInt, UnsignedInt>(id_a, id_b));
							}
						}
					}
				}
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("collect relevant collision indices");
	}
};

}

#endif /* VORONOTALT_PREPARATION_FOR_TESSELLATION_H_ */
