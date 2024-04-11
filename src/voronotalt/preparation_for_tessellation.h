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
		std::vector<int> all_exclusion_statuses;
		std::vector< std::vector<UnsignedInt> > all_colliding_ids;
		std::vector< std::pair<UnsignedInt, UnsignedInt> > relevant_collision_ids;
		UnsignedInt total_spheres;
		UnsignedInt total_collisions;

		Result() : total_spheres(0), total_collisions(0)
		{
		}
	};

	static void prepare_for_tessellation(
			const std::vector<SimpleSphere>& spheres,
			const std::vector<int>& grouping_of_spheres,
			Result& result,
			TimeRecorder& time_recorder)
	{
		time_recorder.reset();

		result=Result();

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
							result.relevant_collision_ids.push_back(std::pair<UnsignedInt, UnsignedInt>(id_a, id_b));
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
