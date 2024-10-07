#ifndef VORONOTALT_SIMPLIFIED_AW_TESSELLATION_H_
#define VORONOTALT_SIMPLIFIED_AW_TESSELLATION_H_

#include <map>

#include "spheres_container.h"
#include "simplified_aw_tessellation_contact_construction.h"
#include "time_recorder.h"

namespace voronotalt
{

class SimplifiedAWTessellation
{
public:
	struct ContactDescriptorSummary
	{
		Float area;
		Float arc_length;
		Float distance;
		UnsignedInt id_a;
		UnsignedInt id_b;

		ContactDescriptorSummary() noexcept :
			area(FLOATCONST(0.0)),
			arc_length(FLOATCONST(0.0)),
			distance(FLOATCONST(0.0)),
			id_a(0),
			id_b(0)
		{
		}

		void set(const SimplifiedAWTessellationContactConstruction::ContactDescriptor& cd) noexcept
		{
			if(cd.area>FLOATCONST(0.0))
			{
				id_a=cd.id_a;
				id_b=cd.id_b;
				area=cd.area;
				arc_length=cd.arc_length;
				distance=cd.distance;
			}
		}

		void ensure_ids_ordered() noexcept
		{
			if(id_a>id_b)
			{
				std::swap(id_a, id_b);
			}
		}
	};

	struct TotalContactDescriptorsSummary
	{
		Float area;
		Float arc_length;
		Float distance;
		UnsignedInt count;

		TotalContactDescriptorsSummary() noexcept :
			area(FLOATCONST(0.0)),
			arc_length(FLOATCONST(0.0)),
			distance(FLOATCONST(-1.0)),
			count(0)
		{
		}

		void add(const ContactDescriptorSummary& cds) noexcept
		{
			if(cds.area>FLOATCONST(0.0))
			{
				count++;
				area+=cds.area;
				arc_length+=cds.arc_length;
				distance=((distance<FLOATCONST(0.0)) ? cds.distance : std::min(distance, cds.distance));
			}
		}
	};

	struct Result
	{
		UnsignedInt total_spheres;
		UnsignedInt total_collisions;
		UnsignedInt total_relevant_collisions;
		std::vector<ContactDescriptorSummary> contacts_summaries;
		TotalContactDescriptorsSummary total_contacts_summary;

		Result() noexcept : total_spheres(0), total_collisions(0), total_relevant_collisions(0)
		{
		}
	};

	struct ResultGraphics
	{
		std::vector<SimplifiedAWTessellationContactConstruction::ContactDescriptorGraphics> contacts_graphics;
	};

	struct GroupedResult
	{
		std::vector<UnsignedInt> grouped_contacts_representative_ids;
		std::vector<TotalContactDescriptorsSummary> grouped_contacts_summaries;
	};

	static void construct_full_tessellation(
			const std::vector<SimpleSphere>& spheres,
			Result& result) noexcept
	{
		TimeRecorder time_recorder;
		ResultGraphics result_graphics;
		construct_full_tessellation(spheres, std::vector<int>(), false, result, result_graphics, time_recorder);
	}

	static void construct_full_tessellation(
			const std::vector<SimpleSphere>& spheres,
			const std::vector<int>& grouping_of_spheres,
			const bool with_graphics,
			Result& result,
			ResultGraphics& result_graphics,
			TimeRecorder& time_recorder) noexcept
	{
		SpheresContainer spheres_container;
		spheres_container.init(spheres, time_recorder);

		SpheresContainer::ResultOfPreparationForTessellation preparation_result;
		spheres_container.prepare_for_tessellation(grouping_of_spheres, preparation_result, time_recorder);

		time_recorder.reset();

		result=Result();
		result_graphics=ResultGraphics();

		result.total_spheres=spheres_container.input_spheres().size();
		result.total_collisions=spheres_container.total_collisions();
		result.total_relevant_collisions=preparation_result.relevant_collision_ids.size();

		std::vector<ContactDescriptorSummary> possible_contacts_summaries(preparation_result.relevant_collision_ids.size());

		std::vector<SimplifiedAWTessellationContactConstruction::ContactDescriptorGraphics> possible_contacts_graphics;
		if(with_graphics)
		{
			possible_contacts_graphics.resize(possible_contacts_summaries.size());
		}

		time_recorder.record_elapsed_miliseconds_and_reset("allocate possible contact summaries");

#ifdef VORONOTALT_OPENMP
#pragma omp parallel
#endif
		{
			SimplifiedAWTessellationContactConstruction::ContactDescriptor cd;
			cd.neighbor_descriptors.reserve(24);

#ifdef VORONOTALT_OPENMP
#pragma omp for
#endif
			for(UnsignedInt i=0;i<preparation_result.relevant_collision_ids.size();i++)
			{
				const UnsignedInt id_a=preparation_result.relevant_collision_ids[i].first;
				const UnsignedInt id_b=preparation_result.relevant_collision_ids[i].second;
				if(SimplifiedAWTessellationContactConstruction::construct_contact_descriptor(spheres_container.populated_spheres(), spheres_container.all_exclusion_statuses(), id_a, id_b, spheres_container.all_colliding_ids()[id_a], 0.2, 5, with_graphics, cd))
				{
					possible_contacts_summaries[i].set(cd);
					if(with_graphics)
					{
						possible_contacts_graphics[i]=cd.graphics;
					}
				}
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("construct contacts");

		UnsignedInt number_of_valid_contact_summaries=0;
		for(UnsignedInt i=0;i<possible_contacts_summaries.size();i++)
		{
			if(possible_contacts_summaries[i].area>FLOATCONST(0.0))
			{
				number_of_valid_contact_summaries++;
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("count valid contact summaries");

		std::vector<UnsignedInt> ids_of_valid_pairs;
		ids_of_valid_pairs.reserve(number_of_valid_contact_summaries);

		for(UnsignedInt i=0;i<possible_contacts_summaries.size();i++)
		{
			if(possible_contacts_summaries[i].area>FLOATCONST(0.0))
			{
				ids_of_valid_pairs.push_back(i);
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("collect indices of valid contact summaries");

		result.contacts_summaries.resize(ids_of_valid_pairs.size());

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<ids_of_valid_pairs.size();i++)
			{
				result.contacts_summaries[i]=possible_contacts_summaries[ids_of_valid_pairs[i]];
				result.contacts_summaries[i].ensure_ids_ordered();
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("copy valid contact summaries");

		for(UnsignedInt i=0;i<result.contacts_summaries.size();i++)
		{
			result.total_contacts_summary.add(result.contacts_summaries[i]);
		}

		time_recorder.record_elapsed_miliseconds_and_reset("accumulate total contacts summary");

		if(with_graphics)
		{
			result_graphics.contacts_graphics.resize(ids_of_valid_pairs.size());

			for(UnsignedInt i=0;i<ids_of_valid_pairs.size();i++)
			{
				result_graphics.contacts_graphics[i]=possible_contacts_graphics[ids_of_valid_pairs[i]];
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("copy valid contacts graphics");
	}

	static bool group_results(
			const Result& full_result,
			const std::vector<int>& grouping_of_spheres,
			GroupedResult& grouped_result) noexcept
	{
		TimeRecorder time_recorder;
		return group_results(full_result, grouping_of_spheres, grouped_result, time_recorder);
	}

	static bool group_results(
			const Result& full_result,
			const std::vector<int>& grouping_of_spheres,
			GroupedResult& grouped_result,
			TimeRecorder& time_recorder) noexcept
	{
		time_recorder.reset();

		grouped_result=GroupedResult();

		if(!grouping_of_spheres.empty() && grouping_of_spheres.size()==full_result.total_spheres)
		{
			grouped_result.grouped_contacts_representative_ids.reserve(full_result.contacts_summaries.size()/5);
			grouped_result.grouped_contacts_summaries.reserve(full_result.contacts_summaries.size()/5);

			std::map< std::pair<int, int>, UnsignedInt > map_of_groups;

			for(UnsignedInt i=0;i<full_result.contacts_summaries.size();i++)
			{
				const ContactDescriptorSummary& cds=full_result.contacts_summaries[i];
				if(cds.id_a<grouping_of_spheres.size() && cds.id_b<grouping_of_spheres.size())
				{
					std::pair<int, int> group_id(grouping_of_spheres[cds.id_a], grouping_of_spheres[cds.id_b]);
					if(group_id.first!=group_id.second)
					{
						if(group_id.first>group_id.second)
						{
							std::swap(group_id.first, group_id.second);
						}
						UnsignedInt group_index=0;
						std::map< std::pair<int, int>, UnsignedInt >::const_iterator it=map_of_groups.find(group_id);
						if(it==map_of_groups.end())
						{
							group_index=grouped_result.grouped_contacts_summaries.size();
							grouped_result.grouped_contacts_representative_ids.push_back(i);
							grouped_result.grouped_contacts_summaries.push_back(TotalContactDescriptorsSummary());
							map_of_groups[group_id]=group_index;
						}
						else
						{
							group_index=it->second;
						}
						grouped_result.grouped_contacts_summaries[group_index].add(cds);
					}
				}
			}

			time_recorder.record_elapsed_miliseconds_and_reset("grouped contacts summaries");
		}

		return (!grouped_result.grouped_contacts_summaries.empty());
	}
};

}

#endif /* VORONOTALT_SIMPLIFIED_AW_TESSELLATION_H_ */
