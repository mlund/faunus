#ifndef VORONOTALT_RADICAL_TESSELLATION_H_
#define VORONOTALT_RADICAL_TESSELLATION_H_

#include <map>

#include "preparation_for_tessellation.h"
#include "radical_tessellation_contact_construction.h"
#include "time_recorder.h"

namespace voronotalt
{

class RadicalTessellation
{
public:
	struct ContactDescriptorSummary
	{
		Float area;
		Float arc_length;
		Float solid_angle_a;
		Float solid_angle_b;
		Float pyramid_volume_a;
		Float pyramid_volume_b;
		Float distance;
		UnsignedInt flags;
		UnsignedInt id_a;
		UnsignedInt id_b;

		ContactDescriptorSummary() :
			area(FLOATCONST(0.0)),
			arc_length(FLOATCONST(0.0)),
			solid_angle_a(FLOATCONST(0.0)),
			solid_angle_b(FLOATCONST(0.0)),
			pyramid_volume_a(FLOATCONST(0.0)),
			pyramid_volume_b(FLOATCONST(0.0)),
			distance(FLOATCONST(0.0)),
			flags(0),
			id_a(0),
			id_b(0)
		{
		}

		void set(const RadicalTessellationContactConstruction::ContactDescriptor& cd)
		{
			if(cd.area>FLOATCONST(0.0))
			{
				id_a=cd.id_a;
				id_b=cd.id_b;
				area=cd.area;
				arc_length=(cd.sum_of_arc_angles*cd.intersection_circle_sphere.r);
				solid_angle_a=cd.solid_angle_a;
				solid_angle_b=cd.solid_angle_b;
				pyramid_volume_a=cd.pyramid_volume_a;
				pyramid_volume_b=cd.pyramid_volume_b;
				distance=cd.distance;
				flags=cd.flags;
			}
		}

		void ensure_ids_ordered()
		{
			if(id_a>id_b)
			{
				std::swap(id_a, id_b);
				std::swap(solid_angle_a, solid_angle_b);
				std::swap(pyramid_volume_a, pyramid_volume_b);
			}
		}
	};

	struct CellContactDescriptorsSummary
	{
		Float area;
		Float arc_length;
		Float explained_solid_angle_positive;
		Float explained_solid_angle_negative;
		Float explained_pyramid_volume_positive;
		Float explained_pyramid_volume_negative;
		Float sas_area;
		Float sas_inside_volume;
		UnsignedInt id;
		UnsignedInt count;
		int stage;

		CellContactDescriptorsSummary() :
			area(FLOATCONST(0.0)),
			arc_length(FLOATCONST(0.0)),
			explained_solid_angle_positive(FLOATCONST(0.0)),
			explained_solid_angle_negative(FLOATCONST(0.0)),
			explained_pyramid_volume_positive(FLOATCONST(0.0)),
			explained_pyramid_volume_negative(FLOATCONST(0.0)),
			sas_area(FLOATCONST(0.0)),
			sas_inside_volume(FLOATCONST(0.0)),
			id(0),
			count(0),
			stage(0)
		{
		}

		void add(const ContactDescriptorSummary& cds)
		{
			if(cds.area>FLOATCONST(0.0) && (cds.id_a==id || cds.id_b==id))
			{
				count++;
				area+=cds.area;
				arc_length+=cds.arc_length;
				explained_solid_angle_positive+=std::max(FLOATCONST(0.0), (cds.id_a==id ? cds.solid_angle_a : cds.solid_angle_b));
				explained_solid_angle_negative+=FLOATCONST(0.0)-std::min(FLOATCONST(0.0), (cds.id_a==id ? cds.solid_angle_a : cds.solid_angle_b));
				explained_pyramid_volume_positive+=std::max(FLOATCONST(0.0), (cds.id_a==id ? cds.pyramid_volume_a : cds.pyramid_volume_b));
				explained_pyramid_volume_negative+=FLOATCONST(0.0)-std::min(FLOATCONST(0.0), (cds.id_a==id ? cds.pyramid_volume_a : cds.pyramid_volume_b));
				stage=1;
			}
		}

		void add(const UnsignedInt new_id, const ContactDescriptorSummary& cds)
		{
			if(cds.area>FLOATCONST(0.0))
			{
				if(stage==0)
				{
					id=new_id;
				}
				add(cds);
			}
		}

		void compute_sas(const Float r)
		{
			if(stage==1)
			{
				sas_area=FLOATCONST(0.0);
				sas_inside_volume=FLOATCONST(0.0);
				if(arc_length>FLOATCONST(0.0) && !equal(explained_solid_angle_positive, explained_solid_angle_negative))
				{
					if(explained_solid_angle_positive>explained_solid_angle_negative)
					{
						sas_area=((FLOATCONST(4.0)*PIVALUE)-std::max(FLOATCONST(0.0), explained_solid_angle_positive-explained_solid_angle_negative))*(r*r);
					}
					else if(explained_solid_angle_negative>explained_solid_angle_positive)
					{
						sas_area=(std::max(FLOATCONST(0.0), explained_solid_angle_negative-explained_solid_angle_positive))*(r*r);
					}
					sas_inside_volume=(sas_area*r/FLOATCONST(3.0))+explained_pyramid_volume_positive-explained_pyramid_volume_negative;
					if(sas_inside_volume>(FLOATCONST(4.0)/FLOATCONST(3.0)*PIVALUE*r*r*r))
					{
						sas_area=FLOATCONST(0.0);
						sas_inside_volume=explained_pyramid_volume_positive-explained_pyramid_volume_negative;
					}
				}
				else
				{
					sas_inside_volume=explained_pyramid_volume_positive-explained_pyramid_volume_negative;
				}
				stage=2;
			}
		}

		void compute_sas_detached(const UnsignedInt new_id, const Float r)
		{
			if(stage==0)
			{
				id=new_id;
				sas_area=(FLOATCONST(4.0)*PIVALUE)*(r*r);
				sas_inside_volume=(sas_area*r/FLOATCONST(3.0));
				stage=2;
			}
		}
	};

	struct TotalContactDescriptorsSummary
	{
		Float area;
		Float arc_length;
		Float distance;
		UnsignedInt count;

		TotalContactDescriptorsSummary() :
			area(FLOATCONST(0.0)),
			arc_length(FLOATCONST(0.0)),
			distance(FLOATCONST(-1.0)),
			count(0)
		{
		}

		void add(const ContactDescriptorSummary& cds)
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

	struct TotalCellContactDescriptorsSummary
	{
		Float sas_area;
		Float sas_inside_volume;
		UnsignedInt count;

		TotalCellContactDescriptorsSummary() :
			sas_area(FLOATCONST(0.0)),
			sas_inside_volume(FLOATCONST(0.0)),
			count(0)
		{
		}

		void add(const CellContactDescriptorsSummary& ccds)
		{
			if(ccds.stage==2)
			{
				count++;
				sas_area+=ccds.sas_area;
				sas_inside_volume+=ccds.sas_inside_volume;
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
		std::vector<CellContactDescriptorsSummary> cells_summaries;
		TotalCellContactDescriptorsSummary total_cells_summary;

		Result() : total_spheres(0), total_collisions(0), total_relevant_collisions(0)
		{
		}
	};

	struct ResultGraphics
	{
		std::vector<RadicalTessellationContactConstruction::ContactDescriptorGraphics> contacts_graphics;
	};

	struct GroupedResult
	{
		std::vector<UnsignedInt> grouped_contacts_representative_ids;
		std::vector<TotalContactDescriptorsSummary> grouped_contacts_summaries;
		std::vector<UnsignedInt> grouped_cells_representative_ids;
		std::vector<TotalCellContactDescriptorsSummary> grouped_cells_summaries;
	};

	static void construct_full_tessellation(
			const std::vector<SimpleSphere>& spheres,
			Result& result)
	{
		TimeRecorder time_recorder;
		ResultGraphics result_graphics;
		construct_full_tessellation(spheres, std::vector<int>(), false, true, result, result_graphics, time_recorder);
	}

	static void construct_full_tessellation(
			const std::vector<SimpleSphere>& spheres,
			const std::vector<int>& grouping_of_spheres,
			const bool with_graphics,
			const bool summarize_cells,
			Result& result,
			ResultGraphics& result_graphics,
			TimeRecorder& time_recorder)
	{
		PreparationForTessellation::Result preparation_result;
		PreparationForTessellation::prepare_for_tessellation(spheres, grouping_of_spheres, preparation_result, time_recorder);
		construct_full_tessellation(spheres, preparation_result, with_graphics, summarize_cells, result, result_graphics, time_recorder);
	}

	static void construct_full_tessellation(
			const std::vector<SimpleSphere>& spheres,
			const PreparationForTessellation::Result& preparation_result,
			const bool with_graphics,
			const bool summarize_cells,
			Result& result,
			ResultGraphics& result_graphics,
			TimeRecorder& time_recorder)
	{
		time_recorder.reset();

		result=Result();
		result_graphics=ResultGraphics();

		result.total_spheres=preparation_result.total_spheres;
		result.total_collisions=preparation_result.total_collisions;
		result.total_relevant_collisions=preparation_result.relevant_collision_ids.size();

		std::vector<ContactDescriptorSummary> possible_contacts_summaries(preparation_result.relevant_collision_ids.size());

		std::vector<RadicalTessellationContactConstruction::ContactDescriptorGraphics> possible_contacts_graphics;
		if(with_graphics)
		{
			possible_contacts_graphics.resize(possible_contacts_summaries.size());
		}

		time_recorder.record_elapsed_miliseconds_and_reset("allocate possible contact summaries");

		#pragma omp parallel
		{
			RadicalTessellationContactConstruction::ContactDescriptor cd;
			cd.contour.reserve(12);

			#pragma omp for
			for(UnsignedInt i=0;i<preparation_result.relevant_collision_ids.size();i++)
			{
				const UnsignedInt id_a=preparation_result.relevant_collision_ids[i].first;
				const UnsignedInt id_b=preparation_result.relevant_collision_ids[i].second;
				if(RadicalTessellationContactConstruction::construct_contact_descriptor(spheres, preparation_result.all_exclusion_statuses, id_a, id_b, preparation_result.all_colliding_ids[id_a], cd))
				{
					possible_contacts_summaries[i].set(cd);
					if(with_graphics)
					{
						RadicalTessellationContactConstruction::construct_contact_descriptor_graphics(cd, 0.2, possible_contacts_graphics[i]);
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

		#pragma omp parallel
		{
			#pragma omp for
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

		if(summarize_cells)
		{
			result.cells_summaries.resize(preparation_result.total_spheres);

			for(UnsignedInt i=0;i<result.contacts_summaries.size();i++)
			{
				const ContactDescriptorSummary& cds=result.contacts_summaries[i];
				result.cells_summaries[cds.id_a].add(cds.id_a, cds);
				result.cells_summaries[cds.id_b].add(cds.id_b, cds);
			}

			time_recorder.record_elapsed_miliseconds_and_reset("accumulate cell summaries");

			for(UnsignedInt i=0;i<result.cells_summaries.size();i++)
			{
				result.cells_summaries[i].compute_sas(spheres[i].r);
				if(result.cells_summaries[i].stage==0 && preparation_result.all_exclusion_statuses[i]==0 && preparation_result.all_colliding_ids[i].empty())
				{
					result.cells_summaries[i].compute_sas_detached(i, spheres[i].r);
				}
			}

			time_recorder.record_elapsed_miliseconds_and_reset("compute sas for cell summaries");

			for(UnsignedInt i=0;i<result.cells_summaries.size();i++)
			{
				result.total_cells_summary.add(result.cells_summaries[i]);
			}

			time_recorder.record_elapsed_miliseconds_and_reset("accumulate total cells summary");
		}
	}

	static bool group_results(
			const Result& full_result,
			const std::vector<int>& grouping_of_spheres,
			GroupedResult& grouped_result)
	{
		TimeRecorder time_recorder;
		return group_results(full_result, grouping_of_spheres, grouped_result, time_recorder);
	}

	static bool group_results(
			const Result& full_result,
			const std::vector<int>& grouping_of_spheres,
			GroupedResult& grouped_result,
			TimeRecorder& time_recorder)
	{
		time_recorder.reset();

		grouped_result=GroupedResult();

		if(!grouping_of_spheres.empty() && grouping_of_spheres.size()==full_result.total_spheres)
		{
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

			if(!full_result.cells_summaries.empty() && full_result.cells_summaries.size()==grouping_of_spheres.size())
			{
				grouped_result.grouped_cells_representative_ids.reserve(full_result.cells_summaries.size()/5);
				grouped_result.grouped_cells_summaries.reserve(full_result.cells_summaries.size()/5);

				std::map< int, UnsignedInt > map_of_groups;

				for(UnsignedInt i=0;i<full_result.cells_summaries.size();i++)
				{
					const CellContactDescriptorsSummary& ccds=full_result.cells_summaries[i];
					if(ccds.stage==2 && ccds.id<grouping_of_spheres.size())
					{
						const int group_id=grouping_of_spheres[ccds.id];
						UnsignedInt group_index=0;
						std::map< int, UnsignedInt >::const_iterator it=map_of_groups.find(group_id);
						if(it==map_of_groups.end())
						{
							group_index=grouped_result.grouped_cells_summaries.size();
							grouped_result.grouped_cells_representative_ids.push_back(i);
							grouped_result.grouped_cells_summaries.push_back(TotalCellContactDescriptorsSummary());
							map_of_groups[group_id]=group_index;
						}
						else
						{
							group_index=it->second;
						}
						grouped_result.grouped_cells_summaries[group_index].add(ccds);
					}
				}

				time_recorder.record_elapsed_miliseconds_and_reset("grouped cells summaries");
			}
		}

		return (!grouped_result.grouped_contacts_summaries.empty());
	}
};

}

#endif /* VORONOTALT_RADICAL_TESSELLATION_H_ */
