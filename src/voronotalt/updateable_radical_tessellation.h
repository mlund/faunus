#ifndef VORONOTALT_UPDATEABLE_RADICAL_TESSELLATION_H_
#define VORONOTALT_UPDATEABLE_RADICAL_TESSELLATION_H_

#include "radical_tessellation.h"

namespace voronotalt
{

class UpdateableRadicalTessellation
{
public:
	struct Result
	{
		std::vector<RadicalTessellation::CellContactDescriptorsSummary> cells_summaries;
		std::vector< std::vector<RadicalTessellation::ContactDescriptorSummary> > contacts_summaries;
		std::vector< std::vector<RadicalTessellation::ContactDescriptorSummary> > contacts_summaries_with_redundancy_in_periodic_box;

		Result() noexcept
		{
		}

		bool empty() noexcept
		{
			return (cells_summaries.empty() || contacts_summaries.empty());
		}
	};

	struct ResultSummary
	{
		RadicalTessellation::TotalContactDescriptorsSummary total_contacts_summary;
		RadicalTessellation::TotalCellContactDescriptorsSummary total_cells_summary;

		ResultSummary() noexcept
		{
		}
	};

	UpdateableRadicalTessellation() noexcept : backup_enabled_(false), in_sync_with_backup_(false)
	{
	}

	explicit UpdateableRadicalTessellation(const bool backup_enabled) noexcept : backup_enabled_(backup_enabled), in_sync_with_backup_(false)
	{
	}

	bool init(const std::vector<SimpleSphere>& input_spheres) noexcept
	{
		TimeRecorder time_recorder;
		return init(input_spheres, time_recorder);
	}

	bool init(const std::vector<SimpleSphere>& input_spheres, TimeRecorder& time_recorder) noexcept
	{
		return init(input_spheres, PeriodicBox(), time_recorder);
	}

	bool init(const std::vector<SimpleSphere>& input_spheres, const PeriodicBox& periodic_box) noexcept
	{
		TimeRecorder time_recorder;
		return init(input_spheres, periodic_box, time_recorder);
	}

	bool init(const std::vector<SimpleSphere>& input_spheres, const PeriodicBox& periodic_box, TimeRecorder& time_recorder) noexcept
	{
		prepare_for_possible_init_or_update(time_recorder);

		state_.spheres_container.init(input_spheres, periodic_box, time_recorder);

		in_sync_with_backup_=false;

		RadicalTessellation::ResultGraphics result_graphics;
		RadicalTessellation::construct_full_tessellation(state_.spheres_container, std::vector<int>(), std::vector<int>(), false, true, FLOATCONST(0.0), std::vector<Float>(), buffered_temporary_storage_.tessellation_result, result_graphics, time_recorder);

		involvement_of_spheres_for_update_.clear();

		init_result_from_tessellation_result();

		return !state_.result.empty();
	}

	bool update(const std::vector<SimpleSphere>& new_input_spheres) noexcept
	{
		TimeRecorder time_recorder;
		return update(new_input_spheres, std::vector<UnsignedInt>(), false, time_recorder);
	}

	bool update(const std::vector<SimpleSphere>& new_input_spheres, TimeRecorder& time_recorder) noexcept
	{
		return update(new_input_spheres, std::vector<UnsignedInt>(), false, time_recorder);
	}

	bool update(const std::vector<SimpleSphere>& new_input_spheres, const std::vector<UnsignedInt>& ids_of_changed_input_spheres) noexcept
	{
		TimeRecorder time_recorder;
		return update(new_input_spheres, ids_of_changed_input_spheres, true, time_recorder);
	}

	bool update(const std::vector<SimpleSphere>& new_input_spheres, const std::vector<UnsignedInt>& ids_of_changed_input_spheres, TimeRecorder& time_recorder) noexcept
	{
		return update(new_input_spheres, ids_of_changed_input_spheres, true, time_recorder);
	}

	bool update(const std::vector<SimpleSphere>& new_input_spheres, const std::vector<UnsignedInt>& provided_ids_of_changed_input_spheres, const bool trust_provided_ids_of_changed_input_spheres, TimeRecorder& time_recorder) noexcept
	{
		prepare_for_possible_init_or_update(time_recorder);

		if(trust_provided_ids_of_changed_input_spheres && provided_ids_of_changed_input_spheres.empty())
		{
			return false;
		}

		if(!state_.spheres_container.update(new_input_spheres, provided_ids_of_changed_input_spheres, trust_provided_ids_of_changed_input_spheres, state_.ids_of_changed_input_spheres, state_.ids_of_affected_input_spheres, time_recorder))
		{
			return false;
		}

		in_sync_with_backup_=false;

		if(state_.ids_of_affected_input_spheres.empty())
		{
			RadicalTessellation::ResultGraphics result_graphics;
			RadicalTessellation::construct_full_tessellation(state_.spheres_container, std::vector<int>(), std::vector<int>(), false, true, FLOATCONST(0.0), std::vector<Float>(), buffered_temporary_storage_.tessellation_result, result_graphics, time_recorder);
			init_result_from_tessellation_result();
		}
		else
		{
			update_using_current_state(time_recorder);
		}

		return true;
	}

	bool update_by_setting_exclusion_mask(const UnsignedInt id_of_targeted_input_sphere, const bool new_exclusion_status) noexcept
	{
		TimeRecorder time_recorder;
		return update_by_setting_exclusion_mask(id_of_targeted_input_sphere, new_exclusion_status, time_recorder);
	}

	bool update_by_setting_exclusion_mask(const UnsignedInt id_of_targeted_input_sphere, const bool new_exclusion_status, TimeRecorder& time_recorder) noexcept
	{
		if(state_.result.empty() || id_of_targeted_input_sphere>=state_.result.contacts_summaries.size()
				|| id_of_targeted_input_sphere>=state_.spheres_container.all_exclusion_statuses().size()
				|| (new_exclusion_status ? state_.spheres_container.all_exclusion_statuses()[id_of_targeted_input_sphere]>0 : state_.spheres_container.all_exclusion_statuses()[id_of_targeted_input_sphere]<1))
		{
			return false;
		}

		prepare_for_possible_init_or_update(time_recorder);

		if(new_exclusion_status)
		{
			if(!state_.spheres_container.update_by_setting_exclusion_status(id_of_targeted_input_sphere, true))
			{
				return false;
			}

			state_.ids_of_changed_input_spheres.push_back(id_of_targeted_input_sphere);
			state_.ids_of_affected_input_spheres.push_back(id_of_targeted_input_sphere);

			for(UnsignedInt i=0;i<state_.result.contacts_summaries[id_of_targeted_input_sphere].size();i++)
			{
				const RadicalTessellation::ContactDescriptorSummary& cds=state_.result.contacts_summaries[id_of_targeted_input_sphere][i];
				if(cds.id_a!=id_of_targeted_input_sphere)
				{
					state_.ids_of_affected_input_spheres.push_back(cds.id_a);
				}
				else if(cds.id_b!=id_of_targeted_input_sphere)
				{
					state_.ids_of_affected_input_spheres.push_back(cds.id_b);
				}
			}
		}
		else
		{
			if(!state_.spheres_container.update_by_setting_exclusion_status(id_of_targeted_input_sphere, false, state_.ids_of_changed_input_spheres, state_.ids_of_affected_input_spheres))
			{
				return false;
			}
		}

		in_sync_with_backup_=false;

		update_using_current_state(time_recorder);

		return true;
	}

	bool restore_from_backup() noexcept
	{
		if(backup_enabled_ && !in_sync_with_backup_)
		{
			state_.assign_to_undo_update(state_backup_);
			in_sync_with_backup_=true;
		}

		return in_sync_with_backup_;
	}

	bool calculate_second_order_cell_volumes(std::vector< std::vector<Float> >& all_result_volumes_for_contacts_summaries) noexcept
	{
		all_result_volumes_for_contacts_summaries.clear();

		if(!backup_enabled_ || state_.result.empty())
		{
			return false;
		}

		const UnsignedInt N=state_.result.contacts_summaries.size();

		std::vector< std::vector<UnsignedInt> > all_neighbor_ids(N);
		std::vector< std::vector<Float> > all_volumes_after_masking(N);

		for(UnsignedInt i=0;i<N;i++)
		{
			std::vector<UnsignedInt>& i_neighbor_ids=all_neighbor_ids[i];

			{
				const std::vector<RadicalTessellation::ContactDescriptorSummary>& cdss=state_.result.contacts_summaries[i];
				i_neighbor_ids.resize(cdss.size());
				for(UnsignedInt j=0;j<cdss.size();j++)
				{
					i_neighbor_ids[j]=(cdss[j].id_a!=i ? cdss[j].id_a : cdss[j].id_b);
				}
			}

			std::vector<Float>& i_volumes_after_masking=all_volumes_after_masking[i];
			i_volumes_after_masking.resize(i_neighbor_ids.size(), FLOATCONST(0.0));

			if(update_by_setting_exclusion_mask(i, true))
			{
				for(UnsignedInt j=0;j<i_neighbor_ids.size();j++)
				{
					i_volumes_after_masking[j]=state_.result.cells_summaries[i_neighbor_ids[j]].sas_inside_volume;
				}
				restore_from_backup();
			}
		}

		all_result_volumes_for_contacts_summaries.resize(N);

		for(UnsignedInt i=0;i<N;i++)
		{
			all_result_volumes_for_contacts_summaries[i].resize(all_neighbor_ids[i].size());
			const Float i_volume_before_masking=state_.result.cells_summaries[i].sas_inside_volume;
			for(UnsignedInt j=0;j<all_neighbor_ids[i].size();j++)
			{
				const UnsignedInt neighbor_id=all_neighbor_ids[i][j];
				const Float neighbor_volume_before_masking=state_.result.cells_summaries[neighbor_id].sas_inside_volume;
				const Float neighbor_volume_after_masking_i=all_volumes_after_masking[i][j];
				Float i_volume_after_masking_neighbor=FLOATCONST(0.0);
				for(UnsignedInt relative_i=0;relative_i<all_volumes_after_masking[neighbor_id].size() && i_volume_after_masking_neighbor<=FLOATCONST(0.0);relative_i++)
				{
					if(all_neighbor_ids[neighbor_id][relative_i]==i)
					{
						i_volume_after_masking_neighbor=all_volumes_after_masking[neighbor_id][relative_i];
					}
				}
				all_result_volumes_for_contacts_summaries[i][j]=(i_volume_after_masking_neighbor+neighbor_volume_after_masking_i)-(i_volume_before_masking+neighbor_volume_before_masking);
			}
		}

		return true;
	}

	bool backup_enabled() const noexcept
	{
		return backup_enabled_;
	}

	bool in_sync_with_backup() const noexcept
	{
		return in_sync_with_backup_;
	}

	const std::vector<SimpleSphere>& input_spheres() const noexcept
	{
		return state_.spheres_container.input_spheres();
	}

	bool exclusion_status_of_input_sphere(const UnsignedInt id_of_input_sphere) const noexcept
	{
		return (id_of_input_sphere<state_.spheres_container.all_exclusion_statuses().size() && state_.spheres_container.all_exclusion_statuses()[id_of_input_sphere]>0);
	}

	const Result& result() const noexcept
	{
		return state_.result;
	}

	ResultSummary result_summary() const noexcept
	{
		ResultSummary result_summary;
		for(UnsignedInt i=0;i<state_.result.contacts_summaries.size();i++)
		{
			for(UnsignedInt j=0;j<state_.result.contacts_summaries[i].size();j++)
			{
				const RadicalTessellation::ContactDescriptorSummary& cds=state_.result.contacts_summaries[i][j];
				if(cds.id_a==i)
				{
					result_summary.total_contacts_summary.add(cds);
				}
			}
		}
		for(UnsignedInt i=0;i<state_.result.cells_summaries.size();i++)
		{
			result_summary.total_cells_summary.add(state_.result.cells_summaries[i]);
		}
		return result_summary;
	}

	const std::vector<UnsignedInt>& last_update_ids_of_changed_input_spheres() const noexcept
	{
		return state_.ids_of_changed_input_spheres;
	}

	const std::vector<UnsignedInt>& last_update_ids_of_affected_input_spheres() const noexcept
	{
		return state_.ids_of_affected_input_spheres;
	}

	bool last_update_was_full_reinit() const noexcept
	{
		return state_.last_update_was_full_reinit;
	}

private:
	struct BufferedTemporaryStorage
	{
		RadicalTessellation::Result tessellation_result;

		void clear() noexcept
		{
			tessellation_result.clear();
		}
	};

	class ConditionToRemoveContact
	{
	public:
		explicit ConditionToRemoveContact(const std::vector<int>& involvement) noexcept : involvement_(involvement)
		{
		}

		bool operator()(const RadicalTessellation::ContactDescriptorSummary& cds) noexcept
		{
			return (involvement_.empty() || (involvement_[cds.id_a%involvement_.size()]>0 && involvement_[cds.id_b%involvement_.size()]>0));
		}

	private:
		const std::vector<int>& involvement_;
	};

	struct State
	{
		SpheresContainer spheres_container;
		Result result;
		std::vector<UnsignedInt> ids_of_changed_input_spheres;
		std::vector<UnsignedInt> ids_of_affected_input_spheres;
		bool last_update_was_full_reinit;

		State() noexcept : last_update_was_full_reinit(true)
		{
		}

		void assign(const State& obj) noexcept
		{
			spheres_container.assign(obj.spheres_container);

			result.cells_summaries.resize(obj.result.cells_summaries.size());
			result.contacts_summaries.resize(obj.result.contacts_summaries.size());
			result.contacts_summaries_with_redundancy_in_periodic_box.resize(obj.result.contacts_summaries_with_redundancy_in_periodic_box.size());
			ids_of_changed_input_spheres.resize(obj.ids_of_changed_input_spheres.size());
			ids_of_affected_input_spheres.resize(obj.ids_of_affected_input_spheres.size());

			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<obj.result.cells_summaries.size();i++)
				{
					result.cells_summaries[i]=obj.result.cells_summaries[i];
				}
			}

			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<obj.result.contacts_summaries.size();i++)
				{
					result.contacts_summaries[i]=obj.result.contacts_summaries[i];
				}
			}

			if(!obj.result.contacts_summaries_with_redundancy_in_periodic_box.empty())
			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<obj.result.contacts_summaries_with_redundancy_in_periodic_box.size();i++)
				{
					result.contacts_summaries_with_redundancy_in_periodic_box[i]=obj.result.contacts_summaries_with_redundancy_in_periodic_box[i];
				}
			}

			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<obj.ids_of_changed_input_spheres.size();i++)
				{
					ids_of_changed_input_spheres[i]=obj.ids_of_changed_input_spheres[i];
				}
			}

			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<obj.ids_of_affected_input_spheres.size();i++)
				{
					ids_of_affected_input_spheres[i]=obj.ids_of_affected_input_spheres[i];
				}
			}

			last_update_was_full_reinit=obj.last_update_was_full_reinit;
		}

		void assign(const State& obj, const bool assign_everything, const std::vector<UnsignedInt>& subset_of_ids_of_spheres) noexcept
		{
			if(!assign_everything && subset_of_ids_of_spheres.empty())
			{
				return;
			}

			if(assign_everything
					|| result.cells_summaries.size()!=obj.result.cells_summaries.size()
					|| result.contacts_summaries.size()!=obj.result.contacts_summaries.size()
					|| result.contacts_summaries_with_redundancy_in_periodic_box.size()!=obj.result.contacts_summaries_with_redundancy_in_periodic_box.size())
			{
				assign(obj);
				return;
			}

			const bool periodic=!obj.result.contacts_summaries_with_redundancy_in_periodic_box.empty();

			for(UnsignedInt i=0;i<subset_of_ids_of_spheres.size();i++)
			{
				const UnsignedInt sphere_id=subset_of_ids_of_spheres[i];
				if(sphere_id>=result.cells_summaries.size() || (periodic && sphere_id>=result.contacts_summaries_with_redundancy_in_periodic_box.size()))
				{
					assign(obj);
					return;
				}
			}

			spheres_container.assign(obj.spheres_container, ids_of_changed_input_spheres);

			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<subset_of_ids_of_spheres.size();i++)
				{
					const UnsignedInt sphere_id=subset_of_ids_of_spheres[i];
					result.cells_summaries[sphere_id]=obj.result.cells_summaries[sphere_id];
					result.contacts_summaries[sphere_id]=obj.result.contacts_summaries[sphere_id];
					if(periodic)
					{
						result.contacts_summaries_with_redundancy_in_periodic_box[sphere_id]=obj.result.contacts_summaries_with_redundancy_in_periodic_box[sphere_id];
					}
				}
			}

			ids_of_changed_input_spheres.resize(obj.ids_of_changed_input_spheres.size());
			ids_of_affected_input_spheres.resize(obj.ids_of_affected_input_spheres.size());

			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<obj.ids_of_changed_input_spheres.size();i++)
				{
					ids_of_changed_input_spheres[i]=obj.ids_of_changed_input_spheres[i];
				}
			}

			{
				#pragma omp parallel for
				for(UnsignedInt i=0;i<obj.ids_of_affected_input_spheres.size();i++)
				{
					ids_of_affected_input_spheres[i]=obj.ids_of_affected_input_spheres[i];
				}
			}

			last_update_was_full_reinit=obj.last_update_was_full_reinit;
		}

		void assign_to_undo_update(const State& obj) noexcept
		{
			assign(obj, last_update_was_full_reinit, ids_of_affected_input_spheres);
		}

		void assign_to_apply_update(const State& obj) noexcept
		{
			assign(obj, obj.last_update_was_full_reinit, obj.ids_of_affected_input_spheres);
		}
	};

	void init_result_from_tessellation_result() noexcept
	{
		state_.ids_of_changed_input_spheres.clear();
		state_.ids_of_affected_input_spheres.clear();
		state_.last_update_was_full_reinit=true;

		state_.result.cells_summaries.swap(buffered_temporary_storage_.tessellation_result.cells_summaries);

		{
			state_.result.contacts_summaries.resize(state_.spheres_container.input_spheres().size());
			for(UnsignedInt i=0;i<state_.result.contacts_summaries.size();i++)
			{
				state_.result.contacts_summaries[i].clear();
			}
			for(UnsignedInt i=0;i<buffered_temporary_storage_.tessellation_result.contacts_summaries.size();i++)
			{
				const RadicalTessellation::ContactDescriptorSummary& cds=buffered_temporary_storage_.tessellation_result.contacts_summaries[i];
				state_.result.contacts_summaries[cds.id_a].push_back(cds);
				state_.result.contacts_summaries[cds.id_b].push_back(cds);
			}
		}

		if(state_.spheres_container.periodic_box().enabled())
		{
			state_.result.contacts_summaries_with_redundancy_in_periodic_box.resize(state_.spheres_container.input_spheres().size());
			for(UnsignedInt i=0;i<state_.result.contacts_summaries_with_redundancy_in_periodic_box.size();i++)
			{
				state_.result.contacts_summaries_with_redundancy_in_periodic_box[i].clear();
			}
			const std::vector<RadicalTessellation::ContactDescriptorSummary>& all_contacts_summaries=(buffered_temporary_storage_.tessellation_result.contacts_summaries_with_redundancy_in_periodic_box.empty() ? buffered_temporary_storage_.tessellation_result.contacts_summaries : buffered_temporary_storage_.tessellation_result.contacts_summaries_with_redundancy_in_periodic_box);
			for(UnsignedInt i=0;i<all_contacts_summaries.size();i++)
			{
				const RadicalTessellation::ContactDescriptorSummary& cds=all_contacts_summaries[i];
				if(cds.id_a<state_.result.contacts_summaries_with_redundancy_in_periodic_box.size())
				{
					state_.result.contacts_summaries_with_redundancy_in_periodic_box[cds.id_a].push_back(cds);
				}
				if(cds.id_b<state_.result.contacts_summaries_with_redundancy_in_periodic_box.size())
				{
					state_.result.contacts_summaries_with_redundancy_in_periodic_box[cds.id_b].push_back(cds);
				}
			}
		}
		else
		{
			state_.result.contacts_summaries_with_redundancy_in_periodic_box.clear();
		}
	}

	void prepare_for_possible_init_or_update(TimeRecorder& time_recorder) noexcept
	{
		time_recorder.reset();

		if(backup_enabled_ && !in_sync_with_backup_)
		{
			state_backup_.assign_to_apply_update(state_);
			in_sync_with_backup_=true;
		}

		time_recorder.record_elapsed_miliseconds_and_reset("backup state");

		state_.ids_of_changed_input_spheres.clear();
		state_.ids_of_affected_input_spheres.clear();
		state_.last_update_was_full_reinit=false;
	}

	void update_using_current_state(TimeRecorder& time_recorder) noexcept
	{
		time_recorder.reset();

		if(involvement_of_spheres_for_update_.size()!=state_.spheres_container.input_spheres().size())
		{
			involvement_of_spheres_for_update_.clear();
			involvement_of_spheres_for_update_.resize(state_.spheres_container.input_spheres().size(), 0);
		}

		for(UnsignedInt i=0;i<state_.ids_of_affected_input_spheres.size();i++)
		{
			involvement_of_spheres_for_update_[state_.ids_of_affected_input_spheres[i]]=1;
		}

		{
			RadicalTessellation::ResultGraphics result_graphics;
			RadicalTessellation::construct_full_tessellation(state_.spheres_container, involvement_of_spheres_for_update_, std::vector<int>(), false, false, FLOATCONST(0.0), std::vector<Float>(), buffered_temporary_storage_.tessellation_result, result_graphics, time_recorder);

			{
				const ConditionToRemoveContact condition_to_remove_contact(involvement_of_spheres_for_update_);

				for(UnsignedInt i=0;i<state_.ids_of_affected_input_spheres.size();i++)
				{
					const UnsignedInt sphere_id=state_.ids_of_affected_input_spheres[i];
					{
						std::vector<RadicalTessellation::ContactDescriptorSummary>& v=state_.result.contacts_summaries[sphere_id];
						std::vector<RadicalTessellation::ContactDescriptorSummary>::iterator it=std::remove_if(v.begin(), v.end(), condition_to_remove_contact);
						v.erase(it, v.end());
					}
					if(!state_.result.contacts_summaries_with_redundancy_in_periodic_box.empty())
					{
						std::vector<RadicalTessellation::ContactDescriptorSummary>& v=state_.result.contacts_summaries_with_redundancy_in_periodic_box[sphere_id];
						std::vector<RadicalTessellation::ContactDescriptorSummary>::iterator it=std::remove_if(v.begin(), v.end(), condition_to_remove_contact);
						v.erase(it, v.end());
					}
				}
			}

			if(!buffered_temporary_storage_.tessellation_result.contacts_summaries.empty())
			{
				for(UnsignedInt i=0;i<buffered_temporary_storage_.tessellation_result.contacts_summaries.size();i++)
				{
					const RadicalTessellation::ContactDescriptorSummary& cds=buffered_temporary_storage_.tessellation_result.contacts_summaries[i];
					state_.result.contacts_summaries[cds.id_a].push_back(cds);
					state_.result.contacts_summaries[cds.id_b].push_back(cds);
				}

				if(!state_.result.contacts_summaries_with_redundancy_in_periodic_box.empty())
				{
					const std::vector<RadicalTessellation::ContactDescriptorSummary>& all_contacts_summaries=(buffered_temporary_storage_.tessellation_result.contacts_summaries_with_redundancy_in_periodic_box.empty() ? buffered_temporary_storage_.tessellation_result.contacts_summaries : buffered_temporary_storage_.tessellation_result.contacts_summaries_with_redundancy_in_periodic_box);

					for(UnsignedInt i=0;i<all_contacts_summaries.size();i++)
					{
						const RadicalTessellation::ContactDescriptorSummary& cds=all_contacts_summaries[i];
						if(cds.id_a<state_.result.contacts_summaries_with_redundancy_in_periodic_box.size())
						{
							state_.result.contacts_summaries_with_redundancy_in_periodic_box[cds.id_a].push_back(cds);
						}
						if(cds.id_b<state_.result.contacts_summaries_with_redundancy_in_periodic_box.size())
						{
							state_.result.contacts_summaries_with_redundancy_in_periodic_box[cds.id_b].push_back(cds);
						}
					}
				}
			}

			time_recorder.record_elapsed_miliseconds_and_reset("update contacts summaries");
		}

		{
			const std::vector< std::vector<RadicalTessellation::ContactDescriptorSummary> >& all_contacts_summaries=(state_.result.contacts_summaries_with_redundancy_in_periodic_box.empty() ? state_.result.contacts_summaries : state_.result.contacts_summaries_with_redundancy_in_periodic_box);

			for(UnsignedInt i=0;i<state_.ids_of_affected_input_spheres.size();i++)
			{
				const UnsignedInt sphere_id=state_.ids_of_affected_input_spheres[i];
				RadicalTessellation::CellContactDescriptorsSummary& cs=state_.result.cells_summaries[sphere_id];
				cs=RadicalTessellation::CellContactDescriptorsSummary();
				const std::vector<RadicalTessellation::ContactDescriptorSummary>& v=all_contacts_summaries[sphere_id];
				for(UnsignedInt j=0;j<v.size();j++)
				{
					cs.add(sphere_id, v[j]);
				}
				cs.compute_sas(state_.spheres_container.input_spheres()[sphere_id].r);
				if(cs.stage==0 && state_.spheres_container.all_exclusion_statuses()[sphere_id]==0 && state_.spheres_container.all_colliding_ids()[sphere_id].empty())
				{
					cs.compute_sas_detached(sphere_id, state_.spheres_container.input_spheres()[sphere_id].r);
				}
			}

			time_recorder.record_elapsed_miliseconds_and_reset("update cell summaries");
		}

		for(UnsignedInt i=0;i<state_.ids_of_affected_input_spheres.size();i++)
		{
			involvement_of_spheres_for_update_[state_.ids_of_affected_input_spheres[i]]=0;
		}
	}

	State state_;
	State state_backup_;
	bool backup_enabled_;
	bool in_sync_with_backup_;
	std::vector<int> involvement_of_spheres_for_update_;
	BufferedTemporaryStorage buffered_temporary_storage_;
};

}

#endif /* VORONOTALT_UPDATEABLE_RADICAL_TESSELLATION_H_ */
