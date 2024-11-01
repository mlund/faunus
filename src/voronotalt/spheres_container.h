#ifndef VORONOTALT_SPHERES_CONTAINER_H_
#define VORONOTALT_SPHERES_CONTAINER_H_

#include "spheres_searcher.h"
#include "periodic_box.h"
#include "time_recorder.h"

namespace voronotalt
{

class SpheresContainer
{
public:
	struct ResultOfPreparationForTessellation
	{
		std::vector< std::pair<UnsignedInt, UnsignedInt> > relevant_collision_ids;

		ResultOfPreparationForTessellation() noexcept
		{
		}
	};

	SpheresContainer() noexcept : total_collisions_(0)
	{
	}

	void init(const std::vector<SimpleSphere>& input_spheres, TimeRecorder& time_recorder) noexcept
	{
		init(input_spheres, PeriodicBox(), time_recorder);
	}

	void init(const std::vector<SimpleSphere>& input_spheres, const PeriodicBox& periodic_box, TimeRecorder& time_recorder) noexcept
	{
		time_recorder.reset();

		periodic_box_=periodic_box;
		input_spheres_=input_spheres;

		if(periodic_box_.enabled())
		{
			populated_spheres_.resize(input_spheres_.size()*27);
			std::vector<UnsignedInt> collected_indices;
			for(UnsignedInt i=0;i<input_spheres_.size();i++)
			{
				set_sphere_periodic_instances(i, false, collected_indices);
			}
		}
		else
		{
			populated_spheres_=input_spheres_;
		}

		all_exclusion_statuses_.resize(populated_spheres_.size(), 0);

		time_recorder.record_elapsed_miliseconds_and_reset("populate spheres");

		spheres_searcher_.init(populated_spheres_);

		time_recorder.record_elapsed_miliseconds_and_reset("init spheres searcher");

		all_colliding_ids_.resize(input_spheres_.size());

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<input_spheres_.size();i++)
			{
				all_colliding_ids_[i].reserve(100);
				spheres_searcher_.find_colliding_ids(i, all_colliding_ids_[i], true, all_exclusion_statuses_[i]);
			}
		}

		if(periodic_box_.enabled())
		{
			for(UnsignedInt i=0;i<input_spheres_.size();i++)
			{
				set_exclusion_status_periodic_instances(i);
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("detect all collisions");

		total_collisions_=0;

		for(UnsignedInt i=0;i<all_colliding_ids_.size();i++)
		{
			total_collisions_+=all_colliding_ids_[i].size();
		}

		total_collisions_=total_collisions_/2;

		time_recorder.record_elapsed_miliseconds_and_reset("count all collisions");
	}

	bool update(
			const std::vector<SimpleSphere>& new_input_spheres,
			const std::vector<UnsignedInt>& provided_ids_of_changed_input_spheres,
			const bool trust_provided_ids_of_changed_input_spheres,
			std::vector<UnsignedInt>& ids_of_changed_input_spheres,
			std::vector<UnsignedInt>& ids_of_affected_input_spheres,
			TimeRecorder& time_recorder) noexcept
	{
		time_recorder.reset();

		ids_of_changed_input_spheres.clear();
		ids_of_affected_input_spheres.clear();

		if(new_input_spheres.size()!=input_spheres_.size())
		{
			reinit(new_input_spheres, ids_of_changed_input_spheres, ids_of_affected_input_spheres, time_recorder);
			return true;
		}

		if(trust_provided_ids_of_changed_input_spheres)
		{
			ids_of_changed_input_spheres=provided_ids_of_changed_input_spheres;
		}
		else
		{
			for(UnsignedInt i=0;i<new_input_spheres.size();i++)
			{
				if(!sphere_equals_sphere(new_input_spheres[i], input_spheres_[i]))
				{
					if(ids_of_changed_input_spheres.size()<size_threshold_for_full_reinit())
					{
						ids_of_changed_input_spheres.push_back(i);
					}
					else
					{
						reinit(new_input_spheres, ids_of_changed_input_spheres, ids_of_affected_input_spheres, time_recorder);
						return true;
					}
				}
			}

			time_recorder.record_elapsed_miliseconds_and_reset("identify changed spheres ids for update");
		}

		if(ids_of_changed_input_spheres.empty())
		{
			return false;
		}

		if(ids_of_changed_input_spheres.size()>size_threshold_for_full_reinit())
		{
			reinit(new_input_spheres, ids_of_changed_input_spheres, ids_of_affected_input_spheres, time_recorder);
			return true;
		}

		for(UnsignedInt i=0;i<ids_of_changed_input_spheres.size();i++)
		{
			if(ids_of_changed_input_spheres[i]>=input_spheres_.size())
			{
				reinit(new_input_spheres, ids_of_changed_input_spheres, ids_of_affected_input_spheres, time_recorder);
				return true;
			}
		}

		{
			bool update_needed=false;
			for(UnsignedInt i=0;!update_needed && i<ids_of_changed_input_spheres.size();i++)
			{
				const UnsignedInt sphere_id=ids_of_changed_input_spheres[i];
				if(!sphere_equals_sphere(new_input_spheres[sphere_id], input_spheres_[sphere_id]))
				{
					update_needed=true;
				}
			}
			if(!update_needed)
			{
				return false;
			}
		}

		{
			ids_of_affected_input_spheres=ids_of_changed_input_spheres;
			std::sort(ids_of_affected_input_spheres.begin(), ids_of_affected_input_spheres.end());

			for(UnsignedInt i=0;i<ids_of_changed_input_spheres.size();i++)
			{
				const UnsignedInt sphere_id=ids_of_changed_input_spheres[i];
				for(UnsignedInt j=0;j<all_colliding_ids_[sphere_id].size();j++)
				{
					if(ids_of_affected_input_spheres.size()<size_threshold_for_full_reinit())
					{
						const UnsignedInt affected_sphere_id=all_colliding_ids_[sphere_id][j].index%input_spheres_.size();
						std::vector<UnsignedInt>::iterator it=std::lower_bound(ids_of_affected_input_spheres.begin(), ids_of_affected_input_spheres.end(), affected_sphere_id);
						if(it==ids_of_affected_input_spheres.end() || (*it)!=affected_sphere_id)
						{
							ids_of_affected_input_spheres.insert(it, affected_sphere_id);
						}
					}
					else
					{
						reinit(new_input_spheres, ids_of_changed_input_spheres, ids_of_affected_input_spheres, time_recorder);
						return true;
					}
				}
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("gather affected spheres ids for update");

		{
			if(periodic_box_.enabled())
			{
				std::vector<UnsignedInt> ids_of_changed_populated_spheres;
				ids_of_changed_populated_spheres.reserve(ids_of_changed_input_spheres.size()*27);
				for(UnsignedInt i=0;i<ids_of_changed_input_spheres.size();i++)
				{
					const UnsignedInt sphere_id=ids_of_changed_input_spheres[i];
					input_spheres_[sphere_id]=new_input_spheres[sphere_id];
					if(periodic_box_.enabled())
					{
						set_sphere_periodic_instances(sphere_id, true, ids_of_changed_populated_spheres);
					}
				}
				spheres_searcher_.update(populated_spheres_, ids_of_changed_populated_spheres);
			}
			else
			{
				for(UnsignedInt i=0;i<ids_of_changed_input_spheres.size();i++)
				{
					const UnsignedInt sphere_id=ids_of_changed_input_spheres[i];
					input_spheres_[sphere_id]=new_input_spheres[sphere_id];
					populated_spheres_[sphere_id]=input_spheres_[sphere_id];
				}
				spheres_searcher_.update(input_spheres_, ids_of_changed_input_spheres);
			}

			time_recorder.record_elapsed_miliseconds_and_reset("update spheres searcher");

			{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
				for(UnsignedInt i=0;i<ids_of_affected_input_spheres.size();i++)
				{
					const UnsignedInt sphere_id=ids_of_affected_input_spheres[i];
					all_colliding_ids_[sphere_id].clear();
					spheres_searcher_.find_colliding_ids(sphere_id, all_colliding_ids_[sphere_id], true, all_exclusion_statuses_[sphere_id]);
				}
			}

			if(periodic_box_.enabled())
			{
				for(UnsignedInt i=0;i<ids_of_affected_input_spheres.size();i++)
				{
					const UnsignedInt sphere_id=ids_of_affected_input_spheres[i];
					set_exclusion_status_periodic_instances(sphere_id);
				}
			}

			buffered_temporary_storage_.clear();

			for(UnsignedInt i=0;i<ids_of_changed_input_spheres.size();i++)
			{
				const UnsignedInt sphere_id=ids_of_changed_input_spheres[i];
				for(UnsignedInt j=0;j<all_colliding_ids_[sphere_id].size();j++)
				{
					const UnsignedInt affected_sphere_id=all_colliding_ids_[sphere_id][j].index%input_spheres_.size();
					std::vector<UnsignedInt>::iterator it=std::lower_bound(buffered_temporary_storage_.more_ids_of_affected_input_spheres.begin(), buffered_temporary_storage_.more_ids_of_affected_input_spheres.end(), affected_sphere_id);
					if(it==buffered_temporary_storage_.more_ids_of_affected_input_spheres.end() || (*it)!=affected_sphere_id)
					{
						if(!std::binary_search(ids_of_affected_input_spheres.begin(), ids_of_affected_input_spheres.end(), affected_sphere_id))
						{
							buffered_temporary_storage_.more_ids_of_affected_input_spheres.insert(it, affected_sphere_id);
						}
					}
				}
			}

			if(!buffered_temporary_storage_.more_ids_of_affected_input_spheres.empty())
			{
				{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
					for(UnsignedInt i=0;i<buffered_temporary_storage_.more_ids_of_affected_input_spheres.size();i++)
					{
						const UnsignedInt sphere_id=buffered_temporary_storage_.more_ids_of_affected_input_spheres[i];
						all_colliding_ids_[sphere_id].clear();
						spheres_searcher_.find_colliding_ids(sphere_id, all_colliding_ids_[sphere_id], true, all_exclusion_statuses_[sphere_id]);
					}
				}

				if(periodic_box_.enabled())
				{
					for(UnsignedInt i=0;i<ids_of_affected_input_spheres.size();i++)
					{
						const UnsignedInt sphere_id=ids_of_affected_input_spheres[i];
						set_exclusion_status_periodic_instances(sphere_id);
					}
				}

				buffered_temporary_storage_.merged_ids_of_affected_input_spheres.resize(ids_of_affected_input_spheres.size()+buffered_temporary_storage_.more_ids_of_affected_input_spheres.size());

				std::merge(ids_of_affected_input_spheres.begin(), ids_of_affected_input_spheres.end(),
						buffered_temporary_storage_.more_ids_of_affected_input_spheres.begin(), buffered_temporary_storage_.more_ids_of_affected_input_spheres.end(),
						buffered_temporary_storage_.merged_ids_of_affected_input_spheres.begin());

				ids_of_affected_input_spheres.swap(buffered_temporary_storage_.merged_ids_of_affected_input_spheres);
			}

			time_recorder.record_elapsed_miliseconds_and_reset("update relevant collisions");

			total_collisions_=0;

			for(UnsignedInt i=0;i<all_colliding_ids_.size();i++)
			{
				total_collisions_+=all_colliding_ids_[i].size();
			}

			total_collisions_=total_collisions_/2;

			time_recorder.record_elapsed_miliseconds_and_reset("recount all collisions");
		}

		return true;
	}

	bool update_by_setting_exclusion_status(const UnsignedInt id_of_masked_input_sphere, const bool new_exclusion_status) noexcept
	{
        if(id_of_masked_input_sphere>=input_spheres_.size() || id_of_masked_input_sphere>=all_exclusion_statuses_.size() || (new_exclusion_status ? all_exclusion_statuses_[id_of_masked_input_sphere]>0 : all_exclusion_statuses_[id_of_masked_input_sphere]<1))
		{
			return false;
		}

        all_exclusion_statuses_[id_of_masked_input_sphere]=(new_exclusion_status ? 1 : 0);

		if(periodic_box_.enabled())
		{
			set_exclusion_status_periodic_instances(id_of_masked_input_sphere);
		}

		return true;
	}

	bool update_by_setting_exclusion_status(const UnsignedInt id_of_masked_input_sphere, const bool new_exclusion_status, std::vector<UnsignedInt>& ids_of_changed_input_spheres, std::vector<UnsignedInt>& ids_of_affected_input_spheres) noexcept
	{
	    ids_of_changed_input_spheres.clear();
	    ids_of_affected_input_spheres.clear();

	    if(id_of_masked_input_sphere>=input_spheres_.size() || id_of_masked_input_sphere>=all_exclusion_statuses_.size() || (new_exclusion_status ? all_exclusion_statuses_[id_of_masked_input_sphere]>0 : all_exclusion_statuses_[id_of_masked_input_sphere]<1))
	    {
	        return false;
	    }

	    ids_of_changed_input_spheres.push_back(id_of_masked_input_sphere);
	    ids_of_affected_input_spheres.push_back(id_of_masked_input_sphere);

	    for(std::size_t j=0;j<all_colliding_ids_[id_of_masked_input_sphere].size();j++)
	    {
	        const UnsignedInt affected_sphere_id=all_colliding_ids_[id_of_masked_input_sphere][j].index%input_spheres_.size();
	        std::vector<UnsignedInt>::iterator it=std::lower_bound(ids_of_affected_input_spheres.begin(), ids_of_affected_input_spheres.end(), affected_sphere_id);
	        if(it==ids_of_affected_input_spheres.end() || (*it)!=affected_sphere_id)
	        {
	            ids_of_affected_input_spheres.insert(it, affected_sphere_id);
	        }
	    }

	    all_exclusion_statuses_[id_of_masked_input_sphere]=(new_exclusion_status ? 1 : 0);

	    if(periodic_box_.enabled())
	    {
	        set_exclusion_status_periodic_instances(id_of_masked_input_sphere);
	    }

	    return true;
	}

	void assign(const SpheresContainer& obj) noexcept
	{
		periodic_box_=obj.periodic_box_;
		total_collisions_=obj.total_collisions_;

		spheres_searcher_.assign(obj.spheres_searcher_);

		input_spheres_.resize(obj.input_spheres_.size());
		populated_spheres_.resize(obj.populated_spheres_.size());
		all_exclusion_statuses_.resize(obj.all_exclusion_statuses_.size());
		all_colliding_ids_.resize(obj.all_colliding_ids_.size());

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<obj.input_spheres_.size();i++)
			{
				input_spheres_[i]=obj.input_spheres_[i];
			}
		}

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<obj.populated_spheres_.size();i++)
			{
				populated_spheres_[i]=obj.populated_spheres_[i];
			}
		}

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<obj.all_exclusion_statuses_.size();i++)
			{
				all_exclusion_statuses_[i]=obj.all_exclusion_statuses_[i];
			}
		}

		{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
			for(UnsignedInt i=0;i<obj.all_colliding_ids_.size();i++)
			{
				all_colliding_ids_[i]=obj.all_colliding_ids_[i];
			}
		}
	}

	void assign(const SpheresContainer& obj, const std::vector<UnsignedInt>& subset_of_ids) noexcept
	{
		if(subset_of_ids.empty()
				|| obj.input_spheres_.empty()
				|| !periodic_box_.equals(obj.periodic_box_)
				|| input_spheres_.size()!=obj.input_spheres_.size()
				|| populated_spheres_.size()!=obj.populated_spheres_.size()
				|| all_exclusion_statuses_.size()!=obj.all_exclusion_statuses_.size()
				|| all_colliding_ids_.size()!=obj.all_colliding_ids_.size()
				|| subset_of_ids.size()>=size_threshold_for_full_reinit())
		{
			assign(obj);
		}
		else
		{
			for(UnsignedInt i=0;i<subset_of_ids.size();i++)
			{
				const UnsignedInt sphere_id=subset_of_ids[i];
				if(sphere_id>=input_spheres_.size())
				{
					assign(obj);
					return;
				}
			}

			periodic_box_=obj.periodic_box_;
			total_collisions_=obj.total_collisions_;

			spheres_searcher_.assign(obj.spheres_searcher_);

			{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
				for(UnsignedInt i=0;i<subset_of_ids.size();i++)
				{
					const UnsignedInt sphere_id=subset_of_ids[i];
					input_spheres_[sphere_id]=obj.input_spheres_[sphere_id];
					populated_spheres_[sphere_id]=obj.populated_spheres_[sphere_id];
					all_exclusion_statuses_[sphere_id]=obj.all_exclusion_statuses_[sphere_id];
					all_colliding_ids_[sphere_id]=obj.all_colliding_ids_[sphere_id];
				}
			}

			if(periodic_box_.enabled() && populated_spheres_.size()==input_spheres_.size()*27 && all_exclusion_statuses_.size()==all_exclusion_statuses_.size()*27)
			{
#ifdef VORONOTALT_OPENMP
#pragma omp parallel for
#endif
				for(UnsignedInt i=0;i<subset_of_ids.size();i++)
				{
					const UnsignedInt sphere_id=subset_of_ids[i];
					for(UnsignedInt m=1;m<27;m++)
					{
						const UnsignedInt shifted_sphere_id=(m*input_spheres_.size()+sphere_id);
						populated_spheres_[shifted_sphere_id]=obj.populated_spheres_[sphere_id];
						all_exclusion_statuses_[shifted_sphere_id]=all_exclusion_statuses_[sphere_id];
					}
				}
			}
		}
	}

	const PeriodicBox& periodic_box() const noexcept
	{
		return periodic_box_;
	}

	const std::vector<SimpleSphere>& input_spheres() const noexcept
	{
		return input_spheres_;
	}

	const std::vector<SimpleSphere>& populated_spheres() const noexcept
	{
		return populated_spheres_;
	}

	const std::vector<int>& all_exclusion_statuses() const noexcept
	{
		return all_exclusion_statuses_;
	}

	const std::vector< std::vector<ValuedID> >& all_colliding_ids() const noexcept
	{
		return all_colliding_ids_;
	}

	UnsignedInt total_collisions() const noexcept
	{
		return total_collisions_;
	}

	bool prepare_for_tessellation(const std::vector<int>& grouping_of_spheres, ResultOfPreparationForTessellation& result, TimeRecorder& time_recorder) const noexcept
	{
		return prepare_for_tessellation(std::vector<int>(), grouping_of_spheres, result, time_recorder);
	}

	bool prepare_for_tessellation(const std::vector<int>& involvement_of_spheres, const std::vector<int>& grouping_of_spheres, ResultOfPreparationForTessellation& result, TimeRecorder& time_recorder) const noexcept
	{
		time_recorder.reset();

		result.relevant_collision_ids.clear();

		if(populated_spheres_.empty())
		{
			return false;
		}

		result.relevant_collision_ids.reserve(total_collisions_);

		for(UnsignedInt id_a=0;id_a<input_spheres_.size();id_a++)
		{
			if((involvement_of_spheres.empty() || id_a>=involvement_of_spheres.size() || involvement_of_spheres[id_a]>0) && all_exclusion_statuses_[id_a]==0)
			{
				for(UnsignedInt j=0;j<all_colliding_ids_[id_a].size();j++)
				{
					const UnsignedInt id_b=all_colliding_ids_[id_a][j].index;
					const UnsignedInt id_b_canonical=(id_b%input_spheres_.size());
					if((involvement_of_spheres.empty() || id_b_canonical>=involvement_of_spheres.size() || involvement_of_spheres[id_b_canonical]>0) && all_exclusion_statuses_[id_b_canonical]==0)
					{
						if(id_b>=input_spheres_.size() || (all_colliding_ids_[id_a].size()<all_colliding_ids_[id_b_canonical].size()) || (id_a<id_b && all_colliding_ids_[id_a].size()==all_colliding_ids_[id_b_canonical].size()))
						{
							if(grouping_of_spheres.empty() || id_a>=grouping_of_spheres.size() || id_b_canonical>=grouping_of_spheres.size() || grouping_of_spheres[id_a]!=grouping_of_spheres[id_b_canonical])
							{
								result.relevant_collision_ids.push_back(std::pair<UnsignedInt, UnsignedInt>(id_a, id_b));
							}
						}
					}
				}
			}
		}

		time_recorder.record_elapsed_miliseconds_and_reset("collect relevant collision indices");

		return true;
	}

private:
	struct BufferedTemporaryStorage
	{
		std::vector<UnsignedInt> more_ids_of_affected_input_spheres;
		std::vector<UnsignedInt> merged_ids_of_affected_input_spheres;

		void clear() noexcept
		{
			more_ids_of_affected_input_spheres.clear();
			merged_ids_of_affected_input_spheres.clear();
		}
	};

	UnsignedInt size_threshold_for_full_reinit() const noexcept
	{
		return static_cast<UnsignedInt>(input_spheres_.size()/2);
	}

	void reinit(const std::vector<SimpleSphere>& new_input_spheres, std::vector<UnsignedInt>& ids_of_changed_input_spheres, std::vector<UnsignedInt>& ids_of_affected_input_spheres, TimeRecorder& time_recorder) noexcept
	{
		init(new_input_spheres, periodic_box_, time_recorder);
		ids_of_changed_input_spheres.clear();
		ids_of_affected_input_spheres.clear();
	}

	void set_sphere_periodic_instances(const UnsignedInt i, const bool collect_indices, std::vector<UnsignedInt>& collected_indices) noexcept
	{
		if(i<input_spheres_.size())
		{
			const SimpleSphere& o=input_spheres_[i];
			if(!periodic_box_.enabled())
			{
				if(populated_spheres_.size()!=input_spheres_.size())
				{
					populated_spheres_.resize(input_spheres_.size());
				}
				populated_spheres_[i]=o;
				if(collect_indices)
				{
					collected_indices.push_back(i);
				}
			}
			else
			{
				if(populated_spheres_.size()!=(input_spheres_.size()*27))
				{
					populated_spheres_.resize(input_spheres_.size()*27);
				}
				populated_spheres_[i]=o;
				if(collect_indices)
				{
					collected_indices.push_back(i);
				}
				UnsignedInt g=1;
				for(int sx=-1;sx<=1;sx++)
				{
					for(int sy=-1;sy<=1;sy++)
					{
						for(int sz=-1;sz<=1;sz++)
						{
							if(sx!=0 || sy!=0 || sz!=0)
							{
								const UnsignedInt mi=(g*input_spheres_.size()+i);
								populated_spheres_[mi]=periodic_box_.shift_by_weighted_directions(o, static_cast<Float>(sx), static_cast<Float>(sy), static_cast<Float>(sz));
								if(collect_indices)
								{
									collected_indices.push_back(mi);
								}
								g++;
							}
						}
					}
				}
			}
		}
	}

	void set_exclusion_status_periodic_instances(const UnsignedInt i) noexcept
	{
		if(i<input_spheres_.size() && all_exclusion_statuses_.size()==input_spheres_.size()*27)
		{
			for(UnsignedInt m=1;m<27;m++)
			{
				all_exclusion_statuses_[m*input_spheres_.size()+i]=all_exclusion_statuses_[i];
			}
		}
	}

	PeriodicBox periodic_box_;
	std::vector<SimpleSphere> input_spheres_;
	std::vector<SimpleSphere> populated_spheres_;
	std::vector<int> all_exclusion_statuses_;
	SpheresSearcher spheres_searcher_;
	std::vector< std::vector<ValuedID> > all_colliding_ids_;
	UnsignedInt total_collisions_;
	BufferedTemporaryStorage buffered_temporary_storage_;
};

}

#endif /* VORONOTALT_SPHERES_CONTAINER_H_ */
