#ifndef VORONOTALT_SIMPLIFIED_AW_TESSELLATION_CONTACT_CONSTRUCTION_H_
#define VORONOTALT_SIMPLIFIED_AW_TESSELLATION_CONTACT_CONSTRUCTION_H_

#include <vector>
#include <list>
#include <algorithm>

#include "basic_types_and_functions.h"

namespace voronotalt
{

class SimplifiedAWTessellationContactConstruction
{
public:
	struct ContourPoint
	{
		SimplePoint p;
		UnsignedInt left_id;
		UnsignedInt right_id;

		ContourPoint(const SimplePoint& p, const UnsignedInt left_id, const UnsignedInt right_id) noexcept : p(p), left_id(left_id), right_id(right_id)
		{
		}
	};

	typedef std::list<ContourPoint> Contour;

	struct NeighborDescriptor
	{
		Float sort_value;
		UnsignedInt neighbor_id;

		NeighborDescriptor() noexcept : sort_value(FLOATCONST(0.0)), neighbor_id(0)
		{
		}

		bool operator<(const NeighborDescriptor& d) const noexcept
		{
			return (sort_value<d.sort_value || (sort_value==d.sort_value && neighbor_id<d.neighbor_id));
		}
	};

	struct ContourGraphics
	{
		std::vector<SimplePoint> outer_points;
		SimplePoint barycenter;
	};

	struct ContactDescriptorGraphics
	{
		std::vector<ContourGraphics> contours_graphics;

		ContactDescriptorGraphics() noexcept
		{
		}

		void clear() noexcept
		{
			contours_graphics.clear();
		}
	};

	struct ContactDescriptor
	{
		std::list<Contour> contours;
		std::vector<NeighborDescriptor> neighbor_descriptors;
		ContactDescriptorGraphics graphics;
		SimpleSphere intersection_circle_sphere;
		SimplePoint intersection_circle_axis;
		Float area;
		Float arc_length;
		Float distance;
		UnsignedInt id_a;
		UnsignedInt id_b;

		ContactDescriptor() noexcept :
			area(FLOATCONST(0.0)),
			arc_length(FLOATCONST(0.0)),
			distance(FLOATCONST(0.0)),
			id_a(0),
			id_b(0)
		{
		}

		void clear() noexcept
		{
			id_a=0;
			id_b=0;
			neighbor_descriptors.clear();
			contours.clear();
			graphics.clear();
			area=FLOATCONST(0.0);
			arc_length=FLOATCONST(0.0);
			distance=FLOATCONST(0.0);
		}
	};

	static bool construct_contact_descriptor(
			const std::vector<SimpleSphere>& spheres,
			const std::vector<int>& spheres_exclusion_statuses,
			const UnsignedInt a_id,
			const UnsignedInt b_id,
			const std::vector<ValuedID>& a_neighbor_collisions,
			const Float step,
			const int projections,
			const bool record_graphics,
			ContactDescriptor& result_contact_descriptor) noexcept
	{
		result_contact_descriptor.clear();
		if(a_id<spheres.size() && b_id<spheres.size())
		{
			result_contact_descriptor.id_a=a_id;
			result_contact_descriptor.id_b=b_id;
			const SimpleSphere& a=spheres[a_id];
			const SimpleSphere& b=spheres[b_id];
			if(sphere_intersects_sphere(a, b) && !sphere_contains_sphere(a, b) && !sphere_contains_sphere(b, a))
			{
				result_contact_descriptor.intersection_circle_sphere=intersection_circle_of_two_spheres(a, b);
				if(result_contact_descriptor.intersection_circle_sphere.r>FLOATCONST(0.0))
				{
					bool discarded=false;
					{
						for(UnsignedInt i=0;i<a_neighbor_collisions.size() && !discarded;i++)
						{
							const UnsignedInt neighbor_id=a_neighbor_collisions[i].index;
							if(neighbor_id!=b_id && (neighbor_id>=spheres_exclusion_statuses.size() || spheres_exclusion_statuses[neighbor_id]==0))
							{
								const SimpleSphere& c=spheres[neighbor_id];
								if(sphere_intersects_sphere(result_contact_descriptor.intersection_circle_sphere, c) && sphere_intersects_sphere(b, c))
								{
									if(sphere_contains_sphere(c, a) || sphere_contains_sphere(c, b))
									{
										discarded=true;
									}
									else
									{
										const SimplePoint nd_ac_plane_center=center_of_intersection_circle_of_two_spheres(a, c);
										const Float dist_from_a_center_to_c_center=distance_from_point_to_point(a.p, c.p);
										const SimplePoint nd_ac_plane_normal=point_and_number_product(sub_of_points(c.p, a.p), FLOATCONST(1.0)/dist_from_a_center_to_c_center);
										const Float dist_from_a_center_to_nd_ac_plane_center=distance_from_point_to_point(a.p, nd_ac_plane_center);
										const Float cos_val=dot_product(unit_point(sub_of_points(result_contact_descriptor.intersection_circle_sphere.p, a.p)), point_and_number_product(sub_of_points(nd_ac_plane_center, a.p), FLOATCONST(1.0)/dist_from_a_center_to_nd_ac_plane_center));
										const int hsi=halfspace_of_point(nd_ac_plane_center, nd_ac_plane_normal, result_contact_descriptor.intersection_circle_sphere.p);
										if(std::abs(cos_val)<FLOATCONST(1.0))
										{
											const Float l=std::abs(signed_distance_from_point_to_plane(nd_ac_plane_center, nd_ac_plane_normal, result_contact_descriptor.intersection_circle_sphere.p));
											const Float xl=l/std::sqrt(1-(cos_val*cos_val));
											const Float dist_from_a_to_midpoint=a.r+((dist_from_a_center_to_c_center-a.r-c.r)/FLOATCONST(2.0));
											const Float xl_buffer=std::abs(dist_from_a_to_midpoint-dist_from_a_center_to_nd_ac_plane_center);
											if(xl>=(result_contact_descriptor.intersection_circle_sphere.r+xl_buffer))
											{
												if(hsi>=0)
												{
													discarded=true;
												}
											}
											else
											{
												NeighborDescriptor nd;
												nd.neighbor_id=neighbor_id;
												nd.sort_value=(hsi>0 ? (FLOATCONST(0.0)-xl) : xl);
												result_contact_descriptor.neighbor_descriptors.push_back(nd);
											}
										}
										else
										{
											if(hsi>0)
											{
												discarded=true;
											}
										}
									}
								}
							}
						}
					}
					if(!discarded)
					{
						result_contact_descriptor.intersection_circle_axis=unit_point(sub_of_points(b.p, a.p));

						{
							Contour initial_contour;
							init_contour_from_base_and_axis(a_id, result_contact_descriptor.intersection_circle_sphere, result_contact_descriptor.intersection_circle_axis, step, initial_contour);
							if(!initial_contour.empty())
							{
								result_contact_descriptor.contours.push_back(initial_contour);
							}
						}

						if(!result_contact_descriptor.neighbor_descriptors.empty())
						{
							if(!result_contact_descriptor.contours.empty())
							{
								std::sort(result_contact_descriptor.neighbor_descriptors.begin(), result_contact_descriptor.neighbor_descriptors.end());
								for(UnsignedInt i=0;i<result_contact_descriptor.neighbor_descriptors.size();i++)
								{
									const UnsignedInt c_id=result_contact_descriptor.neighbor_descriptors[i].neighbor_id;
									const SimpleSphere& c=spheres[c_id];
									std::list<Contour>::iterator jt=result_contact_descriptor.contours.begin();
									while(jt!=result_contact_descriptor.contours.end())
									{
										Contour& contour=(*jt);
										std::list<Contour> segments;
										if(cut_and_split_contour(a, c, c_id, contour, segments))
										{
											if(!contour.empty())
											{
												mend_contour(a, b, c, c_id, step, projections, contour);
												if(check_contour_intersects_sphere(result_contact_descriptor.intersection_circle_sphere, contour))
												{
													++jt;
												}
												else
												{
													jt=result_contact_descriptor.contours.erase(jt);
												}
											}
											else
											{
												if(!segments.empty())
												{
													for(std::list<Contour>::iterator st=segments.begin();st!=segments.end();++st)
													{
														mend_contour(a, b, c, c_id, step, projections, (*st));
													}
													filter_contours_intersecting_sphere(result_contact_descriptor.intersection_circle_sphere, segments);
													if(!segments.empty())
													{
														result_contact_descriptor.contours.splice(jt, segments);
													}
												}
												jt=result_contact_descriptor.contours.erase(jt);
											}
										}
										else
										{
											++jt;
										}
									}
								}
							}

							if(!result_contact_descriptor.contours.empty())
							{
								const Float tolerated_deviation=FLOATCONST(0.5);
								bool strangely_extended=false;
								for(std::list<Contour>::const_iterator it=result_contact_descriptor.contours.begin();it!=result_contact_descriptor.contours.end() && !strangely_extended;++it)
								{
									const Contour& contour=(*it);
									for(Contour::const_iterator jt=contour.begin();jt!=contour.end() && !strangely_extended;++jt)
									{
										strangely_extended=strangely_extended || (distance_from_point_to_point(jt->p, a.p)-a.r>tolerated_deviation);
										strangely_extended=strangely_extended || (distance_from_point_to_point(jt->p, b.p)-b.r>tolerated_deviation);
									}
								}
								if(strangely_extended)
								{
									SimplePoint safe_center=point_and_number_product(sum_of_points(a.p, b.p), FLOATCONST(0.5));
									std::list<Contour> forcibly_shrunk_result=result_contact_descriptor.contours;
									for(std::list<Contour>::iterator it=forcibly_shrunk_result.begin();it!=forcibly_shrunk_result.end();++it)
									{
										Contour& contour=(*it);
										for(Contour::iterator jt=contour.begin();jt!=contour.end();++jt)
										{
											if((distance_from_point_to_point(jt->p, a.p)-a.r>tolerated_deviation) || (distance_from_point_to_point(jt->p, b.p)-b.r>tolerated_deviation))
											{
												jt->p=sum_of_points(safe_center, point_and_number_product(unit_point(sub_of_points(jt->p, safe_center)), std::min(a.r, b.r)));
											}
										}
									}
									result_contact_descriptor.contours.swap(forcibly_shrunk_result);
								}
							}
						}

						if(!result_contact_descriptor.contours.empty())
						{
							for(std::list<Contour>::const_iterator it=result_contact_descriptor.contours.begin();it!=result_contact_descriptor.contours.end();++it)
							{
								const Contour& contour=(*it);
								SimplePoint barycenter;
								result_contact_descriptor.area+=calculate_area_from_contour(contour, a, b, barycenter);
								result_contact_descriptor.arc_length+=calculate_arc_length_from_contour(a_id, contour);
								if(record_graphics && result_contact_descriptor.area>FLOATCONST(0.0))
								{
									result_contact_descriptor.graphics.contours_graphics.push_back(ContourGraphics());
									result_contact_descriptor.graphics.contours_graphics.back().barycenter=barycenter;
									collect_points_from_contour(contour, result_contact_descriptor.graphics.contours_graphics.back().outer_points);
								}
							}
						}

						result_contact_descriptor.distance=distance_from_point_to_point(a.p, b.p);
					}
				}
			}
		}
		return (result_contact_descriptor.area>FLOATCONST(0.0));
	}

private:
	class HyperboloidBetweenTwoSpheres
	{
	public:
		static inline SimplePoint project_point_on_hyperboloid(const SimplePoint& p, const SimpleSphere& s1, const SimpleSphere& s2) noexcept
		{
			if(s1.r>s2.r)
			{
				return project_point_on_hyperboloid(p, s2, s1);
			}
			else
			{
				const SimplePoint dv=point_and_number_product(sub_of_points(s1.p, s2.p), FLOATCONST(0.5));
				const SimplePoint dv_unit=unit_point(dv);
				const SimplePoint c=sum_of_points(s2.p, dv);
				const SimplePoint cp=sub_of_points(p, c);
				const Float lz=dot_product(dv_unit, cp);
				const Float lx=std::sqrt(std::max(squared_point_module(cp)-(lz*lz), FLOATCONST(0.0)));
				const Float z=project_point_on_simple_hyperboloid(lx, 0, point_module(dv), s1.r, s2.r);
				return sum_of_points(sub_of_points(p, point_and_number_product(dv_unit, lz)), point_and_number_product(dv_unit, z));
			}
		}

		static inline Float intersect_vector_with_hyperboloid(const SimplePoint& a, const SimplePoint& b, const SimpleSphere& s1, const SimpleSphere& s2) noexcept
		{
			if(s1.r>s2.r)
			{
				return intersect_vector_with_hyperboloid(a, b, s2, s1);
			}
			else
			{
				const SimplePoint dv=point_and_number_product(sub_of_points(s1.p, s2.p), FLOATCONST(0.5));
				const SimplePoint dv_unit=unit_point(dv);
				const SimplePoint c=sum_of_points(s2.p, dv);

				const SimplePoint ca=sub_of_points(a, c);
				const Float maz=dot_product(dv_unit, ca);
				const SimplePoint cax=sub_of_points(sub_of_points(a, point_and_number_product(dv_unit, maz)), c);
				const SimplePoint cax_unit=unit_point(cax);
				const Float max=dot_product(cax_unit, ca);
				const Float may=FLOATCONST(0.0);

				const SimplePoint cb=sub_of_points(b, c);
				const Float mbz=dot_product(dv_unit, cb);
				const Float mbx=dot_product(cax_unit, cb);
				const Float mby=std::sqrt(std::max(squared_point_module(cb)-mbz*mbz-mbx*mbx, FLOATCONST(0.0)));

				return intersect_vector_with_simple_hyperboloid(SimplePoint(max, may, maz), SimplePoint(mbx, mby, mbz), point_module(dv), s1.r, s2.r);
			}
		}

	private:
		static inline Float project_point_on_simple_hyperboloid(const Float x, const Float y, const Float d, const Float r1, const Float r2) noexcept
		{
			if(r1>r2)
			{
				return project_point_on_simple_hyperboloid(x, y, d, r2, r1);
			}
			else
			{
				const Float r=r2-r1;
				return 2*r*std::sqrt(std::max((0-r*r+4*d*d)*(4*x*x+4*y*y+4*d*d-r*r), FLOATCONST(0.0)))/(0-4*r*r+16*d*d);
			}
		}

		static inline Float intersect_vector_with_simple_hyperboloid(const SimplePoint& a, const SimplePoint& b, const Float d, const Float r1, const Float r2) noexcept
		{
			if(r1>r2)
			{
				return intersect_vector_with_simple_hyperboloid(a, b, d, r2, r1);
			}
			else
			{
				const Float r=r2-r1;
				SimplePoint ab=sub_of_points(b, a);
				SimplePoint v=unit_point(ab);
				const Float k=(4*r*r/((0-4*r*r+16*d*d)*(0-4*r*r+16*d*d))) * (0-r*r+4*d*d) * 4;
				const Float m=(4*d*d-r*r)*k/4;

				const Float x0=a.x;
				const Float y0=a.y;
				const Float z0=a.z;
				const Float vx=v.x;
				const Float vy=v.y;
				const Float vz=v.z;

				const Float t1 =  (std::sqrt((k*vy*vy+k*vx*vx)*z0*z0+(-2*k*vy*vz*y0-2*k*vx*vz*x0)*z0+(k*vz*vz-k*k*vx*vx)*y0*y0+2*k*k*vx*vy*x0*y0+(k*vz*vz-k*k*vy*vy)*x0*x0+m*vz*vz-k*m*vy*vy-k*m*vx*vx)-vz*z0+k*vy*y0+k*vx*x0)/(vz*vz-k*vy*vy-k*vx*vx);

				const Float t2 = -(std::sqrt((k*vy*vy+k*vx*vx)*z0*z0+(-2*k*vy*vz*y0-2*k*vx*vz*x0)*z0+(k*vz*vz-k*k*vx*vx)*y0*y0+2*k*k*vx*vy*x0*y0+(k*vz*vz-k*k*vy*vy)*x0*x0+m*vz*vz-k*m*vy*vy-k*m*vx*vx)+vz*z0-k*vy*y0-k*vx*x0)/(vz*vz-k*vy*vy-k*vx*vx);

				const SimplePoint tp1=sum_of_points(a, point_and_number_product(v, t1));
				const SimplePoint tp2=sum_of_points(a, point_and_number_product(v, t2));
				if(greater(t1, FLOATCONST(0.0)) && less(t1, point_module(ab)) && equal(tp1.z, std::sqrt(k*tp1.x*tp1.x+k*tp1.y*tp1.y+m), FLOATCONST(0.000001)))
				{
					return t1;
				}
				else if(greater(t2, FLOATCONST(0.0)) && less(t2, point_module(ab)) && equal(tp2.z, std::sqrt(k*tp2.x*tp2.x+k*tp2.y*tp2.y+m), FLOATCONST(0.000001)))
				{
					return t2;
				}
				else
				{
					return FLOATCONST(0.0);
				}
			}
		}
	};

	static void init_contour_from_base_and_axis(
			const UnsignedInt a_id,
			const SimpleSphere& base,
			const SimplePoint& axis,
			const Float length_step,
			Contour& result) noexcept
	{
		const Float angle_step=std::max(std::min(length_step/base.r, PIVALUE/FLOATCONST(3.0)), PIVALUE/FLOATCONST(36.0));
		const SimplePoint first_point=point_and_number_product(any_normal_of_vector(axis), base.r);
		result.push_back(ContourPoint(sum_of_points(base.p, first_point), a_id, a_id));
		for(Float rotation_angle=angle_step;rotation_angle<(PIVALUE*FLOATCONST(2.0));rotation_angle+=angle_step)
		{
			result.push_back(ContourPoint(sum_of_points(base.p, rotate_point_around_axis(axis, rotation_angle, first_point)), a_id, a_id));
		}
	}

	static bool cut_and_split_contour(
			const SimpleSphere& a,
			const SimpleSphere& c,
			const UnsignedInt c_id,
			Contour& contour,
			std::list<Contour>& segments) noexcept
	{
		const UnsignedInt outsiders_count=mark_contour(a, c, c_id, contour);
		if(outsiders_count>0)
		{
			if(outsiders_count<contour.size())
			{
				std::list<Contour::iterator> cuts;
				const int cuts_count=cut_contour(a, c, c_id, contour, cuts);
				if(cuts_count>0 && cuts_count%2==0)
				{
					if(cuts_count>2)
					{
						order_cuts(cuts);
						split_contour(contour, cuts, segments);
					}
				}
			}
			else
			{
				contour.clear();
			}
			return true;
		}
		return false;
	}

	static UnsignedInt mark_contour(
			const SimpleSphere& a,
			const SimpleSphere& c,
			const UnsignedInt c_id,
			Contour& contour) noexcept
	{
		UnsignedInt outsiders_count=0;
		for(Contour::iterator it=contour.begin();it!=contour.end();++it)
		{
			if((distance_from_point_to_point(it->p, c.p)-c.r)<(distance_from_point_to_point(it->p, a.p)-a.r))
			{
				it->left_id=c_id;
				it->right_id=c_id;
				outsiders_count++;
			}
		}
		return outsiders_count;
	}

	static UnsignedInt cut_contour(
			const SimpleSphere& a,
			const SimpleSphere& c,
			const UnsignedInt c_id,
			Contour& contour,
			std::list<Contour::iterator>& cuts) noexcept
	{
		UnsignedInt cuts_count=0;
		Contour::iterator it=contour.begin();
		while(it!=contour.end())
		{
			if(it->left_id==c_id && it->right_id==c_id)
			{
				const Contour::iterator left_it=get_left_iterator(contour, it);
				const Contour::iterator right_it=get_right_iterator(contour, it);

				if(left_it->right_id!=c_id)
				{
					const SimplePoint& p0=it->p;
					const SimplePoint& p1=left_it->p;
					const Float l=HyperboloidBetweenTwoSpheres::intersect_vector_with_hyperboloid(p0, p1, a, c);
					cuts.push_back(contour.insert(it, ContourPoint(sum_of_points(p0, point_and_number_product(unit_point(sub_of_points(p1, p0)), l)), left_it->right_id, it->left_id)));
					cuts_count++;
				}

				if(right_it->left_id!=c_id)
				{
					const SimplePoint& p0=it->p;
					const SimplePoint& p1=right_it->p;
					const Float l=HyperboloidBetweenTwoSpheres::intersect_vector_with_hyperboloid(p0, p1, a, c);
					cuts.push_back(contour.insert(right_it, ContourPoint(sum_of_points(p0, point_and_number_product(unit_point(sub_of_points(p1, p0)), l)), it->right_id, right_it->left_id)));
					cuts_count++;
				}

				it=contour.erase(it);
			}
			else
			{
				++it;
			}
		}
		return cuts_count;
	}

	static void order_cuts(std::list<Contour::iterator>& cuts) noexcept
	{
		Float sums[2]={FLOATCONST(0.0), FLOATCONST(0.0)};
		for(int i=0;i<2;i++)
		{
			if(i==1)
			{
				shift_list(cuts, false);
			}
			std::list<Contour::iterator>::const_iterator it=cuts.begin();
			while(it!=cuts.end())
			{
				std::list<Contour::iterator>::const_iterator next=it;
				++next;
				if(next!=cuts.end())
				{
					sums[i]+=distance_from_point_to_point((*it)->p, (*next)->p);
					it=next;
					++it;
				}
				else
				{
					it=cuts.end();
				}
			}
		}
		if(sums[0]<sums[1])
		{
			shift_list(cuts, true);
		}
	}

	static UnsignedInt split_contour(
			Contour& contour,
			const std::list<Contour::iterator>& ordered_cuts,
			std::list<Contour>& segments) noexcept
	{
		UnsignedInt segments_count=0;
		std::list<Contour::iterator>::const_iterator it=ordered_cuts.begin();
		while(it!=ordered_cuts.end())
		{
			std::list<Contour::iterator>::const_iterator next=it;
			++next;
			if(next!=ordered_cuts.end())
			{
				if((*next)!=get_right_iterator(contour, (*it)))
				{
					segments_count++;
					Contour segment;
					Contour::iterator jt=(*it);
					do
					{
						segment.push_back(*jt);
						jt=get_right_iterator(contour, jt);
					}
					while(jt!=(*next));
					segments.push_back(segment);
				}
				it=next;
				++it;
			}
			else
			{
				it=ordered_cuts.end();
			}
		}
		if(segments_count>0)
		{
			contour.clear();
		}
		return segments_count;
	}

	static void mend_contour(
			const SimpleSphere& a,
			const SimpleSphere& b,
			const SimpleSphere& c,
			const UnsignedInt c_id,
			const Float step,
			const int projections,
			Contour& contour) noexcept
	{
		Contour::iterator it=contour.begin();
		while(it!=contour.end())
		{
			if(it->left_id!=c_id && it->right_id==c_id)
			{
				const Contour::iterator jt=get_right_iterator(contour, it);
				if(jt->left_id==c_id)
				{
					SimplePoint& p0=it->p;
					SimplePoint& p1=jt->p;
					for(int e=0;e<projections;e++)
					{
						p0=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p0, b, c);
						p0=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p0, a, c);
						p0=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p0, a, b);
						p1=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p1, b, c);
						p1=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p1, a, c);
						p1=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p1, a, b);
					}
					const Float distance=distance_from_point_to_point(p0, p1);
					if(distance>step)
					{
						const int leap_distance=static_cast<int>(std::floor(distance/step+FLOATCONST(0.5)));
						const Float leap_size=distance/static_cast<Float>(leap_distance);
						const SimplePoint direction=unit_point(sub_of_points(p1, p0));
						for(int leap=1;leap<leap_distance;leap++)
						{
							SimplePoint p=sum_of_points(p0, point_and_number_product(direction, (leap_size*static_cast<Float>(leap))));
							for(int e=0;e<projections;e++)
							{
								p=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p, b, c);
								p=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p, a, c);
								p=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(p, a, b);
							}
							contour.insert(jt, ContourPoint(p, c_id, c_id));
						}
					}
				}
			}
			++it;
		}
	}

	static bool check_contour_intersects_sphere(const SimpleSphere& shell, const Contour& contour) noexcept
	{
		for(Contour::const_iterator it=contour.begin();it!=contour.end();++it)
		{
			if(distance_from_point_to_point(shell.p, it->p)<=shell.r)
			{
				return true;
			}
		}
		return false;
	}

	static void filter_contours_intersecting_sphere(const SimpleSphere& shell, std::list<Contour>& contours) noexcept
	{
		std::list<Contour>::iterator it=contours.begin();
		while(it!=contours.end())
		{
			if(!check_contour_intersects_sphere(shell, (*it)))
			{
				it=contours.erase(it);
			}
			else
			{
				++it;
			}
		}
	}

	template<typename List, typename Iterator>
	static Iterator get_left_iterator(List& container, const Iterator& iterator) noexcept
	{
		Iterator left_it=iterator;
		if(left_it==container.begin())
		{
			left_it=container.end();
		}
		--left_it;
		return left_it;
	}

	template<typename List, typename Iterator>
	static Iterator get_right_iterator(List& container, const Iterator& iterator) noexcept
	{
		Iterator right_it=iterator;
		++right_it;
		if(right_it==container.end())
		{
			right_it=container.begin();
		}
		return right_it;
	}

	template<typename List>
	static void shift_list(List& list, const bool reverse) noexcept
	{
		if(!reverse)
		{
			list.push_front(*list.rbegin());
			list.pop_back();
		}
		else
		{
			list.push_back(*list.begin());
			list.pop_front();
		}
	}

	static bool collect_points_from_contour(const Contour& contour, std::vector<SimplePoint>& contour_points) noexcept
	{
		contour_points.clear();
		if(!contour.empty())
		{
			contour_points.reserve(contour.size());
			for(Contour::const_iterator jt=contour.begin();jt!=contour.end();++jt)
			{
				contour_points.push_back(jt->p);
			}
		}
		return (!contour_points.empty());
	}

	static Float calculate_area_from_contour(const Contour& contour, const SimpleSphere& sphere1, const SimpleSphere& sphere2, SimplePoint& contour_barycenter) noexcept
	{
		Float area=FLOATCONST(0.0);
		if(!contour.empty())
		{
			contour_barycenter.x=FLOATCONST(0.0);
			contour_barycenter.y=FLOATCONST(0.0);
			contour_barycenter.z=FLOATCONST(0.0);
			for(Contour::const_iterator jt=contour.begin();jt!=contour.end();++jt)
			{
				contour_barycenter.x+=jt->p.x;
				contour_barycenter.y+=jt->p.y;
				contour_barycenter.z+=jt->p.z;
			}
			contour_barycenter.x/=static_cast<Float>(contour.size());
			contour_barycenter.y/=static_cast<Float>(contour.size());
			contour_barycenter.z/=static_cast<Float>(contour.size());

			contour_barycenter=HyperboloidBetweenTwoSpheres::project_point_on_hyperboloid(contour_barycenter, sphere1, sphere2);

			for(Contour::const_iterator jt=contour.begin();jt!=contour.end();++jt)
			{
				Contour::const_iterator jt_next=jt;
				++jt_next;
				if(jt_next==contour.end())
				{
					jt_next=contour.begin();
				}
				area+=triangle_area(contour_barycenter, jt->p, jt_next->p);
			}
		}
		return area;
	}

	static Float calculate_arc_length_from_contour(const UnsignedInt a_id, const Contour& contour) noexcept
	{
		Float arc_length=FLOATCONST(0.0);
		for(Contour::const_iterator jt=contour.begin();jt!=contour.end();++jt)
		{
			Contour::const_iterator jt_next=jt;
			++jt_next;
			if(jt_next==contour.end())
			{
				jt_next=contour.begin();
			}
			if(jt->right_id==a_id && jt_next->left_id==a_id)
			{
				arc_length+=distance_from_point_to_point(jt->p, jt_next->p);
			}
		}
		return arc_length;
	}
};

}

#endif /* VORONOTALT_SIMPLIFIED_AW_TESSELLATION_CONTACT_CONSTRUCTION_H_ */
