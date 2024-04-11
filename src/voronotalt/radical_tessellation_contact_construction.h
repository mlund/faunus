#ifndef VORONOTALT_RADICAL_TESSELLATION_CONTACT_CONSTRUCTION_H_
#define VORONOTALT_RADICAL_TESSELLATION_CONTACT_CONSTRUCTION_H_

#include <vector>

#include "basic_types_and_functions.h"

namespace voronotalt
{

class RadicalTessellationContactConstruction
{
public:
	struct ContourPoint
	{
		SimplePoint p;
		Float angle;
		UnsignedInt left_id;
		UnsignedInt right_id;
		int indicator;


		ContourPoint(const SimplePoint& p, const UnsignedInt left_id, const UnsignedInt right_id) :
			p(p),
			angle(FLOATCONST(0.0)),
			left_id(left_id),
			right_id(right_id),
			indicator(0)
		{
		}
	};

	typedef std::vector<ContourPoint> Contour;

	struct ContactDescriptor
	{
		Contour contour;
		SimpleSphere intersection_circle_sphere;
		SimplePoint intersection_circle_axis;
		SimplePoint contour_barycenter;
		Float sum_of_arc_angles;
		Float area;
		Float solid_angle_a;
		Float solid_angle_b;
		Float pyramid_volume_a;
		Float pyramid_volume_b;
		Float distance;
		UnsignedInt flags;
		UnsignedInt id_a;
		UnsignedInt id_b;

		ContactDescriptor() :
			sum_of_arc_angles(FLOATCONST(0.0)),
			area(FLOATCONST(0.0)),
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

		void clear()
		{
			id_a=0;
			id_b=0;
			contour.clear();
			sum_of_arc_angles=FLOATCONST(0.0);
			area=FLOATCONST(0.0);
			solid_angle_a=FLOATCONST(0.0);
			solid_angle_b=FLOATCONST(0.0);
			pyramid_volume_a=FLOATCONST(0.0);
			pyramid_volume_b=FLOATCONST(0.0);
			distance=FLOATCONST(0.0);
			flags=0;
		}
	};

	struct ContactDescriptorGraphics
	{
		std::vector<SimplePoint> outer_points;
		SimplePoint barycenter;

		ContactDescriptorGraphics()
		{
		}

		void clear()
		{
			outer_points.clear();
		}
	};

	static bool construct_contact_descriptor(
			const std::vector<SimpleSphere>& spheres,
			const std::vector<int>& spheres_exclusion_statuses,
			const UnsignedInt a_id,
			const UnsignedInt b_id,
			const std::vector<UnsignedInt>& a_neighbor_collisions,
			ContactDescriptor& result_contact_descriptor)
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
					bool contour_initialized=false;
					bool nostop=true;
					{
						for(UnsignedInt i=0;i<a_neighbor_collisions.size() && !discarded && nostop;i++)
						{
							const UnsignedInt neighbor_id=a_neighbor_collisions[i];
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
										const SimplePoint neighbor_ac_plane_center=center_of_intersection_circle_of_two_spheres(a, c);
										const SimplePoint neighbor_ac_plane_normal=unit_point(sub_of_points(c.p, a.p));
										const Float cos_val=dot_product(unit_point(sub_of_points(result_contact_descriptor.intersection_circle_sphere.p, a.p)), unit_point(sub_of_points(neighbor_ac_plane_center, a.p)));
										if(std::abs(cos_val)<FLOATCONST(1.0))
										{
											const Float l=std::abs(signed_distance_from_point_to_plane(neighbor_ac_plane_center, neighbor_ac_plane_normal, result_contact_descriptor.intersection_circle_sphere.p));
											const Float xl=l/std::sqrt(1-(cos_val*cos_val));
											if(xl>=result_contact_descriptor.intersection_circle_sphere.r)
											{
												if(halfspace_of_point(neighbor_ac_plane_center, neighbor_ac_plane_normal, result_contact_descriptor.intersection_circle_sphere.p)>=0)
												{
													discarded=true;
												}
											}
											else
											{
												if(!contour_initialized)
												{
													result_contact_descriptor.intersection_circle_axis=unit_point(sub_of_points(b.p, a.p));
													init_contour_from_base_and_axis(a_id, result_contact_descriptor.intersection_circle_sphere, result_contact_descriptor.intersection_circle_axis, result_contact_descriptor.contour);
													contour_initialized=true;
												}
												else
												{
													nostop=(test_if_contour_is_still_cuttable(a.p, neighbor_ac_plane_center, result_contact_descriptor.contour));
												}
												if(nostop)
												{
													mark_and_cut_contour(neighbor_ac_plane_center, neighbor_ac_plane_normal, neighbor_id, result_contact_descriptor.contour);
													if(contour_initialized && result_contact_descriptor.contour.empty())
													{
														discarded=true;
													}
												}
											}
										}
										else
										{
											if(halfspace_of_point(neighbor_ac_plane_center, neighbor_ac_plane_normal, result_contact_descriptor.intersection_circle_sphere.p)>0)
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
						if(!contour_initialized)
						{
							result_contact_descriptor.intersection_circle_axis=unit_point(sub_of_points(b.p, a.p));
							result_contact_descriptor.contour_barycenter=result_contact_descriptor.intersection_circle_sphere.p;
							result_contact_descriptor.sum_of_arc_angles=(PIVALUE*FLOATCONST(2.0));
							result_contact_descriptor.area=result_contact_descriptor.intersection_circle_sphere.r*result_contact_descriptor.intersection_circle_sphere.r*PIVALUE;
						}
						else
						{
							if(!result_contact_descriptor.contour.empty())
							{
								restrict_contour_to_circle(result_contact_descriptor.intersection_circle_sphere, result_contact_descriptor.intersection_circle_axis, a_id, result_contact_descriptor.contour, result_contact_descriptor.sum_of_arc_angles);
								if(!result_contact_descriptor.contour.empty())
								{
									result_contact_descriptor.area=calculate_contour_area(result_contact_descriptor.intersection_circle_sphere, result_contact_descriptor.contour, result_contact_descriptor.contour_barycenter);
								}
							}
						}

						if(result_contact_descriptor.area>FLOATCONST(0.0))
						{
							result_contact_descriptor.solid_angle_a=calculate_contour_solid_angle(a, b, result_contact_descriptor.intersection_circle_sphere, result_contact_descriptor.contour);
							result_contact_descriptor.solid_angle_b=calculate_contour_solid_angle(b, a, result_contact_descriptor.intersection_circle_sphere, result_contact_descriptor.contour);
							result_contact_descriptor.pyramid_volume_a=distance_from_point_to_point(a.p, result_contact_descriptor.intersection_circle_sphere.p)*result_contact_descriptor.area/FLOATCONST(3.0)*(result_contact_descriptor.solid_angle_a<FLOATCONST(0.0) ? FLOATCONST(-1.0) : FLOATCONST(1.0));
							result_contact_descriptor.pyramid_volume_b=distance_from_point_to_point(b.p, result_contact_descriptor.intersection_circle_sphere.p)*result_contact_descriptor.area/FLOATCONST(3.0)*(result_contact_descriptor.solid_angle_b<FLOATCONST(0.0) ? FLOATCONST(-1.0) : FLOATCONST(1.0));
							result_contact_descriptor.flags=(check_if_contour_is_central(result_contact_descriptor.intersection_circle_sphere.p, result_contact_descriptor.contour, result_contact_descriptor.contour_barycenter) ? 1 : 0);
						}

						result_contact_descriptor.distance=distance_from_point_to_point(a.p, b.p);
					}
				}
			}
		}
		return (result_contact_descriptor.area>FLOATCONST(0.0));
	}

	static bool construct_contact_descriptor_graphics(const ContactDescriptor& contact_descriptor, const Float length_step, ContactDescriptorGraphics& result_contact_descriptor_graphics)
	{
		result_contact_descriptor_graphics.clear();
		if(contact_descriptor.area>FLOATCONST(0.0))
		{
			const Float angle_step=std::max(std::min(length_step/contact_descriptor.intersection_circle_sphere.r, PIVALUE/FLOATCONST(3.0)), PIVALUE/FLOATCONST(36.0));
			if(contact_descriptor.contour.empty())
			{
				const SimplePoint first_point=point_and_number_product(any_normal_of_vector(contact_descriptor.intersection_circle_axis), contact_descriptor.intersection_circle_sphere.r);
				result_contact_descriptor_graphics.outer_points.reserve(static_cast<int>((PIVALUE*FLOATCONST(2.0))/angle_step)+2);
				result_contact_descriptor_graphics.outer_points.push_back(sum_of_points(contact_descriptor.intersection_circle_sphere.p, first_point));
				for(Float rotation_angle=angle_step;rotation_angle<(PIVALUE*FLOATCONST(2.0));rotation_angle+=angle_step)
				{
					result_contact_descriptor_graphics.outer_points.push_back(sum_of_points(contact_descriptor.intersection_circle_sphere.p, rotate_point_around_axis(contact_descriptor.intersection_circle_axis, rotation_angle, first_point)));
				}
				result_contact_descriptor_graphics.barycenter=contact_descriptor.intersection_circle_sphere.p;
			}
			else
			{
				if(contact_descriptor.sum_of_arc_angles>FLOATCONST(0.0))
				{
					result_contact_descriptor_graphics.outer_points.reserve(static_cast<UnsignedInt>(contact_descriptor.sum_of_arc_angles/angle_step)+contact_descriptor.contour.size()+4);
				}
				else
				{
					result_contact_descriptor_graphics.outer_points.reserve(contact_descriptor.contour.size());
				}
				for(UnsignedInt i=0;i<contact_descriptor.contour.size();i++)
				{
					const ContourPoint& pr=contact_descriptor.contour[i];
					result_contact_descriptor_graphics.outer_points.push_back(pr.p);
					if(pr.angle>FLOATCONST(0.0))
					{
						if(pr.angle>angle_step)
						{
							const SimplePoint first_v=sub_of_points(pr.p, contact_descriptor.intersection_circle_sphere.p);
							for(Float rotation_angle=angle_step;rotation_angle<pr.angle;rotation_angle+=angle_step)
							{
								result_contact_descriptor_graphics.outer_points.push_back(sum_of_points(contact_descriptor.intersection_circle_sphere.p, rotate_point_around_axis(contact_descriptor.intersection_circle_axis, rotation_angle, first_v)));
							}
						}
					}
				}
				result_contact_descriptor_graphics.barycenter.x=FLOATCONST(0.0);
				result_contact_descriptor_graphics.barycenter.y=FLOATCONST(0.0);
				result_contact_descriptor_graphics.barycenter.z=FLOATCONST(0.0);
				if(!result_contact_descriptor_graphics.outer_points.empty())
				{
					for(UnsignedInt i=0;i<result_contact_descriptor_graphics.outer_points.size();i++)
					{
						set_close_to_equal(result_contact_descriptor_graphics.outer_points[i].x, FLOATCONST(0.0));
						set_close_to_equal(result_contact_descriptor_graphics.outer_points[i].y, FLOATCONST(0.0));
						set_close_to_equal(result_contact_descriptor_graphics.outer_points[i].z, FLOATCONST(0.0));
						result_contact_descriptor_graphics.barycenter.x+=result_contact_descriptor_graphics.outer_points[i].x;
						result_contact_descriptor_graphics.barycenter.y+=result_contact_descriptor_graphics.outer_points[i].y;
						result_contact_descriptor_graphics.barycenter.z+=result_contact_descriptor_graphics.outer_points[i].z;
					}
					result_contact_descriptor_graphics.barycenter.x/=static_cast<Float>(result_contact_descriptor_graphics.outer_points.size());
					result_contact_descriptor_graphics.barycenter.y/=static_cast<Float>(result_contact_descriptor_graphics.outer_points.size());
					result_contact_descriptor_graphics.barycenter.z/=static_cast<Float>(result_contact_descriptor_graphics.outer_points.size());
				}
			}
		}
		return (!result_contact_descriptor_graphics.outer_points.empty());
	}

private:
	static void init_contour_from_base_and_axis(
			const UnsignedInt a_id,
			const SimpleSphere& base,
			const SimplePoint& axis,
			Contour& result)
	{
		const SimplePoint first_point=point_and_number_product(any_normal_of_vector(axis), base.r*FLOATCONST(1.19));
		const Float angle_step=PIVALUE/FLOATCONST(3.0);
		result.push_back(ContourPoint(sum_of_points(base.p, first_point), a_id, a_id));
		for(Float rotation_angle=angle_step;rotation_angle<(PIVALUE*FLOATCONST(2.0));rotation_angle+=angle_step)
		{
			result.push_back(ContourPoint(sum_of_points(base.p, rotate_point_around_axis(axis, rotation_angle, first_point)), a_id, a_id));
		}
	}

	static bool mark_and_cut_contour(
			const SimplePoint& ac_plane_center,
			const SimplePoint& ac_plane_normal,
			const UnsignedInt c_id,
			Contour& contour)
	{
		const UnsignedInt outsiders_count=mark_contour(ac_plane_center, ac_plane_normal, c_id, contour);
		if(outsiders_count>0)
		{
			if(outsiders_count<contour.size())
			{
				cut_contour(ac_plane_center, ac_plane_normal, c_id, contour);
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
			const SimplePoint& ac_plane_center,
			const SimplePoint& ac_plane_normal,
			const UnsignedInt c_id,
			Contour& contour)
	{
		UnsignedInt outsiders_count=0;
		for(Contour::iterator it=contour.begin();it!=contour.end();++it)
		{
			if(halfspace_of_point(ac_plane_center, ac_plane_normal, it->p)>=0)
			{
				it->left_id=c_id;
				it->right_id=c_id;
				outsiders_count++;
			}
		}
		return outsiders_count;
	}

	static void cut_contour(
			const SimplePoint& ac_plane_center,
			const SimplePoint& ac_plane_normal,
			const UnsignedInt c_id,
			Contour& contour)
	{
		if(contour.size()<3)
		{
			return;
		}

		UnsignedInt i_start=0;
		while(i_start<contour.size() && !(contour[i_start].left_id==c_id && contour[i_start].right_id==c_id))
		{
			i_start++;
		}

		if(i_start>=contour.size())
		{
			return;
		}

		UnsignedInt i_end=(contour.size()-1);
		while(i_end>0 && !(contour[i_end].left_id==c_id && contour[i_end].right_id==c_id))
		{
			i_end--;
		}

		if(i_start==0 && i_end==(contour.size()-1))
		{
			i_end=0;
			while((i_end+1)<contour.size() && contour[i_end+1].left_id==c_id && contour[i_end+1].right_id==c_id)
			{
				i_end++;
			}

			i_start=(contour.size()-1);
			while(i_start>0 && contour[i_start-1].left_id==c_id && contour[i_start-1].right_id==c_id)
			{
				i_start--;
			}
		}

		if(i_start==i_end)
		{
			contour.insert(contour.begin()+i_start, contour[i_start]);
			i_end=i_start+1;
		}
		else if(i_start<i_end)
		{
			if(i_start+1<i_end)
			{
				contour.erase(contour.begin()+i_start+1, contour.begin()+i_end);
			}
			i_end=i_start+1;
		}
		else if(i_start>i_end)
		{
			if(i_start+1<contour.size())
			{
				contour.erase(contour.begin()+i_start+1, contour.end());
			}
			if(i_end>0)
			{
				contour.erase(contour.begin(), contour.begin()+i_end);

			}
			i_start=contour.size()-1;
			i_end=0;
		}

		{
			const UnsignedInt i_left=((i_start)>0 ? (i_start-1) : (contour.size()-1));
			contour[i_start]=ContourPoint(intersection_of_plane_and_segment(ac_plane_center, ac_plane_normal, contour[i_start].p, contour[i_left].p), contour[i_left].right_id, contour[i_start].left_id);
		}

		{
			const UnsignedInt i_right=((i_end+1)<contour.size() ? (i_end+1) : 0);
			contour[i_end]=ContourPoint(intersection_of_plane_and_segment(ac_plane_center, ac_plane_normal, contour[i_end].p, contour[i_right].p), contour[i_end].right_id, contour[i_right].left_id);
		}

		if(!greater(squared_distance_from_point_to_point(contour[i_start].p, contour[i_end].p), FLOATCONST(0.0)))
		{
			contour[i_start].right_id=contour[i_end].right_id;
			contour.erase(contour.begin()+i_end);
		}
	}

	static bool test_if_contour_is_still_cuttable(const SimplePoint& a_center, const SimplePoint& closest_possible_cut_point, const Contour& contour)
	{
		bool cuttable=false;
		const Float dist_threshold=squared_distance_from_point_to_point(a_center, closest_possible_cut_point);
		for(UnsignedInt i=0;!cuttable && i<contour.size();i++)
		{
			cuttable=squared_distance_from_point_to_point(a_center, contour[i].p)>=dist_threshold;
		}
		return cuttable;
	}

	static bool restrict_contour_to_circle(
			const SimpleSphere& ic_sphere,
			const SimplePoint& ic_axis,
			const UnsignedInt a_id,
			Contour& contour,
			Float& sum_of_arc_angles)
	{
		UnsignedInt outsiders_count=0;
		for(UnsignedInt i=0;i<contour.size();i++)
		{
			if(squared_distance_from_point_to_point(contour[i].p, ic_sphere.p)<=(ic_sphere.r*ic_sphere.r))
			{
				contour[i].indicator=0;
			}
			else
			{
				contour[i].indicator=1;
				outsiders_count++;
			}
		}

		if(outsiders_count>0)
		{
			UnsignedInt insertions_count=0;
			{
				UnsignedInt i=0;
				while(i<contour.size())
				{
					ContourPoint& pr1=contour[i];
					ContourPoint& pr2=contour[(i+1)<contour.size() ? (i+1) : 0];
					if(pr1.indicator==1 || pr2.indicator==1)
					{
						if(pr1.indicator==1 && pr2.indicator==1)
						{
							SimplePoint mp;
							if(project_point_inside_line(ic_sphere.p, pr1.p, pr2.p, mp))
							{
								if(squared_distance_from_point_to_point(mp, ic_sphere.p)<=(ic_sphere.r*ic_sphere.r))
								{
									SimplePoint ip1;
									SimplePoint ip2;
									if(intersect_segment_with_circle(ic_sphere, mp, pr1.p, ip1) && intersect_segment_with_circle(ic_sphere, mp, pr2.p, ip2))
									{
										const UnsignedInt pr2_left_id=pr2.left_id;
										contour.insert(((i+1)<contour.size() ? (contour.begin()+(i+1)) : contour.end()), 2, ContourPoint(ip1, a_id, pr1.right_id));
										contour[i+2]=ContourPoint(ip2, pr2_left_id, a_id);
										insertions_count+=2;
										i+=2;
									}
								}
							}
						}
						else
						{
							if(pr1.indicator==1 && pr2.indicator!=1)
							{
								SimplePoint ip;
								if(intersect_segment_with_circle(ic_sphere, pr2.p, pr1.p, ip))
								{
									contour.insert(((i+1)<contour.size() ? (contour.begin()+(i+1)) : contour.end()), ContourPoint(ip, a_id, pr1.right_id));
									insertions_count++;
									i++;
								}
								else
								{
									pr2.left_id=a_id;
									pr2.right_id=pr1.right_id;
								}
							}
							else if(pr1.indicator!=1 && pr2.indicator==1)
							{
								SimplePoint ip;
								if(intersect_segment_with_circle(ic_sphere, pr1.p, pr2.p, ip))
								{
									contour.insert(((i+1)<contour.size() ? (contour.begin()+(i+1)) : contour.end()), ContourPoint(ip, pr2.left_id, a_id));
									insertions_count++;
									i++;
								}
								else
								{
									pr1.left_id=pr2.left_id;
									pr1.right_id=a_id;
								}
							}
						}
					}
					i++;
				}
			}
			if(insertions_count==0)
			{
				contour.clear();
			}
			else if(!contour.empty())
			{
				{
					UnsignedInt i=(contour.size()-1);
					while(i<contour.size())
					{
						if(contour[i].indicator==1)
						{
							contour.erase(contour.begin()+i);
						}
						i=(i>0 ? (i-1) : contour.size());
					}
				}
				if(contour.size()<2)
				{
					contour.clear();
				}
				{
					UnsignedInt i=0;
					while(i<contour.size())
					{
						ContourPoint& pr1=contour[i];
						ContourPoint& pr2=contour[(i+1)<contour.size() ? (i+1) : 0];
						if(pr1.right_id==a_id && pr2.left_id==a_id)
						{
							pr1.angle=directed_angle(ic_sphere.p, pr1.p, pr2.p, sum_of_points(ic_sphere.p, ic_axis));
							sum_of_arc_angles+=pr1.angle;
							pr1.indicator=2;
							pr2.indicator=3;
						}
						i++;
					}
				}
			}
		}

		return (outsiders_count>0);
	}

	static Float calculate_contour_area(const SimpleSphere& ic_sphere, const Contour& contour, SimplePoint& contour_barycenter)
	{
		Float area=FLOATCONST(0.0);

		contour_barycenter.x=FLOATCONST(0.0);
		contour_barycenter.y=FLOATCONST(0.0);
		contour_barycenter.z=FLOATCONST(0.0);
		for(UnsignedInt i=0;i<contour.size();i++)
		{
			contour_barycenter.x+=contour[i].p.x;
			contour_barycenter.y+=contour[i].p.y;
			contour_barycenter.z+=contour[i].p.z;
		}
		contour_barycenter.x/=static_cast<Float>(contour.size());
		contour_barycenter.y/=static_cast<Float>(contour.size());
		contour_barycenter.z/=static_cast<Float>(contour.size());

		for(UnsignedInt i=0;i<contour.size();i++)
		{
			const ContourPoint& pr1=contour[i];
			const ContourPoint& pr2=contour[(i+1)<contour.size() ? (i+1) : 0];
			area+=triangle_area(contour_barycenter, pr1.p, pr2.p);
			if(pr1.angle>FLOATCONST(0.0))
			{
				area+=ic_sphere.r*ic_sphere.r*(pr1.angle-std::sin(pr1.angle))*0.5;
			}
		}

		return area;
	}

	static Float calculate_contour_solid_angle(const SimpleSphere& a, const SimpleSphere& b, const SimpleSphere& ic_sphere, const Contour& contour)
	{
		Float turn_angle=FLOATCONST(0.0);

		if(!contour.empty())
		{
			for(UnsignedInt i=0;i<contour.size();i++)
			{
				const ContourPoint& pr0=contour[(i>0) ? (i-1) : (contour.size()-1)];
				const ContourPoint& pr1=contour[i];
				const ContourPoint& pr2=contour[(i+1)<contour.size() ? (i+1) : 0];

				if(pr0.angle>FLOATCONST(0.0))
				{
					SimplePoint d=cross_product(sub_of_points(b.p, a.p), sub_of_points(pr1.p, ic_sphere.p));
					if((pr0.angle<PIVALUE && dot_product(d, sub_of_points(pr0.p, pr1.p))<FLOATCONST(0.0)) || (pr0.angle>PIVALUE && dot_product(d, sub_of_points(pr0.p, pr1.p))>FLOATCONST(0.0)))
					{
						d=point_and_number_product(d, FLOATCONST(-1.0));
					}
					turn_angle+=(PIVALUE-min_dihedral_angle(a.p, pr1.p, sum_of_points(pr1.p, d), pr2.p));
				}
				else if(pr1.angle>FLOATCONST(0.0))
				{
					SimplePoint d=cross_product(sub_of_points(b.p, a.p), sub_of_points(pr1.p, ic_sphere.p));
					if((pr1.angle<PIVALUE && dot_product(d, sub_of_points(pr2.p, pr1.p))<FLOATCONST(0.0)) || (pr1.angle>PIVALUE && dot_product(d, sub_of_points(pr2.p, pr1.p))>FLOATCONST(0.0)))
					{
						d=point_and_number_product(d, FLOATCONST(-1.0));
					}
					turn_angle+=(PIVALUE-min_dihedral_angle(a.p, pr1.p, pr0.p, sum_of_points(pr1.p, d)));
					turn_angle+=pr1.angle*(distance_from_point_to_point(a.p, ic_sphere.p)/a.r);
				}
				else
				{
					turn_angle+=(PIVALUE-min_dihedral_angle(a.p, pr1.p, pr0.p, pr2.p));
				}
			}
		}
		else
		{
			turn_angle=(PIVALUE*FLOATCONST(2.0))*(distance_from_point_to_point(a.p, ic_sphere.p)/a.r);
		}

		Float solid_angle=((PIVALUE*FLOATCONST(2.0))-turn_angle);

		if(dot_product(sub_of_points(ic_sphere.p, a.p), sub_of_points(ic_sphere.p, b.p))>FLOATCONST(0.0) && squared_distance_from_point_to_point(ic_sphere.p, a.p)<squared_distance_from_point_to_point(ic_sphere.p, b.p))
		{
			solid_angle=(FLOATCONST(0.0)-solid_angle);
		}

		return solid_angle;
	}

	static bool check_if_contour_is_central(const SimplePoint& center, const Contour& contour, const SimplePoint& contour_barycenter)
	{
		bool central=true;
		for(UnsignedInt i=0;i<contour.size() && central;i++)
		{
			const UnsignedInt j=((i+1)<contour.size() ? (i+1) : 0);
			const SimplePoint sidepoint=point_and_number_product(sum_of_points(contour[i].p, contour[j].p), FLOATCONST(0.5));
			if(dot_product(sub_of_points(contour_barycenter, sidepoint), sub_of_points(center, sidepoint))<FLOATCONST(0.0))
			{
				central=false;
			}
		}
		return central;
	}
};

}

#endif /* VORONOTALT_RADICAL_TESSELLATION_CONTACT_CONSTRUCTION_H_ */
