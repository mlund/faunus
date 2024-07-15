#ifndef VORONOTALT_BASIC_TYPES_AND_FUNCTIONS_H_
#define VORONOTALT_BASIC_TYPES_AND_FUNCTIONS_H_

#include <cmath>

//#define VORONOTALT_FP32
#define VORONOTALT_UI32

#ifdef VORONOTALT_FP32
#define FLOATCONST( v ) v ## f
#define FLOATEPSILON 0.000001f
#define PIVALUE 3.14159265358979323846f
#else
#define FLOATCONST( v ) v
#define FLOATEPSILON 0.0000000001
#define PIVALUE 3.14159265358979323846
#endif

namespace voronotalt
{

#ifdef VORONOTALT_FP32
typedef float Float;
#else
typedef double Float;
#endif

#ifdef VORONOTALT_UI32
typedef unsigned int UnsignedInt;
#else
typedef std::size_t UnsignedInt;
#endif

struct SimplePoint
{
	Float x;
	Float y;
	Float z;

	SimplePoint() : x(FLOATCONST(0.0)), y(FLOATCONST(0.0)), z(FLOATCONST(0.0))
	{
	}

	SimplePoint(const Float x, const Float y, const Float z) : x(x), y(y), z(z)
	{
	}
};

struct SimpleSphere
{
	SimplePoint p;
	Float r;

	SimpleSphere() : r(FLOATCONST(0.0))
	{
	}

	SimpleSphere(const SimplePoint& p, const Float r) : p(p), r(r)
	{
	}
};

struct SimpleQuaternion
{
	Float a;
	Float b;
	Float c;
	Float d;

	SimpleQuaternion(const Float a, const Float b, const Float c, const Float d) : a(a), b(b), c(c), d(d)
	{
	}

	SimpleQuaternion(const Float a, const SimplePoint& p) : a(a), b(p.x), c(p.y), d(p.z)
	{
	}
};

inline bool equal(const Float a, const Float b, const Float e)
{
	return (((a-b)<=e) && ((b-a)<=e));
}

inline bool equal(const Float a, const Float b)
{
	return equal(a, b, FLOATEPSILON);
}

inline bool less(const Float a, const Float b)
{
	return ((a+FLOATEPSILON)<b);
}

inline bool greater(const Float a, const Float b)
{
	return ((a-FLOATEPSILON)>b);
}

inline bool less_or_equal(const Float a, const Float b)
{
	return (less(a, b) || equal(a, b));
}

inline bool greater_or_equal(const Float a, const Float b)
{
	return (greater(a, b) || equal(a, b));
}

inline void set_close_to_equal(Float& a, const Float b)
{
	if(equal(a, b))
	{
		a=b;
	}
}

inline Float squared_distance_from_point_to_point(const SimplePoint& a, const SimplePoint& b)
{
	const Float dx=(a.x-b.x);
	const Float dy=(a.y-b.y);
	const Float dz=(a.z-b.z);
	return (dx*dx+dy*dy+dz*dz);
}

inline Float distance_from_point_to_point(const SimplePoint& a, const SimplePoint& b)
{
	return sqrt(squared_distance_from_point_to_point(a, b));
}

inline Float squared_point_module(const SimplePoint& a)
{
	return (a.x*a.x+a.y*a.y+a.z*a.z);
}

inline Float point_module(const SimplePoint& a)
{
	return sqrt(squared_point_module(a));
}

inline Float dot_product(const SimplePoint& a, const SimplePoint& b)
{
	return (a.x*b.x+a.y*b.y+a.z*b.z);
}

inline SimplePoint cross_product(const SimplePoint& a, const SimplePoint& b)
{
	return SimplePoint(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

inline SimplePoint point_and_number_product(const SimplePoint& a, const Float k)
{
	return SimplePoint(a.x*k, a.y*k, a.z*k);
}

inline SimplePoint unit_point(const SimplePoint& a)
{
	return ((equal(squared_point_module(a), FLOATCONST(1.0))) ? a : point_and_number_product(a, FLOATCONST(1.0)/point_module(a)));
}

inline SimplePoint sum_of_points(const SimplePoint& a, const SimplePoint& b)
{
	return SimplePoint(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline SimplePoint sub_of_points(const SimplePoint& a, const SimplePoint& b)
{
	return SimplePoint(a.x-b.x, a.y-b.y, a.z-b.z);
}

inline Float signed_distance_from_point_to_plane(const SimplePoint& plane_point, const SimplePoint& plane_normal, const SimplePoint& x)
{
	return dot_product(unit_point(plane_normal), sub_of_points(x, plane_point));
}

inline int halfspace_of_point(const SimplePoint& plane_point, const SimplePoint& plane_normal, const SimplePoint& x)
{
	const Float sd=signed_distance_from_point_to_plane(plane_point, plane_normal, x);
	if(greater(sd, FLOATCONST(0.0)))
	{
		return 1;
	}
	else if(less(sd, FLOATCONST(0.0)))
	{
		return -1;
	}
	return 0;
}

inline SimplePoint intersection_of_plane_and_segment(const SimplePoint& plane_point, const SimplePoint& plane_normal, const SimplePoint& a, const SimplePoint& b)
{
	const Float da=signed_distance_from_point_to_plane(plane_point, plane_normal, a);
	const Float db=signed_distance_from_point_to_plane(plane_point, plane_normal, b);
	if((da-db)==0)
	{
		return a;
	}
	else
	{
		const Float t=da/(da-db);
		return sum_of_points(a, point_and_number_product(sub_of_points(b, a), t));
	}
}

inline Float triangle_area(const SimplePoint& a, const SimplePoint& b, const SimplePoint& c)
{
	return (point_module(cross_product(sub_of_points(b, a), sub_of_points(c, a)))/FLOATCONST(2.0));
}

inline Float min_angle(const SimplePoint& o, const SimplePoint& a, const SimplePoint& b)
{
	Float cos_val=dot_product(unit_point(sub_of_points(a, o)), unit_point(sub_of_points(b, o)));
	if(cos_val<FLOATCONST(-1.0))
	{
		cos_val=FLOATCONST(-1.0);
	}
	else if(cos_val>FLOATCONST(1.0))
	{
		cos_val=FLOATCONST(1.0);
	}
	return std::acos(cos_val);
}

inline Float directed_angle(const SimplePoint& o, const SimplePoint& a, const SimplePoint& b, const SimplePoint& c)
{
	const Float angle=min_angle(o, a, b);
	const SimplePoint n=cross_product(unit_point(sub_of_points(a, o)), unit_point(sub_of_points(b, o)));
	if(dot_product(sub_of_points(c, o), n)>=0)
	{
		return angle;
	}
	else
	{
		return (PIVALUE*FLOATCONST(2.0)-angle);
	}
}

SimplePoint any_normal_of_vector(const SimplePoint& a)
{
	SimplePoint b=a;
	if(!equal(b.x, FLOATCONST(0.0)) && (!equal(b.y, FLOATCONST(0.0)) || !equal(b.z, FLOATCONST(0.0))))
	{
		b.x=FLOATCONST(0.0)-b.x;
		return unit_point(cross_product(a, b));
	}
	else if(!equal(b.y, FLOATCONST(0.0)) && (!equal(b.x, FLOATCONST(0.0)) || !equal(b.z, FLOATCONST(0.0))))
	{
		b.y=FLOATCONST(0.0)-b.y;
		return unit_point(cross_product(a, b));
	}
	else if(!equal(b.x, FLOATCONST(0.0)))
	{
		return SimplePoint(FLOATCONST(0.0), FLOATCONST(1.0), FLOATCONST(0.0));
	}
	else
	{
		return SimplePoint(FLOATCONST(1.0), FLOATCONST(0.0), FLOATCONST(0.0));
	}
}

inline bool sphere_intersects_sphere(const SimpleSphere& a, const SimpleSphere& b)
{
	return less(squared_distance_from_point_to_point(a.p, b.p), (a.r+b.r)*(a.r+b.r));
}

inline bool sphere_equals_sphere(const SimpleSphere& a, const SimpleSphere& b)
{
	return (equal(a.r, b.r) && equal(a.p.x, b.p.x) && equal(a.p.y, b.p.y) && equal(a.p.z, b.p.z));
}

inline bool sphere_contains_sphere(const SimpleSphere& a, const SimpleSphere& b)
{
	return (greater_or_equal(a.r, b.r) && less_or_equal(squared_distance_from_point_to_point(a.p, b.p), (a.r-b.r)*(a.r-b.r)));
}

inline Float distance_to_center_of_intersection_circle_of_two_spheres(const SimpleSphere& a, const SimpleSphere& b)
{
	const Float cm=distance_from_point_to_point(a.p, b.p);
	const Float cos_g=(a.r*a.r+cm*cm-b.r*b.r)/(2*a.r*cm);
	return (a.r*cos_g);
}

inline SimplePoint center_of_intersection_circle_of_two_spheres(const SimpleSphere& a, const SimpleSphere& b)
{
	const SimplePoint cv=sub_of_points(b.p, a.p);
	const Float cm=point_module(cv);
	const Float cos_g=(a.r*a.r+cm*cm-b.r*b.r)/(2*a.r*cm);
	return sum_of_points(a.p, point_and_number_product(cv, a.r*cos_g/cm));
}

inline SimpleSphere intersection_circle_of_two_spheres(const SimpleSphere& a, const SimpleSphere& b)
{
	const SimplePoint cv=sub_of_points(b.p, a.p);
	const Float cm=point_module(cv);
	const Float cos_g=(a.r*a.r+cm*cm-b.r*b.r)/(2*a.r*cm);
	const Float sin_g=std::sqrt(1-cos_g*cos_g);
	return SimpleSphere(sum_of_points(a.p, point_and_number_product(cv, a.r*cos_g/cm)), a.r*sin_g);
}

inline bool project_point_inside_line(const SimplePoint& o, const SimplePoint& a, const SimplePoint& b, SimplePoint& result)
{
	const SimplePoint v=unit_point(sub_of_points(b, a));
	const Float l=dot_product(v, sub_of_points(o, a));
	if(l>FLOATCONST(0.0) && (l*l)<=squared_distance_from_point_to_point(a, b))
	{
		result=sum_of_points(a, point_and_number_product(v, l));
		return true;
	}
	return false;
}

inline bool intersect_segment_with_circle(const SimpleSphere& circle, const SimplePoint& p_in, const SimplePoint& p_out, SimplePoint& result)
{
	const Float dist_in_to_out=distance_from_point_to_point(p_in, p_out);
	if(dist_in_to_out>FLOATCONST(0.0))
	{
		const SimplePoint v=point_and_number_product(sub_of_points(p_in, p_out), FLOATCONST(1.0)/dist_in_to_out);
		const SimplePoint u=sub_of_points(circle.p, p_out);
		const SimplePoint s=sum_of_points(p_out, point_and_number_product(v, dot_product(v, u)));
		const Float ll=(circle.r*circle.r)-squared_distance_from_point_to_point(circle.p, s);
		if(ll>=FLOATCONST(0.0))
		{
			result=sum_of_points(s, point_and_number_product(v, FLOATCONST(0.0)-std::sqrt(ll)));
			return true;
		}
	}
	return false;
}

inline Float min_dihedral_angle(const SimplePoint& o, const SimplePoint& a, const SimplePoint& b1, const SimplePoint& b2)
{
	const SimplePoint oa=unit_point(sub_of_points(a, o));
	const SimplePoint d1=sub_of_points(b1, sum_of_points(o, point_and_number_product(oa, dot_product(oa, sub_of_points(b1, o)))));
	const SimplePoint d2=sub_of_points(b2, sum_of_points(o, point_and_number_product(oa, dot_product(oa, sub_of_points(b2, o)))));
	const Float cos_val=dot_product(unit_point(d1), unit_point(d2));
	return std::acos(std::max(FLOATCONST(-1.0), std::min(cos_val, FLOATCONST(1.0))));
}

inline SimpleQuaternion quaternion_product(const SimpleQuaternion& q1, const SimpleQuaternion& q2)
{
	return SimpleQuaternion(
			q1.a*q2.a - q1.b*q2.b - q1.c*q2.c - q1.d*q2.d,
			q1.a*q2.b + q1.b*q2.a + q1.c*q2.d - q1.d*q2.c,
			q1.a*q2.c - q1.b*q2.d + q1.c*q2.a + q1.d*q2.b,
			q1.a*q2.d + q1.b*q2.c - q1.c*q2.b + q1.d*q2.a);
}

inline SimpleQuaternion inverted_quaternion(const SimpleQuaternion& q)
{
	return SimpleQuaternion(q.a, FLOATCONST(0.0)-q.b, FLOATCONST(0.0)-q.c, FLOATCONST(0.0)-q.d);
}

inline SimplePoint rotate_point_around_axis(const SimplePoint& axis, const Float angle, const SimplePoint& p)
{
	if(squared_point_module(axis)>0)
	{
		const Float radians_angle_half=(angle*FLOATCONST(0.5));
		const SimpleQuaternion q1=SimpleQuaternion(std::cos(radians_angle_half), point_and_number_product(unit_point(axis), std::sin(radians_angle_half)));
		const SimpleQuaternion q2=SimpleQuaternion(FLOATCONST(0.0), p);
		const SimpleQuaternion q3=quaternion_product(quaternion_product(q1, q2), inverted_quaternion(q1));
		return SimplePoint(q3.b, q3.c, q3.d);
	}
	else
	{
		return p;
	}
}

}

#endif /* VORONOTALT_BASIC_TYPES_AND_FUNCTIONS_H_ */
