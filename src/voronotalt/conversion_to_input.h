#ifndef VORONOTALT_CONVERSION_TO_INPUT_H_
#define VORONOTALT_CONVERSION_TO_INPUT_H_

#include <vector>

#include "basic_types_and_functions.h"

namespace voronotalt
{

template<class Ball>
inline SimpleSphere get_sphere_from_ball(const Ball& ball, const Float probe) noexcept
{
	return SimpleSphere(SimplePoint(static_cast<Float>(ball.x), static_cast<Float>(ball.y), static_cast<Float>(ball.z)), static_cast<Float>(ball.r)+probe);
}

template<class Ball>
inline void fill_sphere_from_ball(const Ball& ball, const Float probe, SimpleSphere& sphere) noexcept
{
	sphere.p.x=static_cast<Float>(ball.x);
	sphere.p.y=static_cast<Float>(ball.y);
	sphere.p.z=static_cast<Float>(ball.z);
	sphere.r=static_cast<Float>(ball.r)+probe;
}

template<class BallsContainer>
inline std::vector<SimpleSphere> get_spheres_from_balls(const BallsContainer& balls, const Float probe) noexcept
{
	std::vector<SimpleSphere> result;
	result.reserve(balls.size());
	for(typename BallsContainer::const_iterator it=balls.begin();it!=balls.end();++it)
	{
		result.push_back(get_sphere_from_ball(*it, probe));
	}
	return result;
}

template<class BallsContainer>
inline void fill_spheres_from_balls(const BallsContainer& balls, const Float probe, std::vector<SimpleSphere>& spheres) noexcept
{
	spheres.resize(balls.size());
	std::size_t i=0;
	for(typename BallsContainer::const_iterator it=balls.begin();it!=balls.end();++it)
	{
		fill_sphere_from_ball(*it, probe, spheres[i]);
		i++;
	}
}

template<class Point>
inline SimplePoint get_simple_point_from_point(const Point& point) noexcept
{
	return SimplePoint(static_cast<Float>(point.x), static_cast<Float>(point.y), static_cast<Float>(point.z));
}

template<class PointsContainer>
inline std::vector<SimplePoint> get_simple_points_from_points(const PointsContainer& points) noexcept
{
	std::vector<SimplePoint> result;
	result.reserve(points.size());
	for(typename PointsContainer::const_iterator it=points.begin();it!=points.end();++it)
	{
		result.push_back(get_simple_point_from_point(*it));
	}
	return result;
}

}

#endif /* VORONOTALT_CONVERSION_TO_INPUT_H_ */
