#pragma once

#include <Eigen/Core>

namespace Faunus {

using Point = Eigen::Vector3d;

/**
 * @brief Convert cartesian- to cylindrical-coordinates
 * @note Input (x,y,z), output \f$ (r,\theta, h) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in
 * [-\pi,\pi) \f$, and \f$ h\in (-\infty,\infty] \f$.
 */
Point xyz2rth(const Point&, const Point& origin = {0, 0, 0}, const Point& dir = {0, 0, 1},
              const Point& dir2 = {1, 0, 0});

/**
 * @brief Convert cartesian- to spherical-coordinates
 * @note Input (x,y,z), output \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in
 * [-\pi,\pi) \f$, and \f$ \phi\in [0,\pi] \f$.
 */
Point xyz2rtp(const Point&, const Point& origin = {0, 0, 0});

/**
 * @brief Convert spherical- to cartesian-coordinates
 * @param origin The origin to be added (optional)
 * @note Input \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [0,2\pi) \f$,
 * and \f$ \phi\in [0,\pi] \f$, and output (x,y,z).
 */
Point rtp2xyz(const Point& rtp, const Point& origin = {0, 0, 0});

} // namespace Faunus
