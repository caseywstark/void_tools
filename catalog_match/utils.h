

#ifndef __util_h__
#define __util_h__

#include <cmath>

struct Point { double x; double y; double z; double value; };
struct Void { Point center; double radius; };

inline double
periodic_distance(const double x0, const double x1, const double l)
{
    double dx = x1 - x0;
    if (dx < -l/2) { dx += l; }
    else if (dx > l/2) { dx -= l; }
    return dx;
}

inline double
periodic_distance(const double x0, const double y0, const double z0,
    const double x1, const double y1, const double z1, const double l)
{
    double dx = periodic_distance(x0, x1, l);
    double dy = periodic_distance(y0, y1, l);
    double dz = periodic_distance(z0, z1, l);
    return sqrt(dx*dx + dy*dy + dz*dz);
}

inline double
periodic_distance(const Point pi, const Point pj, const double l)
{
    return periodic_distance(pi.x, pi.y, pi.z, pj.x, pj.y, pj.z, l);
}

#endif
