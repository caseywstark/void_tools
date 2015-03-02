#include <cmath>

#include <string>

// GP short for grid point
struct Point { double x; double y; double z; double value; };
struct Void { Point center; double radius; };

inline int
periodic_index(const int i, const int n)
{
    int r = i - n * (int)( (double)i / n );
    int r_neg = r < 0;
    r += r_neg * n;
    return r;
}

inline int
periodic_grid_index(const int i, const int j, const int k, const int n)
{
    int iw = periodic_index(i, n);
    int jw = periodic_index(j, n);
    int kw = periodic_index(k, n);
    int ii = (iw * n + jw) * n + kw;
    return ii;
}

inline double
grid_value(const int n, const double * const field,
    const int i, const int j, const int k)
{
    return field[(i * n + j) * n + k];
}

inline double
grid_value_periodic(const int n, const double * const field,
    const int i, const int j, const int k)
{
    int iw = periodic_index(i, n);
    int jw = periodic_index(j, n);
    int kw = periodic_index(k, n);
    return field[ (iw * n + jw) * n + kw ];
}

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

inline double
grid_spherical_average(const double l, const int n, const double * const field,
    const double cx, const double cy, const double cz, const double r)
{
    const double dx = l / n;

    double fw_sum = 0.0;
    double w_sum = 0.0;

    for (int i = floor(cx - r); i < ceil(cx + r); ++i) {
        int iw = periodic_index(i, n);
        double x = dx * (i + 0.5);
        for (int j = floor(cy - r); j < ceil(cy + r); ++j) {
            int jw = periodic_index(j, n);
            double y = dx * (j + 0.5);
            for (int k = floor(cz - r); k < ceil(cz + r); ++k) {
                int kw = periodic_index(k, n);
                double z = dx * (k + 0.5);

                double ri = periodic_distance(cx, cy, cz, x, y, z, l);
                // crappy approx to the fraction of the cell covered
                // by the sphere
                double dr = (r - ri) / dx + 0.5;
                int intersect = dr > 0.0;
                double wi = fmin(1.0, intersect * dr);

                double fi = grid_value(n, field, iw, jw, kw);
                fw_sum += fi * wi;
                w_sum += wi;
            }
        }
    }

    return fw_sum / w_sum;
}
