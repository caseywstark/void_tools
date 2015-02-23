#include <cmath>

#include <string>

struct Point { int x; int y; int z; double value;};
struct Void { Point center; double radius; };

inline std::vector<double>
read_grid_data(const std::string grid_path, const size_t nn)
{
    // allocate array
    std::vector<double> d(nn);

    // open file and read in directly
    FILE *f = fopen(grid_path.c_str(), "r");
    fread(&d[0], sizeof(double), nn, f);
    fclose(f);

    return d;
}

inline int
wrap_index(const int i, const int n)
{
    int r = i - n * (int)( (double)i / n );
    int r_neg = r < 0;
    r += r_neg * n;
    return r;
}

inline int
periodic_grid_index(const int i, const int j, const int k, const int n)
{
    int iw = wrap_index(i, n);
    int jw = wrap_index(j, n);
    int kw = wrap_index(k, n);
    int ii = (iw * n + jw) * n + kw;
    return ii;
}

inline double
wrap_distance(const double x0, const double x1, const double l)
{
    double dx = x1 - x0;
    if (dx < -l/2) { dx += l; }
    else if (dx > l/2) { dx -= l; }
    return dx;
}

inline double
periodic_distance(const Point pi, const Point pj, const double l)
{
    double dx = wrap_distance(pi.x, pj.x, l);
    double dy = wrap_distance(pi.y, pj.y, l);
    double dz = wrap_distance(pi.z, pj.z, l);
    return sqrt(dx*dx + dy*dy + dz*dz);
}

inline double
grid_value(const int ng, const std::vector<double> &field,
    const int i, const int j, const int k)
{
    return field[(i * ng + j) * ng + k];
}

inline double
grid_value_periodic(const int ng, const std::vector<double> &field,
    const int i, const int j, const int k)
{
    int iw = wrap_index(i, ng);
    int jw = wrap_index(j, ng);
    int kw = wrap_index(k, ng);
    return field[ (iw * ng + jw) * ng + kw ];
}

template<int ng>
inline double
grid_spherical_average(const std::vector<double> &field,
    const Point center, const double r)
{
    double f_sum = 0.0;
    int f_count = 0;

    Point pi;
    for (int i = floor(center.x - r); i < ceil(center.x + r); ++i) {
        int iw = wrap_index(i, ng);
        pi.x = i;
        for (int j = floor(center.y - r); j < ceil(center.y + r); ++j) {
            int jw = wrap_index(j, ng);
            pi.y = j;
            for (int k = floor(center.z - r); k < ceil(center.z + r); ++k) {
                int kw = wrap_index(k, ng);
                pi.z = k;

                double ri = periodic_distance(center, pi, (double)ng);

                if (ri < r) {
                    f_sum += grid_value(ng, field, iw, jw, kw);
                    f_count += 1;
                }
            }
        }
    }

    return f_sum / f_count;
}
