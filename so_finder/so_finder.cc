/*

Simple void finder.
2015, Casey W. Stark.

- read grid data
- find points above/below threshold
- grow spheres until they hit the target average.
- remove overlapping voids

*/

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include "utils.h"

using namespace std;

// grid resolution
// define at compile time for unrolling loops.
const int n = 256;
const int nn = n * n * n;

const double l = 256.0;
const double dx = l / n;

inline bool
compare_by_mode(const bool comp_lt, const double value, const double ref)
{
    if (comp_lt) { return value < ref; }
    return value >= ref;
}

inline bool
void_radius_ge(const Void &a, const Void &b)
{
    return b.radius < a.radius;
}

inline bool
void_value_lt(const Void &a, const Void &b)
{
    return a.center.value < b.center.value;
}

std::vector<Void>
find_sos(const double * const field,
    const double thresh_value, const double avg_target,
    const bool overlap_mode_radius)
{
    puts("Thresholding");

    vector<Point> thresh_points;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                // grab this value
                double fi = grid_value(n, field, i, j, k);

                if (fi < thresh_value) {
                    Point p = {dx * (i + 0.5), dx * (j + 0.5), dx * (k + 0.5), fi};
                    thresh_points.push_back(p);
                }
            }
        }
    }

    const int num_points = thresh_points.size();
    printf("Finding radii of %i points.\n", num_points);

    vector<Void> candidates;

    // now for each point, we grow a spherical region until the average gets
    // above the tresh
    for (int i = 0; i < num_points; ++i) {
        if (i % 10000 == 0) { printf("progress %f\n", (double)i/num_points); }

        // grab center
        Point c = thresh_points[i];

        // just to be careful about noise -- test r = 1.1, 1.5, 2.0
        // start at 2.0 Mpc/h
        double avg;
        avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, 1.01);
        if (avg >= avg_target) { continue; }
        avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, 1.5);
        if (avg >= avg_target) { continue; }
        avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, 2.01);
        if (avg >= avg_target) { continue; }

        double r = 2.0;
        const double delta_r = 0.05;
        const double r_max = 50.0;

        while (avg < avg_target) {
            // update radius
            r += delta_r;
            // check for crazy scales...
            if (r > r_max) {
                puts("[ERROR] Exceeded r_max. Use more restrictive values for threshold and average target.");
                exit(1);
            }
            // compute new avg.
            avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, r);
        }

        Void v = {c, r};
        candidates.push_back(v);
    }

    // sort cands by radius
    const int num_candidates = candidates.size();
    printf("Sorting %i candidates\n", num_candidates);

    if (overlap_mode_radius) {
        std::sort(candidates.begin(), candidates.end(), void_radius_ge);
    }
    else {
        std::sort(candidates.begin(), candidates.end(), void_value_lt);
    }

    puts("Removing overlaps");

    // remove overlapping regions
    vector<bool> processed_indexes(num_candidates, false);
    vector<Void> voids;

    // dumb n^2 loop is fine for this many elements.
    for (int i = 0; i < num_candidates; ++i) {
        if (i % 10000 == 0) { printf("progress %f\n", (double)i/num_candidates); }

        // have we already handled this one?
        if (processed_indexes[i]) { continue; }
        // grab this void and mark it as processed
        Void vi = candidates[i];
        processed_indexes[i] = true;
        // the void to save, use vi as default values.
        Void v = vi;

        for (int j = i + 1; j < num_candidates; ++j) {
            // have we already handled this one?
            if (processed_indexes[j]) { continue; }
            Void vj = candidates[j];

            // compare distance.
            double delta_center = periodic_distance(vi.center, vj.center, l);
            bool overlap = delta_center < vi.radius + vj.radius;

            if (overlap) {
                // first mark this one as processed
                processed_indexes[j] = true;
                // now figure out which one to save.
                if (overlap_mode_radius) {
                    // replace if larger, or in the case of a tie, has a smaller center val.
                    if (vj.radius > v.radius
                        || (vj.radius == v.radius
                            && vj.center.value < v.center.value) ) {
                        v = vj;
                    }
                }
                else {
                    if (vj.center.value < v.center.value
                        || (vj.center.value == v.center.value
                            && vj.radius > v.radius)) {
                        v = vj;
                    }
                }
            }
        }

        voids.push_back(v);
    }

    printf("Pruned down to %lu voids\n", voids.size());

    return voids;
}

int
main(int argc, char **argv)
{
    if (argc != 7) {
        puts("[ERROR] bad number of args.");
        exit(1);
    }

    const string grid_path = argv[1];
    const bool thresh_mode_ge = atoi(argv[2]);
    const double in_threshold = atof(argv[3]);
    const double in_avg_target = atof(argv[4]);
    const bool overlap_mode_radius = atoi(argv[5]);
    const string output_path = argv[6];

    puts("Reading");

    // read the grid data
    double *field = new double[nn];
    FILE *f = fopen(grid_path.c_str(), "r");
    fread(&field[0], sizeof(double), nn, f);
    fclose(f);

    // check if we need to change for the threshold direction
    double threshold = in_threshold;
    double avg_target = in_avg_target;
    if (thresh_mode_ge) {
        for (int i = 0; i < nn; ++i) { field[i] *= -1.0; }
        threshold *= -1.0;
        avg_target *= -1.0;
    }

    vector<Void> voids = find_sos(field, threshold, avg_target, overlap_mode_radius);

    // check if we need to fix void values.
    if (thresh_mode_ge) {
        for (auto &v : voids) {
            v.center.value *= -1.0;
        }
    }

    puts("Writing output.");

    // write out the final voids.
    FILE *output_file = fopen(output_path.c_str(), "w");
    fprintf(output_file, "x,y,z,value,radius\n");
    for (const auto &v : voids) {
        fprintf(output_file, "%f,%f,%f,%f,%f\n",
            v.center.x, v.center.y, v.center.z, v.center.value, v.radius);
    }
    fclose(output_file);

    delete [] field;

    return 0;
}
