/*

Simple void finder.
2015, Casey W. Stark.

- read grid data
- find points below threshold
- grow spheres until they hit the threshold
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

bool
comp_void_radii_desc(const Void &a, const Void &b)
{
    return a.radius >= b.radius;
}

int
main(int argc, char **argv)
{
    if (argc != 4) {
        puts("[ERROR] bad number of args.");
        exit(1);
    }

    const string grid_path = argv[1];
    const double thresh_value = atof(argv[2]);
    const string output_path = argv[3];

    // so we can test modifying this later.
    const double avg_target_value = 2.0 * thresh_value;

    puts("Reading");

    // read the grid data
    double *field = new double[nn];
    FILE *f = fopen(grid_path.c_str(), "r");
    fread(&field[0], sizeof(double), nn, f);
    fclose(f);

    puts("Thresholding");

    vector<Point> under_thresh_points;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                // grab this value
                double fi = grid_value(n, field, i, j, k);

                if (fi < thresh_value) {
                    Point p = {dx * (i + 0.5), dx * (j + 0.5), dx * (k + 0.5), fi};
                    under_thresh_points.push_back(p);
                }
            }
        }
    }

    puts("Finding radii");

    const int num_points = under_thresh_points.size();
    vector<Void> candidates;

    // now for each point, we grow a spherical region until the average gets
    // above the tresh
    for (int i = 0; i < num_points; ++i) {
        if (i % 10000 == 0) { printf("progress %f\n", (double)i/num_points); }

        // grab center
        Point c = under_thresh_points[i];

        // just to be careful about noise -- test r = 1.1, 1.5, 2.0
        // start at 2.0 Mpc/h
        double avg;
        avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, 1.01);
        if (avg >= avg_target_value) { continue; }
        avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, 1.5);
        if (avg >= avg_target_value) { continue; }
        avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, 2.01);
        if (avg >= avg_target_value) { continue; }

        double r = 2.0;
        double delta_r = 0.1;
        while (avg < avg_target_value) {
            // update radius
            r += delta_r;
            // compute new avg.
            avg = grid_spherical_average(l, n, field, c.x, c.y, c.z, r);
        }

        Void v = {c, r};
        candidates.push_back(v);
    }

    // sort cands by radius
    std::sort(candidates.begin(), candidates.end(), comp_void_radii_desc);

    const int num_candidates = candidates.size();
    printf("Found %i large voids\n", num_candidates);
    puts("Removing overlaps");

    // remove overlapping regions
    vector<bool> processed_indexes(num_candidates, false);
    vector<Void> voids;

    // dumb n^2 loop is fine for this many elements.
    for (int i = 0; i < num_candidates; ++i) {
        // have we already handled this one?
        if (processed_indexes[i]) { continue; }
        // grab this void and mark it as processed
        Void vi = candidates[i];
        processed_indexes[i] = true;
        // if we haven't processed this one, we must be adding a final void
        // use the i'th min as default values.
        Void v = vi;

        for (int j = 0; j < num_candidates; ++j) {
            // have we already handled this one?
            if (processed_indexes[j]) { continue; }
            Void vj = candidates[j];

            // compare distance.
            double delta_center = periodic_distance(vi.center, vj.center, l);
            bool overlap = delta_center < vi.radius + vj.radius;

            if (overlap) {
                // DEBUG
                if (vj.radius > 11.0) {
                    printf("%i %f %f %f %f - %i %f %f %f %f, dist %f\n",
                        i, vi.center.x, vi.center.y, vi.center.z, vi.radius,
                        j, vj.center.x, vj.center.y, vj.center.z, vj.radius,
                        delta_center);
                }

                // first mark this one as processed
                processed_indexes[j] = true;
                // now figure out which one to save.
                if (vj.radius > v.radius) {
                    v = vj;
                }
            }
        }

        voids.push_back(v);
    }

    printf("Pruned down to %lu voids\n", voids.size());
    puts("Writing output.");

    // write out the final voids.
    FILE *output_file = fopen(output_path.c_str(), "w");
    fprintf(output_file, "x,y,z,value,radius\n");
    for (const auto &v : voids) {
        fprintf(output_file, "%f,%f,%f,%f,%f\n",
            v.center.x, v.center.y, v.center.z, v.center.value, v.radius);
    }
    fclose(output_file);

    return 0;
}
