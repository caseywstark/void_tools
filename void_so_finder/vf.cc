/*

Simple void finder.
2015, Casey W. Stark.

- read grid data
- find local mins
- grow mins until they hit the threshold
- remove overlapping mins

*/

#include <array>
#include <string>
#include <vector>

#include "utils.h"

using namespace std;

// grid resolution
// define at compile time for unrolling loops.
const int ng = 256;
const int nn = ng * ng * ng;

int
main(int argc, char **argv)
{
    const string grid_path = argv[1];
    const double thresh_value = atof(argv[2]);
    const string output_path = argv[3];

    // so we can test modifying this later.
    const double avg_target_value = thresh_value;

    puts("Reading");

    // read the grid data
    vector<double> field = read_grid_data(grid_path, nn);

    puts("Thresholding");

    vector<Point> under_thresh_points;

    for (int i = 0; i < ng; ++i) {
        for (int j = 0; j < ng; ++j) {
            for (int k = 0; k < ng; ++k) {
                // grab this value
                double fi = grid_value(ng, field, i, j, k);

                if (fi < thresh_value) {
                    Point p = {i, j, k, fi};
                    under_thresh_points.push_back(p);
                }
            }
        }
    }

    puts("Finding mins");

    // find local mins, store in vector of points
    vector<Void> local_min_voids;
    for (const auto &p : under_thresh_points) {
        const int i = p.x;
        const int j = p.y;
        const int k = p.z;
        const double fi = p.value;

        // now check if this point is smaller than neighbors.
        bool is_min = true;

        is_min &= fi < grid_value_periodic(ng, field, i, j, k - 1);
        is_min &= fi < grid_value_periodic(ng, field, i, j, k + 1);
        is_min &= fi < grid_value_periodic(ng, field, i, j - 1, k);
        is_min &= fi < grid_value_periodic(ng, field, i, j + 1, k);
        is_min &= fi < grid_value_periodic(ng, field, i - 1, j, k);
        is_min &= fi < grid_value_periodic(ng, field, i + 1, j, k);

        /*
        is_min &= fi < grid_value_periodic(ng, field, i, j - 1, k - 1);
        is_min &= fi < grid_value_periodic(ng, field, i, j - 1, k + 1);
        is_min &= fi < grid_value_periodic(ng, field, i, j + 1, k - 1);
        is_min &= fi < grid_value_periodic(ng, field, i, j + 1, k + 1);

        is_min &= fi < grid_value_periodic(ng, field, i - 1, j, k - 1);
        is_min &= fi < grid_value_periodic(ng, field, i - 1, j, k + 1);
        is_min &= fi < grid_value_periodic(ng, field, i + 1, j, k - 1);
        is_min &= fi < grid_value_periodic(ng, field, i + 1, j, k + 1);

        is_min &= fi < grid_value_periodic(ng, field, i - 1, j - 1, k);
        is_min &= fi < grid_value_periodic(ng, field, i - 1, j + 1, k);
        is_min &= fi < grid_value_periodic(ng, field, i + 1, j - 1, k);
        is_min &= fi < grid_value_periodic(ng, field, i + 1, j + 1, k);
        */

        if (is_min) {
            // we have a local min under the threshold.
            Void v = {p, 0.0};
            local_min_voids.push_back(v);
        }
    }

    puts("Finding radii");

    vector<Void> large_voids;

    // now for each min point, we grow a spherical region until the average
    for (size_t i = 0; i < local_min_voids.size(); ++i) {
        Void v = local_min_voids[i];
        Point center = v.center;

        // before we try to find the size, make sure it's at least 2 mpc/h across
        const double r_min = 2.0;
        double avg_min = grid_spherical_average<ng>(field, center, r_min + 0.01);
        if (avg_min >= avg_target_value) {
            continue;
        }

        // Bisection is fairly slow, but this is fast enough.
        // bisection params
        double r = 0.0;
        double r_lo = r_min;
        double r_hi = 50.0;
        const double r_tol = 0.02;
        const int max_iter = 100;

        int iter = 0;
        for (iter = 0; iter < max_iter; ++iter) {
            // check for stop
            double delta_r_range = r_hi - r_lo;
            if (delta_r_range < r_tol) { break; }

            // new radius.
            r = 0.5 * (r_lo + r_hi);
            // get spherical average
            double avg = grid_spherical_average<ng>(field, center, r);
            // update bounds
            if (avg < avg_target_value) { r_lo = r; }
            else { r_hi = r; }
        }

        printf("frac %f radius %f\n", (double)i / local_min_voids.size(), r);

        if (iter == max_iter) {
            // Did not converge in max_iter!
            printf("[ERROR] radius bisection did not converge.");
            exit(1);
        }

        // success.
        v.radius = r;
        large_voids.push_back(v);
    }

    printf("Found %lu large voids\n", large_voids.size());
    puts("Removing overlaps");

    // remove overlapping regions
    vector<bool> processed(large_voids.size(), false);
    vector<Void> voids;

    // dumb n^2 loop is fine for this many elements.
    for (size_t i = 0; i < large_voids.size(); ++i) {

        // have we already handled this one?
        if (processed[i]) { continue; }

        Void vi = large_voids[i];
        processed[i] = true;

        // also check if the radius was never set.
        if (vi.radius == 0.0) { continue; }

        // if we haven't processed this one, we must be adding a final void
        // use the i'th min as default values.
        Void v = vi;

        for (size_t j = 0; j < large_voids.size(); ++j) {
            // have we already handled this one?
            if (processed[j]) { continue; }

            Void vj = large_voids[j];

            // compare distance.
            double delta_center = periodic_distance(vi.center, vj.center, (double)ng);
            bool overlap = delta_center < vi.radius + vj.radius;
            if (overlap) {
                // first mark this one as processed
                processed[j] = true;
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
    for (const auto &v : voids) {
        fprintf(output_file, "%i %i %i %f %f\n",
            v.center.x, v.center.y, v.center.z, v.center.value, v.radius);
    }
    fclose(output_file);

    return 0;
}
