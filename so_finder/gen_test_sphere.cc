
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <string>
#include <vector>

#include "utils.h"

using namespace std;

// grid resolution
// define at compile time for unrolling loops.
const int n = 256;
const int nn = n * n * n;

const double l = 256.0;
const double grid_dx = l / n;


int
main()
{
    const string output_path = "test_sphere.bin";

    double *field = new double[nn];

    const double r = 10.0;
    const double c = l / 2;

    for (int i = 0; i < n; ++i) {
        double x = grid_dx * (i + 0.5);
        double dx = x - c;
        for (int j = 0; j < n; ++j) {
            double y = grid_dx * (j + 0.5);
            double dy = y - c;
            for (int k = 0; k < n; ++k) {
                double z = grid_dx * (k + 0.5);
                double dz = z - c;

                double ri = sqrt(dx*dx + dy*dy + dz*dz);
                int ii = (i * n + j) * n + k;

                if (ii == 128 * n*n + 128 * n + 128) {
                    field[ii] = 0;
                }
                else if (ri < r) {
                    field[ii] = 0.5;
                }
                else {
                    field[ii] = 1.5;
                }
            }
        }
    }

    // write out the final voids.
    FILE *output_file = fopen(output_path.c_str(), "w");
    fwrite(field, sizeof(double), nn, output_file);
    fclose(output_file);

    delete [] field;

    return 0;
}
