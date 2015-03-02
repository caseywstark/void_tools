
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

#include <mpi.h>

// simulation constants
const double l = 256.0;
const long n = 2560;
const long np = n*n*n;
const double dx = l / n;

std::vector<std::string>
split_string(const std::string str, const char split_char)
{
    // make return obj
    std::vector<std::string> result;

    // split and push each back to result.
    std::stringstream str_stream(str);
    std::string element;
    while(std::getline(str_stream, element, split_char))
    {
        result.push_back(element);
    }

    return result;
}

std::vector< std::array<double, 3> >
read_centers(const std::string path)
{
    // Prep return value.
    std::vector< std::array<double, 3> > centers;

    // Open and check file.
    std::ifstream centers_file(path.c_str());
    if (centers_file.is_open()) {
        // line obj to read into.
        std::string line;

        // read past headers
        std::getline(centers_file, line);

        // now read for centers and push onto return vec.
        while (std::getline(centers_file, line)) {
            std::vector<std::string> t = split_string(line, ',');

            // TODO
            // hard coding position indexes for now.
            std::array<double, 3> c = {stod(t[0]), stod(t[1]), stod(t[2])};
            // DEBUG
            // just to be safe...
            if (c[0] < 0.0 || c[0] >= l || c[1] < 0.0 || c[1] >= l
                || c[2] < 0.0 || c[2] >= l) {
                std::cerr << "bad position: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
                exit(1);
            }

            centers.push_back(c);
        }
    }
    else {
        std::cerr << "Could not read " << path << std::endl;
        exit(1);
    }

    return centers;
}


int
main(int argc, char **argv)
{
    // Initialize MPI stuff
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    int master_rank = 0;

    // command line
    const std::string centers_path = argv[1];
    const std::string field_path = argv[2];
    const std::string pre_path = argv[3];

    // figure out decomp.
    if (n % mpi_size != 0) {
        printf("Bad mpi size %i for grid %li\n", mpi_size, n);
    }

    //
    // Read halo positions.
    //

    if (mpi_rank == 0) {
        printf("# Reading centers %s\n", centers_path.c_str());
    }

    vector< array<double, 3> > centers = read_centers(centers_path);

    if (mpi_rank == 0) {
        printf("# Found %lu centers\n", centers.size());
    }


    //
    // Read field values.
    //

    if (mpi_rank == 0) {
        printf("# Reading %s\n", field_path.c_str());
    }

    const long n_local = n / mpi_size;
    const long np_local = n_local * n * n;
    const long local_ix0 = mpi_rank * n_local;

    double *f = new double[np_local];

    FILE *field_file = fopen(field_path.c_str(), "r");
    fseek(field_file, sizeof(double) * local_ix0 * n * n, SEEK_SET);
    fread(f, sizeof(double), np_local, field_file);
    fclose(field_file);

    //
    // Compute field mean and variance.
    //

    double f_sum_local = 0.0;
    double f2_sum_local = 0.0;

    for (long i = 0; i < np_local; ++i) {
        f_sum_local += f[i];
        f2_sum_local += f[i]*f[i];
    }

    double f_sum;
    MPI_Allreduce(&f_sum_local, &f_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    double f2_sum;
    MPI_Allreduce(&f2_sum_local, &f2_sum, 1, MPI_DOUBLE, MPI_SUM, comm);

    double f_mean = f_sum / np;
    double f2_mean = f2_sum / np;

    if (mpi_rank == master_rank) {
        char output_path[1000];
        sprintf(output_path, "%sfield_var.txt", pre_path.c_str());
        FILE *outfile = fopen(output_path, "w");
        fprintf(outfile, "%e %e\n", f_mean, sqrt(f2_mean - f_mean * f_mean));
        fclose(outfile);
    }

    //
    // start accumulating radial bins
    //

    if (mpi_rank == 0) { puts("# Finding local halos."); }

    // 10 Mpc/h in 256 Mpc/h
    const double r_max = 20.0;
    const int num_bins = 100;
    const double bin_dr = r_max / num_bins;

    long *bin_counts_local = new long[num_bins];
    long *bin_counts = new long[num_bins];
    double *bin_f_sums_local = new double[num_bins];
    double *bin_f_sums = new double[num_bins];
    double *bin_r_sums_local = new double[num_bins];
    double *bin_r_sums = new double[num_bins];

    // iterate over groups
    for (size_t i_cent = 0; i_cent < centers.size(); ++i_cent) {
        // reset local stuff
        for (int i = 0; i < num_bins; ++i) {
            bin_counts_local[i] = 0;
            bin_f_sums_local[i] = 0.0;
            bin_r_sums_local[i] = 0.0;
        }

        // get pc position
        double cx = centers[i_cent][0];
        double cy = centers[i_cent][1];
        double cz = centers[i_cent][2];

        long ix0 = floor( (cx - r_max) / dx );
        long ix1 = ceil( (cx + r_max) / dx );
        long iy0 = floor( (cy - r_max) / dx );
        long iy1 = ceil( (cy + r_max) / dx );
        long iz0 = floor( (cz - r_max) / dx );
        long iz1 = ceil( (cz + r_max) / dx );

        // iterate over points.
        for (long ix = ix0; ix <= ix1; ++ix) {
            // periodic fix
            long ixw = ix;
            if (ixw < 0) { ixw += n; }
            if (ixw >= n) { ixw -= n; }

            long ix_local = ixw - local_ix0;

            // check if this is even local before continuing.
            if (ix_local >= 0 && ix_local < n_local) {
                double x = dx * (ix + 0.5);
                double drx = x - cx;
                double drx2 = drx * drx;

                for (long iy = iy0; iy <= iy1; ++iy) {
                    long iyw = iy;
                    if (iyw < 0) { iyw += n; }
                    if (iyw >= n) { iyw -= n; }

                    double y = dx * (iy + 0.5);
                    double dry = y - cy;
                    double drxy2 = drx2 + dry * dry;

                    for (long iz = iz0; iz <= iz1; ++iz) {
                        long izw = iz;
                        if (izw < 0) { izw += n; }
                        if (izw >= n) { izw -= n; }

                        double z = dx * (iz + 0.5);
                        double drz = z - cz;
                        double r = sqrt(drxy2 + drz * drz);
                        int ibin = r / bin_dr;

                        if (ibin < num_bins) {
                            // got a bin point
                            long ii = (ix_local * n + iyw) * n + izw;
                            bin_counts_local[ibin] += 1;
                            bin_f_sums_local[ibin] += f[ii];
                            bin_r_sums_local[ibin] += r;
                        }
                    }
                }
            }
        }

        // done with binning for this PC.
        MPI_Allreduce(bin_counts_local, bin_counts, num_bins, MPI_LONG, MPI_SUM, comm);
        MPI_Allreduce(bin_f_sums_local, bin_f_sums, num_bins, MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(bin_r_sums_local, bin_r_sums, num_bins, MPI_DOUBLE, MPI_SUM, comm);

        // write this pc.
        // make sure file is empty
        if (mpi_rank == 0) {
            char output_path[1000];
            sprintf(output_path, "%sprofile_%06d.txt", pre_path.c_str(), (int)i_cent);
            FILE *outfile = fopen(output_path, "w");
            fprintf(outfile, "radius,field_value\n");
            for (int i = 0; i < num_bins; ++i) {
                if (bin_counts[i] > 0) {
                    fprintf(outfile, "%e,%e\n",
                        bin_r_sums[i] / bin_counts[i],
                        bin_f_sums[i] / bin_counts[i]);
                }
                else {
                    fprintf(outfile, "0.0,0.0\n");
                }
            }
            fclose(outfile);
        }

        // done with this center
    }

    MPI_Finalize();
    return 0;
}
