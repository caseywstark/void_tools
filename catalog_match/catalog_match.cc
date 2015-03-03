/*
2015, Casey W. Stark <caseywstark@gmail.com>
Compare two void catalogs.
- Read catalog A and B.
- For each A void,
  - Find nearest B void.
  - If B void within r_A, save it.
- Write joined catalog with format like:
  ax, ar, av, b_id, bx, br, bv, delta_x.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "utils.h"

using namespace std;

const double l = 256.0;

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

std::vector<Void>
read_voids(const std::string path)
{
    // Prep return value.
    std::vector<Void> voids;

    // Open and check file.
    std::ifstream input_file(path.c_str());
    if (input_file.is_open()) {
        // line obj to read into.
        std::string line;
        // read past headers
        std::getline(input_file, line);

        // now read for centers and push onto return vec.
        while (std::getline(input_file, line)) {
            std::vector<std::string> t = split_string(line, ',');

            // TODO
            // hard coding position indexes for now.
            Point p = {stod(t[0]), stod(t[1]), stod(t[2]), stod(t[3])};
            Void v = {p, stod(t[4])};
            voids.push_back(v);
        }
    }
    else {
        std::cerr << "Could not read " << path << std::endl;
        exit(1);
    }

    return voids;
}


int
main(int argc, char **argv)
{
    if (argc != 4) {
        puts("[ERROR] bad number of args.");
        exit(1);
    }

    const string cat_a_path = argv[1];
    const string cat_b_path = argv[2];
    const string output_path = argv[3];

    puts("Reading");

    vector<Void> voids_a = read_voids(cat_a_path);
    vector<Void> voids_b = read_voids(cat_b_path);

    puts("Matching");

    // default case for no match.
    Void no_match_void;

    // vectors to save matches
    // parallel array style isn't the best, but it's quicker for now...
    vector<int> match_ids;
    vector<double> match_delta_xs;
    vector<Void> match_voids;

    for (const auto &void_a : voids_a ) {

        bool found_a_match = false;
        int match_id = 0;
        double match_delta_x = 100.0; // some large value.

        // dumb N^2 loop
        // explicit loop since the index is the ID.
        for (size_t j = 0; j < voids_b.size(); ++j) {
            // get separation.
            double delta_x = periodic_distance(void_a.center, voids_b[j].center, l);
            if (delta_x < void_a.radius) {
                found_a_match = true;
                if (delta_x < match_delta_x) {
                    match_id = j;
                    match_delta_x = delta_x;
                }
            }
        }

        // did we get a match?
        if (found_a_match) {
            match_ids.push_back(match_id);
            match_delta_xs.push_back(match_delta_x);
            match_voids.push_back( voids_b[match_id] );
        }
        else {
            match_ids.push_back(-1);
            match_delta_xs.push_back(-1);
            match_voids.push_back(no_match_void);
        }
    }

    puts("Writing output.");

    // write out the final voids.
    FILE *output_file = fopen(output_path.c_str(), "w");
    fprintf(output_file, "x,y,z,value,radius,match_id,delta_x,mx,my,mz,mvalue,mradius\n");
    for (size_t i = 0; i < voids_a.size(); ++i) {
        Void va = voids_a[i];
        Void vb = match_voids[i];
        fprintf(output_file, "%f,%f,%f,%f,%f,%i,%f,%f,%f,%f,%f,%f\n",
            va.center.x, va.center.y, va.center.z, va.center.value, va.radius,
            match_ids[i], match_delta_xs[i],
            vb.center.x, vb.center.y, vb.center.z, vb.center.value, vb.radius);
    }
    fclose(output_file);

    return 0;
}
