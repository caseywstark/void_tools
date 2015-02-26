"""
This script reads the ZOBOV input positions, the output volumes, zones, void
hierarchy, and catalog, and adds the effective radius, core particle position,
and volume-weighted centroid, writing a nice new catalog file.

Note that there are some hard-coded factors and paths specific to my analysis
setup. Be sure to edit for your own use.

"""

print "importing"

import os
from collections import defaultdict

import numpy as np
import pandas as pd

def read_positions(path):
    f = open(path, "r")
    n = np.fromfile(f, count=1, dtype="l")
    return np.fromfile(f, dtype="d").reshape(n, 3)

def read_volumes(path):
    f = open(path, "r")
    n = np.fromfile(f, count=1, dtype="l")
    return np.fromfile(f, dtype="d").reshape(n)

def compute_centroid(positions, volumes, zone_pid_map, zone_ids):
    # grab particle ids over all zones
    pids = []
    for z in zone_ids:
        pids.extend( zone_pid_map[z] )

    # grab volumes
    v_void = volumes[pids]
    # grab positions
    x_void = positions[pids]
    # compute centroid
    centroid = (v_void[:, np.newaxis] * x_void).sum(axis=0) / v_void.sum()

    return centroid

def read_zone_pid_map(path):
    zones = np.fromfile(path, dtype="i")

    # create default dict going from zone id (int) to list of particle indexes (list)
    zone_pid_map = defaultdict(list)
    for i, z in enumerate(zones):
        zone_pid_map[z].append(i)

    return zone_pid_map

# path setup
ds = 10
base_path = "catalog_ds%i/default" % ds

rt_cat_path = os.path.join(base_path, "ds%i.0.200000.rthresh.2.0.txt" % ds)
rt_void_zone_path = os.path.join(base_path, "ds%i.0.200000.rthresh.2.0.void" % ds)

pos_path = os.path.join(base_path, "../pos_ds%i.bin" % ds)
vol_path = os.path.join(base_path, "ds%i.vol" % ds)
zone_path = os.path.join(base_path, "zones.bin")

output_path = "ds%i_catalog.txt" % ds

print "reading rthresh catalog"
voids = pd.read_csv(rt_cat_path, sep=" ", header=None)

print "fixing column names"
voids.drop(0, axis=1, inplace=True)
voids.columns = ["void_id", "core_particle_id", "core_density", "zone_volume",
    "zone_num_particles", "num_zones", "volume", "num_particles",
    "density_contrast", "prob"]

print "fixing volume and adding radius"
num_particles = np.fromfile(pos_path, count=1, dtype="l")
voids["volume"] *= 256.0**3 / num_particles
voids["r_eff"] = ( 3.0 / (4.0 * np.pi) * voids["volume"] )**(1.0/3)

print "Cutting voids down"
orig_num_voids = len(voids)
print "Original number:", orig_num_voids

cut = voids["density_contrast"] > 2.0
voids = voids[cut]
print "Fraction left after density contrast cut:", (float(len(voids)) / orig_num_voids)

cut = voids["core_density"] < 0.2
voids = voids[cut]
print "Fraction left after density cut:", (float(len(voids)) / orig_num_voids)

cut = voids["r_eff"] > 2.0
voids = voids[cut]
print "Fraction left after r_eff cut:", (float(len(voids)) / orig_num_voids)

num_voids = len(voids)
print "Volume covering fraction:", ( (voids["volume"]).sum() / 256.0**3 )
print "Number of voids:", num_voids

print "reading positions"
positions = read_positions(pos_path)

print "Saving core positions"
void_positions = positions[ voids["core_particle_id"] ]
voids["core_x"] = void_positions[:, 0]
voids["core_y"] = void_positions[:, 1]
voids["core_z"] = void_positions[:, 2]

print "read zone list into hash"
void_zones_file = open(rt_void_zone_path, "r")
void_zones_map = {}
for line in void_zones_file:
    # split line into chunks, remove \n, and cast to ints
    zones = map(int, line[:-2].split(" "))
    void_zones_map[zones[0]] = zones

print "get void zone lists in same order"
void_zones_list = []
for void_id in voids["void_id"].values:
    void_zones_list.append( void_zones_map[void_id] )

print "creating zone pid map"
zone_pid_map = read_zone_pid_map(zone_path)

print "reading volumes"
volumes = read_volumes(vol_path)

print "computing centroids"
centroids = np.zeros((num_voids, 3))
for i, zones in enumerate(void_zones_list):
    if i % 100 == 0:
        print "%f" % (float(i) / num_voids)
    c = compute_centroid(positions, volumes, zone_pid_map, zones)
    centroids[i, :] = c

print "Adding centroids to catalog"
voids["cent_x"] = centroids[:, 0]
voids["cent_y"] = centroids[:, 1]
voids["cent_z"] = centroids[:, 2]

print "Writing catalog"
voids.to_csv(output_path, index=False)

