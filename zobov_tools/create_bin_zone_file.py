"""
This script reads the ASCII formatted ZOBOV zone file and writes out a binary
version to match the other particle record files.

I should have modified the ZOBOV jozov app source to write in binary in
the first place, but I had already run it on several datasets by the time I
needed to use the zone file.

"""

import os
import numpy as np

# path setup
ds = 10
base_path = "catalog_ds10/default"

zone_path = os.path.join(base_path, "ds%i.0.200000.zone" % ds)
bin_zone_path = os.path.join(base_path, "zones.bin")

f = open(zone_path, "r")
counts = map(int, f.readline().split(" "))
zones = np.array([int(line) for line in f], dtype="i")

zones.tofile(bin_zone_path)

