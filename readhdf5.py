import matplotlib.pyplot as plt
import h5py

# A Sample for handling HDF5 Files with Python

snapname = "./data/snapshot_l6n36_031.hdf5"

hf = h5py.File(snapname, "r")
print "File Read."

header = hf['Header'] # Read the header

attrs = header.attrs
attrs.keys() # The information contained in the header

# Load the three kinds of particles separately:
gp = hf['PartType0'] # Gas Particles
gp.keys()
dp = hf['PartType1'] # Dark Matter Particles
dp.keys()
sp = hf['PartType4'] # Star Particles
sp.keys()

# Get gas properties
pos = gp['Coordinates'] # Position of each particle
vel = gp['Velocities'] # Velocity
mass = gp['Masses'] # Mass
temp = gp['InternalEnergy'] # Temperature in K
rho = gp['Density'] # Density
ne = gp['ElectronAbundance']
# ...
print len(vel) # The total number of particles in the file

def show_snap(): # Show the x-y projection of the snapshot
    x, y = [], []
    for p in pos:
        x.append(p[0])
        y.append(p[1])
    plt.plot(x, y, "b.", markersize=2)
    print len(x), max(x), max(y)
    plt.show()
