import pygad
import matplotlib.pyplot as plt

# snapname = "snap_m12.5n128_125.hdf5"
snapname = "/scratch/shuiyao/data/m25n256/snap_m25n256_125.hdf5"
s = pygad.Snap(snapname)
print s.loadable_blocks()
# snap, halo = pygad.tools.prepare_zoom(s)
# R200, M200 = pygad.analysis.virial_info(snap)
# fig, ax, cbar = pygad.plotting.image(snap.gas, extent='5 Mpc')
# ax.add_artist(plot.Circle([0,0], R200, facecolor='none', edgecolor='w'))
# plt.draw()
s.cosmology
s.redshift
s.scale_factor
s.parts
s['ID']
s.gas['rho']
s.gas['temp'] # Will be derived automatically
