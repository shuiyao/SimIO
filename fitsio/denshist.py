import matplotlib.pyplot as plt
from numpy import array
import ioformat
from numpy import sqrt, log, pi
from astroconst import pc,ac
from numpy import logspace, histogram


# fbase = "/proj/shuiyao/p50n576fi/"
# fname = fbase+"snap_p50n576zw_108.rho"
fbase = "/proj/shuiyao/comaII/snapdir_201/"
fname = fbase+"snap_gw_201.rho"

class tipsy_units():
    def __init__(self, ascale, boxsize, hubble=0.7, unit_system = 'default'):
        # TIPSY conversion factors
        # Comoving, cgs
        self.time = sqrt(8.0 * pi / 3.0) * ac.mpc / (100. * hubble * 1.e5)
        self.density= 1.8791e-29 * hubble ** 2
        self.length= boxsize * ac.kpc
        self.mass = self.density * self.length ** 3 / hubble ** 2
        self.velocity = self.length / self.time
        if(unit_system == 'default'):
            self.density = self.density * hubble ** 2 / ascale ** 3 # pcgs
            self.length = self.length / ac.kpc / hubble * ascale # pkpc
            self.velocity = self.velocity * hubble / ascale ** 1.5 / 1.e5 # peculiar velocity, km/s
            self.mass = self.mass / hubble / ac.msolar # cMsolar

# unit = tipsy_units(1.0, 500, 0.72, unit_system='default')
# print "Reading density from: ", fname
# rho = ioformat.rcol(fname, [0])
# print "Convert units ... "
# rho = array(rho) * unit.density
# rho = array(rho) * unit.density / 1.e9 # Coma II
# hist, edges = histogram(rho, bins=logspace(-30., -22., 100))
# mid = 0.5 * (edges[1:] + edges[:-1])

def show():
    plt.plot(mid, hist, "b.-")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Density [g/cm^3]")
    plt.ylabel("Ncount")
    plt.title("comaII, z=0")
    plt.show()
