from astropy.table import Table
from astropy.io import fits
from ctypes import *
from astroconst import pc, ac
from numpy import sqrt, pi, array
import random

# fname = "snap_p6n36fof_031.bin"
# faux = "snap_p6n36fof_031.aux"
# fname = "snap_p50n288gw_108.bin"
# faux = "snap_p50n288gw_108.aux"
MODEL = "p50n576"
if(MODEL == "p50n576"):
    fbase = "/proj/shuiyao/p50n576fi/"
    fname = fbase+"snap_p50n576zw_108.bin"
    faux = fbase+"snap_p50n576zw_108.aux"
    outbase = "/proj/shuiyao/fits/p50n576fi/"
    fout = outbase+"snap_p50n576_z0.fits"
    LBOX = 50000.0
    HUBBLEPARAM = 0.7
    ASCALE = 1.0
    fac_density = 1.0
    fac_length = 1.0
if(MODEL == "comaII"):
    fbase = "/proj/shuiyao/comaII/snapdir_201/"
    fname = fbase+"snap_gw_201.bin"
    faux = fbase+"snap_gw_201.aux"
    outbase = "/proj/shuiyao/fits/comaII/"
    fout = outbase+"snap_comaII_z0.fits"
    LBOX = 500.0
    HUBBLEPARAM = 0.72
    ASCALE = 1.0
    fac_density = 1.e-9
    fac_length = 1.e3

def convert_unit_to_physical(arr, fac):
    arr = array(arr)
    return arr * fac

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

class tipsy_header(Structure):
    _fields_ = [('time', c_double),\
                ('nbodies', c_int),\
                ('ndim', c_int),\
                ('nsph', c_int),\
                ('ndark', c_int),\
                ('nstar', c_int),\
                ('dummy', c_int)]

class tipsy_position(Structure):
    _fields_ = [('x', c_float), ('y', c_float), ('z', c_float)]
class tipsy_velocity(Structure):
    _fields_ = [('x', c_float), ('y', c_float), ('z', c_float)]
class tipsy_metallicity(Structure):
    _fields_ = [('C', c_float), ('O', c_float), ('Si', c_float), ('Fe', c_float)]

class tipsy_gas(Structure):
    _fields_ = [('mass', c_float),\
                ('pos', c_float * 3),\
                ('vel', c_float * 3),\
                ('rho', c_float),\
                ('temp', c_float),\
                ('hsmooth', c_float),\
                ('metals', c_float),\
                ('phi', c_float)]
    def unit_conversion(self, unit):
        self.mass *= unit.mass
        for i in [0,1,2]: self.pos[i] = (self.pos[i] + 0.5) * unit.length
        for i in [0,1,2]: self.vel[i] = self.vel[i] * unit.velocity
        self.rho *= unit.density
        self.hsmooth *= unit.length

class tipsy_dark(Structure):
    _fields_ = [('mass', c_float),\
                ('pos', c_float * 3),\
                ('vel', c_float * 3),\
                ('eps', c_float),\
                ('phi', c_float)]

class tipsy_aux_gas(Structure):
    _fields_ = [('metals', c_float * 4),\
                ('sfr', c_float),\
                ('tmax', c_float),\
                ('delayt', c_float),\
                ('ne', c_float),\
                ('nh', c_float),\
                ('nspawn', c_int),\
                ('nrec', c_short)]

class tipsy_aux_star(Structure):
    _fields_ = [('metals', c_float * 4),\
                ('age', c_float),\
                ('tmax', c_float),\
                ('nspawn', c_int),\
                ('nrec', c_short)]


class tipsy_star(Structure):
    _fields_ = [('mass', c_float),\
                ('pos', c_float * 3),\
                ('vel', c_float * 3),\
                ('metals', c_float),\
                ('tform', c_float),\
                ('eps', c_float),\
                ('phi', c_float)]
    def unit_conversion(self, unit):
        self.mass *= unit.mass
        for i in [0,1,2]: self.pos[i] = (self.pos[i] + 0.5) * unit.length
        for i in [0,1,2]: self.vel[i] = self.vel[i] * unit.velocity

class tipsy_binary():
    def __init__(self, fname, auxname=''):
        self.file_name = fname
        if(auxname == ''):
            self.aux_name = fname[:-4]+".aux"
    def load_header(self):
        with open(self.file_name, 'rb') as fin:
            h = tipsy_header()
            fin.readinto(h)
        self.header = h
        return self.header
    def load_data(self, load_star = False, \
                  isph_min=0, isph_max=576**3, \
                  istar_min=0, istar_max=576**3):
        print "Loading TIPSY file: ", self.file_name
        print "Including SPH from %d to %d." % (isph_min, isph_max)
        print "Including STAR from %d to %d." % (istar_min, istar_max)        
        fin = open(self.file_name, 'rb')
        faux = open(self.aux_name, 'rb')
        h = tipsy_header()
        fin.readinto(h)
        self.header = h
        self.gas = []
        self.star = []
        self.aux_star = []
        for i in range(self.header.nsph):
            gp = tipsy_gas()
            fin.readinto(gp)
            # For efficiency, we should do unit conversion at the end
            # gp.unit_conversion(unit)
            auxgp = tipsy_aux_gas()
            faux.readinto(auxgp)
            if(i < isph_min): continue
            if(i >= isph_max): continue
            if(random.random() >= 1./2.): continue
            self.gas.append(gp)
        print "GAS Loaded."
        if(load_star == True):
            for i in range(self.header.ndark):
                dp = tipsy_dark()
                fin.readinto(dp)
            print "DARK Loaded."            
            for i in range(self.header.nstar):
                sp = tipsy_star()
                fin.readinto(sp)
                auxsp = tipsy_aux_star()
                faux.readinto(auxsp)
                # sp.unit_conversion(unit)
                if(i < istar_min): continue                                
                if(i >= istar_max): continue 
                if(random.random() >= 1./2.): continue            
                self.star.append(sp)
                self.aux_star.append(auxsp)
            print "STAR Loaded."
        fin.close()
        faux.close()
    def write_to_fits(self, fitsname=""):
        if(fitsname == ""):
            fitsname = outbase+self.file_name.split('.')[0]+".fits"
        print "Writing FITS file: ", fitsname
        x, y, z, vx, vy, vz, mass, rho, temp, metals = \
        [], [], [], [], [], [], [], [], [], []
        for i in range(len(self.gas)):
            x.append(self.gas[i].pos[0])
            y.append(self.gas[i].pos[1])
            z.append(self.gas[i].pos[2])
            vx.append(self.gas[i].vel[0])
            vy.append(self.gas[i].vel[1])
            vz.append(self.gas[i].vel[2])
            mass.append(self.gas[i].mass)
            rho.append(self.gas[i].rho * unit.density * fac_density)
            temp.append(self.gas[i].temp)
            metals.append(self.gas[i].metals)
        for i in range(len(self.star)):
            x.append(self.star[i].pos[0])
            y.append(self.star[i].pos[1])
            z.append(self.star[i].pos[2])
            vx.append(self.star[i].vel[0])
            vy.append(self.star[i].vel[1])
            vz.append(self.star[i].vel[2])
            mass.append(self.star[i].mass)
            rho.append(self.aux_star[i].age / 1.e9) # Formation time for stars
            temp.append(-1.0)
            metals.append(self.star[i].metals)
        x = convert_unit_to_physical(x, unit.length) + 0.5 * unit.length
        y = convert_unit_to_physical(y, unit.length) + 0.5 * unit.length
        z = convert_unit_to_physical(z, unit.length) + 0.5 * unit.length
        x = x * fac_length
        y = y * fac_length
        z = z * fac_length        
        vx = convert_unit_to_physical(vx, unit.velocity)
        vy = convert_unit_to_physical(vy, unit.velocity)
        vz = convert_unit_to_physical(vz, unit.velocity)
        mass = convert_unit_to_physical(mass, unit.mass)
        # rho = convert_unit_to_physical(rho, unit.density)
        tab = Table([x, y, z, vx, vy, vz, mass, rho, temp, metals], \
                    names = ('x[pkpc]', 'y[pkpc]', 'z[pkpc]', 'vx[km/s]', 'vy[km/s]', 'vz[km/s]',
                             'mass[Msolar]', 'density[g/cm^3](formation time[Gyr])', 'temperature[K]', 'metallicity'))
        tab.write(fitsname, format='fits', overwrite=True)

unit = tipsy_units(ASCALE, LBOX, HUBBLEPARAM, unit_system='default')
NFILES = 8
tb = tipsy_binary(fname)
h = tb.load_header()
sph_step = h.nsph / NFILES
star_step = h.nstar / NFILES
print "NSPH = ", h.nsph, "NSTAR = ", h.nstar
isph_range = [0, sph_step-1]
istar_range = [0, star_step-1]
for nread in range(NFILES):
    dat = tb.load_data(load_star = True, \
                     isph_min=isph_range[0], isph_max=isph_range[1], \
                     istar_min=istar_range[0], istar_max=istar_range[1])    
    tb.write_to_fits(fitsname = fout+"."+str(nread))
    isph_range[0] += sph_step
    isph_range[1] += sph_step
    istar_range[0] += star_step
    istar_range[1] += star_step
    tb = tipsy_binary(fname)
    h = tb.load_header()

# ffits = "snap_p50n288gw_108.fits"
# ffits = "snap_p6n36fof_031.fits"
# with fits.open(ffits, 'update') as f:
#     f[1].header['TTYPE1'] = 'x[pkpc]'
#     f[1].header['TTYPE2'] = 'y[pkpc]'
#     f[1].header['TTYPE3'] = 'z[pkpc]'
# print(repr(fits.getheader(ffits), 1))    
