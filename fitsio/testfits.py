from astropy.table import Table
from astropy.io import fits

fname = "table1.fits"
# t = Table([[1, 2], [4, 5], [7, 8]], names=('a', 'b', 'c'))
# t.write('table1.fits', format='fits')
t = fits.open(fname)
h = fits.getheader(fname)
t.info()
tab = t[1]
# tab.dump("table1.txt")
d = fits.getdata(fname, 0)
d['a']
