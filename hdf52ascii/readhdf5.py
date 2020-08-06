import matplotlib.pyplot as plt
import h5py

folder = "/proj/shuiyao/"
modelname = "l25n144-phew"
snapnum = 58
snapstr = ("000"+str(snapnum))[-3:]
fbase = folder + modelname + "/snapshot_"+snapstr
snapname = fbase + ".hdf5"
fkey = fbase + ".key"
fmcloud = fbase + ".mcloud"
flastsftime = fbase + ".lastsftime"
fvinit = fbase + ".vinit"

hf = h5py.File(snapname, "r")
print "File Read."

header = hf['Header']
attrs = header.attrs
attrs.keys()
gp = hf['PartType0']
gp.keys()
# vel = gp['Velocities']
# print len(vel)
PhEWLastSFTime = gp['PhEWLastSFTime']
PhEWKey = gp['PhEWKey']
PhEWMcloud = gp['PhEWMcloud']
PhEWVinit = gp['PhEWVinit']

print "Writing: ", fkey
f = open(fkey, "w")
for x in PhEWKey:
    f.write("%d\n" %  int(x))
f.close()
print "Writing: ", flastsftime
f = open(flastsftime, "w")
for x in PhEWLastSFTime:
    f.write("%7.5f\n" %  x)
f.close()
print "Writing: ", fmcloud
f = open(fmcloud, "w")
for x in PhEWMcloud:
    f.write("%7.5f\n" %  x)
f.close()
print "Writing: ", fvinit
f = open(fvinit, "w")
for x in PhEWVinit:
    f.write("%7.5e\n" %  x)
f.close()
