import struct

fbase = "/proj/shuiyao/l12n144-phew-movie-200/"
#fin = open("snap_p6n36phew_000.0")
fin = open(fbase+"snapshot_200.gad")
print "BlockSize = ", struct.unpack('i', fin.read(4))
print "Npart: "
for i in range(6):
    print struct.unpack('i', fin.read(4))
print "Mass: "
for i in range(6):
    print struct.unpack('d', fin.read(8))
print "Time = ", struct.unpack('d', fin.read(8))
fin.read(256 - 24 - 48 - 8)
print "BlockSize = ", struct.unpack('i', fin.read(4))

# POS
blksize = struct.unpack('i', fin.read(4))[0]
print "[POS] BlockSize = ", blksize
fin.read(blksize)
blksize = struct.unpack('i', fin.read(4))[0]
# VEL
blksize = struct.unpack('i', fin.read(4))[0]
print "[VEL] BlockSize = ", blksize
fin.read(blksize)
blksize = struct.unpack('i', fin.read(4))[0]
# ID
blksize = struct.unpack('i', fin.read(4))[0]
print "[ID] BlockSize = ", blksize
fin.read(blksize)
blksize = struct.unpack('i', fin.read(4))[0]
# Mass
blksize = struct.unpack('i', fin.read(4))[0]
print "[Mass] BlockSize = ", blksize
fin.read(blksize)
blksize = struct.unpack('i', fin.read(4))[0]
# U
blksize = struct.unpack('i', fin.read(4))[0]
fin.read(blksize)
# U = []
# for i in range(blksize/4):
#     U.append(struct.unpack("f", fin.read(4))[0])
blksize = struct.unpack('i', fin.read(4))[0]
# Rho
blksize = struct.unpack('i', fin.read(4))[0]
# Rho = []
# for i in range(blksize/4):
#     Rho.append(struct.unpack("f", fin.read(4))[0])
fin.read(blksize)
print "BlockSize = ", struct.unpack('i', fin.read(4))[0]

# NE
blksize = struct.unpack('i', fin.read(4))[0]
fin.read(blksize)
print "BlockSize = ", struct.unpack('i', fin.read(4))[0]
# NH
blksize = struct.unpack('i', fin.read(4))[0]
fin.read(blksize)
print "BlockSize = ", struct.unpack('i', fin.read(4))[0]
# HSML
blksize = struct.unpack('i', fin.read(4))[0]
fin.read(blksize-4)
print "Hsml[Last] = ", struct.unpack('f', fin.read(4))[0]
print "BlockSize = ", struct.unpack('i', fin.read(4))[0]
# SFR
blksize = struct.unpack('i', fin.read(4))[0]
fin.read(blksize-4)
print "Sfr[Last] = ", struct.unpack('f', fin.read(4))[0]
print "BlockSize = ", struct.unpack('i', fin.read(4))[0]
# DELAYT
blksize = struct.unpack('i', fin.read(4))[0]
fin.read(blksize-4)
print "DelayTime[Last] = ", struct.unpack('f', fin.read(4))[0]
print "BlockSize = ", struct.unpack('i', fin.read(4))[0]
# AGE? WARNING: npart[4] > 0 ONLY
blksize = struct.unpack('i', fin.read(4))[0]
print "BlockSize1 = ", blksize
fin.read(blksize-4)
print "Age[Last] = ", struct.unpack('f', fin.read(4))[0]
print "BlockSize2 = ", struct.unpack('i', fin.read(4))[0]
blksize = struct.unpack('i', fin.read(4))[0]
print "BlockSize1 = ", blksize
fin.read(blksize-4)
print "Metal[Last] = ", struct.unpack('f', fin.read(4))[0]
print "BlockSize2 = ", struct.unpack('i', fin.read(4))[0]
