# To get the HDF5 dependencies, use h5cc -v

# EAGLE
HDF5INCL = -I/home/shuiyao/hdf5/include
HDF5LIB = -L/home/shuiyao/hdf5/lib -lhdf5 -lz

#MIES
HDF5INCL = -I/home/shuiyao/include
HDF5LIB = -L/home/shuiyao/lib -lhdf5 -lz


OPT += -DH5_USE_16_API
OPT += -DNMETALS=11 # How many species are tracked in the HDF5 "Metallicity"  block? If not set correctly, code returns segmentation fault after the output "Reading GAS for file #%d"

CC= gcc
CLINK=gcc
CFLAGS= -O2 -g -Wall $(OPT) $(HDF5INCL)
#CFLAGS= -g -Wall $(OPT) $(HDF5INCL)
CLIB= -lm $(HDF5LIB)

all: hdf52gadget

OBJS=	hdf52gadget.o 
hdf52gadget:  $(OBJS) Makefile
	 $(CLINK) $(CFLAGS) -o hdf52gadget  $(OBJS) $(CLIB)

clean:
	rm -f $(OBJS)
	rm -f hdf52gadget

