#HDF5INCL = -I/Users/bob/local/hdf5/include
# HDF5INCL = -I/usr/include/hdf5/serial
# HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a -lpthread -lz -ldl

HDF5INCL = -I/home/shuiyao/hdf5/include
HDF5LIB = -L/home/shuiyao/hdf5/lib -lhdf5 -lz

OPT += -DH5_USE_16_API
OPT += -DNMETALS=11 # How many species are tracked in the HDF5 "Metallicity"  block? If not set correctly, code returns segmentation fault after the output "Reading GAS for file #%d"
OPT += -DPHEW_EXTRA_OUTPUT

CC= gcc
CLINK=gcc
#CFLAGS= -O2 -g -Wall $(OPT) $(HDF5INCL)
CFLAGS= -g -Wall $(OPT) $(HDF5INCL)
CLIB= -lm $(HDF5LIB)

all: hdf52ascii

OBJS=	hdf52ascii.o 
hdf52ascii:  $(OBJS) Makefile
	 $(CLINK) $(CFLAGS) -o hdf52ascii  $(OBJS) $(CLIB)

clean:
	rm -f $(OBJS)
	rm -f hdf52ascii

