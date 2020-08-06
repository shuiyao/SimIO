#HDF5INCL = -I/Users/bob/local/hdf5/include
HDF5INCL = -I/usr/include/hdf5/serial
HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a -lpthread -lz -ldl

OPT += -DH5_USE_16_API

OPT += -DINTERNAL_ENERGY # Convert T (tipsy default) to u (Gadget default)

OPT += -DBOXSIZE=6000 # kpc
OPT += -DOMEGA_0=0.30
OPT += -DOMEGA_L=0.70
OPT += -DHUBBLE_PARAM=0.70

CC= gcc
CLINK=gcc
#CFLAGS= -O2 -g -Wall $(OPT) $(HDF5INCL)
CFLAGS= -g -Wall $(OPT) $(HDF5INCL)
CLIB= -lm $(HDF5LIB)

all: tipsy2hdf5

OBJS=	tipsy2hdf5.o 
tipsy2hdf5:  $(OBJS) Makefile
	 $(CLINK) $(CFLAGS) -o tipsy2hdf5  $(OBJS) $(CLIB)

clean:
	rm -f $(OBJS)
	rm -f tipsy2hdf5

