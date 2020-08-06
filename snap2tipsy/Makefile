CC= cc
FC= g77
CLINK=cc
FLINK=g77
CFLAGS= -O2 -Wall -g
FFLAGS= -O
CLIB= -lm 
FLIB= 

DEFS =		 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MATH_H=1 -DHAVE_LIMITS_H=1 -DHAVE_FTRUNCATE=1 -DHAVE_LONGLONG=1 -DHAVE_SHMEM_SERVICES=1 -DHAVE_NET_SERVICES=1 

#DEFS += -DMARK64
DEFS += -DOUTPUTDELAYTIME
#DEFS += -DAGE_BLOCK_FORWARD
DEFS += -DPOT_BLOCK_BACKWARD
#DEFS += -DDEBUG_NOTE
DEFS += -DWRITE_NREC
DEFS += -DOUTPUT_DEBUG
DEFS += -DVELCORR

SOURCES = ${CORE_SOURCES} ${FITSIO_SRC}

OBJECTS = 	${SOURCES:.c=.o}

CORE_OBJECTS = 	${CORE_SOURCES:.c=.o}

snap2tipsy:	snap2tipsy.o snap2tipsy_functions.o load_snapshot.o  ${OBJECTS}
	$(CLINK) $(CFLAGS) $(DEFS) $(OPT) -o snap2tipsy snap2tipsy.o snap2tipsy_functions.o load_snapshot.o $(CLIB) $(LDFLAGS)

snap2tipsy_metal:	snap2tipsy_metal.o snap2tipsy_functions.o load_snapshot.o  ${OBJECTS}
	$(CLINK) $(CFLAGS) $(DEFS) -o snap2tipsy_metal snap2tipsy_metal.o snap2tipsy_functions.o load_snapshot.o $(CLIB) $(LDFLAGS)

bin2ascii:     bin2ascii.o ${OBJECTS}
	$(CLINK) $(CFLAGS) $(DEFS) -o bin2ascii bin2ascii.o $(CLIB) $(LDFLAGS)
snap2tipsy.o:	snap2tipsy.c ${OBJECTS}
	$(CLINK) $(CFLAGS) $(DEFS) -c snap2tipsy.c $(CLIB) $(LDFLAGS)
snap2tipsy_functions.o:	snap2tipsy_functions.c ${OBJECTS}
	$(CLINK) $(CFLAGS) $(DEFS) -c snap2tipsy_functions.c $(CLIB) $(LDFLAGS)
load_snapshot.o:	load_snapshot.c ${OBJECTS}
	$(CLINK) $(CFLAGS) $(DEFS) -c load_snapshot.c $(CLIB) $(LDFLAGS)



clean:
	rm -f *.o snap2tipsy
