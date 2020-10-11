#adopted from tipsy-tools
CC = gcc
CFLAGS = -O2 -g
#OPT += -DMARK64
OPT += -DGASONLY
#OPT += -DOUTPUT_ALPHA

ifeq (OUTPUT_ALPHA,$(findstring OUTPUT_ALPHA,$(OPT)))
AUX_CONVERTER = aux2ascii_alpha
else
AUX_CONVERTER = aux2ascii
endif

TOOLS = bin2ascii $(AUX_CONVERTER) idnum2ascii $(AW_CONVERTER)

all: $(TOOLS) 

idnum2ascii: idnum2ascii.o
	$(CC) $(CFLAGS) $(OPT) -o idnum2ascii idnum2ascii.o $(LIBS)

bin2ascii: bin2ascii.o
	$(CC) $(CFLAGS) $(OPT) -o bin2ascii bin2ascii.o $(LIBS)

$(AUX_CONVERTER): aux2ascii.o
	$(CC) $(CFLAGS) $(OPT) -o $(AUX_CONVERTER) aux2ascii.o $(LIBS)

bin2ascii.o: bin2ascii.c /usr/include/stdio.h tipsydefs.h
	$(CC) $(CFLAGS) $(OPT) -c bin2ascii.c $(LIBS)

aux2ascii.o: aux2ascii.c /usr/include/stdio.h tipsydefs.h
	$(CC) $(CFLAGS) $(OPT) -c aux2ascii.c $(LIBS)

idnum2ascii.o: idnum2ascii.c /usr/include/stdio.h 
	$(CC) $(CFLAGS) $(OPT) -c idnum2ascii.c $(LIBS)

clean:
	rm -f *.o
	rm -f $(TOOLS)


