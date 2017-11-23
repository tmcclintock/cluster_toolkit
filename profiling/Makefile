CC = gcc
SDIR = ../src
IDIR = ../include
CFLAGS = -I$(IDIR) `gsl-config --libs --cflags` 
LFLAGS = `gsl-config --libs --cflags` -lm

_DEPS =
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

OBJ = profile.o C_boostfactors.o 

%.o : $(SDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS)

profile : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

clean :
	rm -f *~ *.o profile
