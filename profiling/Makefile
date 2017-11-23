CC = gcc
SDIR = ../src
IDIR = ../include
CFLAGS = -I$(IDIR) `gsl-config --libs --cflags` 
LFLAGS = `gsl-config --libs --cflags` -lm

_DEPS =
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

OBJ = profile.o C_averaging.o C_bias.o C_boostfactors.o C_concentration.o C_deltasigma.o C_density.o C_massfunction.o C_miscentering.o C_xi.o

%.o : $(SDIR)/%.c
	$(CC) -c -o $@ $< $(CFLAGS)

profile : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

clean :
	rm -f *~ *.o profile
