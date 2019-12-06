# This is a makefile.

PROG = generate_mass_spectrum

OBJ = main.o

CC = g++ -std=c++11 -O3 -g -Wall

LFLAG = -lgsl -lgslcblas 

CFLAG = -c

OFLAG = -o

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LIBS) $(LFLAG)

main.o: oligomer.h vector3d.h bead.h subunit.h


clean:
	rm -f *.o
