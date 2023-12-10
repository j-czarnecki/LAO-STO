TARGET = LAO_STO.x
SHELL = /bin/sh
F90 = ifort
F90FLAGS = -O3
LIBS = -llapack

OBJS = 	main.o \
		mod_hamiltonians.o \
		mod_parameters.o \
		mod_utilities.o \
		mod_writers.o

$(TARGET): $(OBJS)
	$(F90) -o $(TARGET) $(F90FLAGS) $^ $(LIBS)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)
	rm -f *.mod
	rm -f $(TARGET)


main.o:	mod_hamiltonians.o \
		mod_parameters.o \
		mod_utilities.o \
		mod_writers.o

mod_utilities.o: mod_parameters.o

mod_parameters.o:

mod_hamiltonians.o:	mod_utilities.o \
					mod_parameters.o 

mod_writers.o: mod_parameters.o



