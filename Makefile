TARGET = LAO_STO.x
SHELL = /bin/sh
F90 = ifort
F90FLAGS = -O3 -ipo
LIBS = -llapack -lblas

OBJS = 	main.o \
		mod_hamiltonians.o \
		mod_parameters.o \
		mod_utilities.o \
		mod_writers.o \
		mod_reader.o \
		mod_broydenV2.o \
		mod_compute_hamiltonians.o \
		mod_integrate.o

.PHONY: all debug clean

$(TARGET): $(OBJS)
	$(F90) -o $(TARGET) $(F90FLAGS) $^ $(LIBS)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

all: $(TARGET)

gnu: F90 = gfortran
gnu: F90FLAGS = -O3 -Wall -Wextra -ffree-line-length-none -DTRACE_DEBUG
gnu: $(TARGET)

debug: F90FLAGS = -O0 -g -check bounds -debug all -DTRACE_DEBUG
debug: $(TARGET)

clean:
	rm -f $(OBJS)
	rm -f *.mod
	rm -f $(TARGET)

main.o:	mod_hamiltonians.o \
		mod_parameters.o \
		mod_utilities.o \
		mod_writers.o \
		mod_reader.o \
		mod_broydenV2.o \
		mod_compute_hamiltonians.o \
		mod_integrate.o

mod_utilities.o: mod_parameters.o \
				 mod_reader.o

mod_parameters.o:

mod_hamiltonians.o:	mod_utilities.o \
					mod_parameters.o \
					mod_reader.o

mod_writers.o: mod_parameters.o

mod_reader.o: mod_parameters.o

mod_compute_hamiltonians.o: mod_parameters.o \
							mod_utilities.o \
							mod_hamiltonians.o 

mod_integrate.o: mod_parameters.o \
				 mod_compute_hamiltonians.o



