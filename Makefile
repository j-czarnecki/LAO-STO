TARGET = LAO_STO.x
SHELL = /bin/sh
F90 = ifort
F90FLAGS = -O3 -ipo
LIBS = -llapack -lblas

TEST_TARGET = TEST_LAO_STO.x
DISPERSION_TARGET = DISPERSION_LAO_STO.x

.PHONY: all debug clean test dispersion gnu


####### Full-performance build (default)
OBJS = 	main.o \
		mod_hamiltonians.o \
		mod_parameters.o \
		mod_utilities.o \
		mod_writers.o \
		mod_reader.o \
		mod_broydenV2.o \
		mod_compute_hamiltonians.o \
		mod_integrate.o

$(TARGET): $(OBJS)
	$(F90) -o $(TARGET) $(F90FLAGS) $^ $(LIBS)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@


###### Running tests
TEST_OBJS = tests.o \
			mod_hamiltonians.o \
			mod_parameters.o \
			mod_utilities.o \
			mod_writers.o \
			mod_reader.o \
			mod_broydenV2.o \
			mod_compute_hamiltonians.o \
			mod_integrate.o

$(TEST_TARGET): $(TEST_OBJS)
	$(F90) -o $(TEST_TARGET) $(F90FLAGS) $^ $(LIBS)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

###### Plotting dispersion relation, calculating DOS
DISPERSION_OBJS = 	dispersion.o \
					mod_hamiltonians.o \
					mod_parameters.o \
					mod_utilities.o \
					mod_writers.o \
					mod_reader.o \
					mod_compute_hamiltonians.o

$(DISPERSION_TARGET): $(DISPERSION_OBJS)
	$(F90) -o $(DISPERSION_TARGET) $(F90FLAGS) $^ $(LIBS)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

all: $(TARGET)

gnu: F90 = gfortran
gnu: F90FLAGS = -O3 -Wall -Wextra -ffree-line-length-none
gnu: $(TARGET)

debug: F90FLAGS = -O0 -g -check bounds -debug all
debug: $(TARGET)

test: F90FLAGS = -O0 -g -check bounds -debug all
test: $(TEST_TARGET)

dispersion: $(DISPERSION_TARGET)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.x




########### Dependencies

main.o:	mod_hamiltonians.o \
		mod_parameters.o \
		mod_utilities.o \
		mod_writers.o \
		mod_reader.o \
		mod_broydenV2.o \
		mod_compute_hamiltonians.o \
		mod_integrate.o

tests.o:	mod_hamiltonians.o \
			mod_parameters.o \
			mod_writers.o

dispersion.o: 	mod_hamiltonians.o \
				mod_parameters.o \
				mod_utilities.o \
				mod_writers.o \
				mod_reader.o \
				mod_compute_hamiltonians.o

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



