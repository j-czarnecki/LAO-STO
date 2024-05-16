TARGET = LAO_STO.x
SHELL = /bin/sh
F90 = ifx
F90FLAGS = -O3 -ipo -check all -warn all
LIBS = -llapack -lblas

TEST_TARGET = TEST_LAO_STO.x
POSTPROCESSING_TARGET = POSTPROCESSING_LAO_STO.x
CHERN_TARGET = CHERN_LAO_STO.x

.PHONY: all debug clean test postprocessing chern gnu


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

###### Plotting dispersion relation, calculating DOS and Chern number
POSTPROCESSING_OBJS = 	main_postprocessing.o \
						mod_postprocessing.o \
						mod_hamiltonians.o \
						mod_parameters.o \
						mod_utilities.o \
						mod_writers.o \
						mod_reader.o \
						mod_compute_hamiltonians.o

$(POSTPROCESSING_TARGET): $(POSTPROCESSING_OBJS)
	$(F90) -o $(POSTPROCESSING_TARGET) $(F90FLAGS) $^ $(LIBS)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

all: $(TARGET)

gnu: F90 = gfortran
gnu: F90FLAGS = -O3 -Wall -Wextra -ffree-line-length-none
gnu: $(TARGET)

debug: F90FLAGS = -O0 -g -check all -debug all -warn all -diag-enable sc
debug: $(TARGET)

test: F90FLAGS = -O0 -g -check all -debug all -warn all -diag-enable sc
test: $(TEST_TARGET)

postprocessing: F90FLAGS = -O3
postprocessing: $(POSTPROCESSING_TARGET)

postprocessing_debug: F90 = gfortran
postprocessing_debug: F90FLAGS = -O0 -g -fcheck=array-temps,bounds,do,mem -fsanitize=address -ffree-line-length-none #-check all -warn all #-g flag necessary to generate debug symbols
postprocessing_debug: $(POSTPROCESSING_TARGET)


clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.x
	rm -rf *.i90




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

main_postprocessing.o: 	mod_hamiltonians.o \
						mod_parameters.o \
						mod_utilities.o \
						mod_writers.o \
						mod_reader.o \
						mod_compute_hamiltonians.o \
						mod_postprocessing.o

chern.o: 		mod_hamiltonians.o \
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
							mod_hamiltonians.o \
							mod_writers.o

mod_integrate.o: mod_parameters.o \
				 mod_compute_hamiltonians.o

mod_postprocessing.o: 	mod_hamiltonians.o \
						mod_parameters.o \
						mod_utilities.o \
						mod_writers.o \
						mod_reader.o \
						mod_compute_hamiltonians.o



