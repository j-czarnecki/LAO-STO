TARGET = ./bin/LAO_STO.x
POSTPROCESSING_TARGET = ./bin/POST_LAO_STO.x
SRC_DIR = SRC
OBJ_DIR = OBJ
MOD_DIR = MOD


ifeq ($(READ_OLD), TRUE)
READ_OLD_FLAG = -DREAD_OLD
else
READ_OLD_FLAG =
endif

# the command shell
SHELL = /bin/sh
F90 = ifx
CC = gcc
CXX = g++
LIB_OPENMP = -qopenmp -qmkl
F90FLAGS = -Ofast -fpp -g -ipo $(LIB_OPENMP) -module $(MOD_DIR) $(READ_OLD_FLAG)
LIBS = -llapack -lblas
LIBS_MKL = -I${MKLROOT}/include \
					 -I/opt/intel/mkl/include \
					 -Wl,--start-group \
           ${MKLROOT}/lib/intel64/libmkl_gf_lp64.so \
           ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so \
           ${MKLROOT}/lib/intel64/libmkl_core.so \
           -Wl,--end-group \
           -lpthread -lm -ldl #-lgomp


.PHONY: all  ares_all  ares_post gnu tsan debug clean test post post_debug analyze


####### Full-performance build (default)
OBJS = 	$(OBJ_DIR)/main.o \
				$(OBJ_DIR)/mod_hamiltonians.o \
				$(OBJ_DIR)/mod_parameters.o \
				$(OBJ_DIR)/mod_utilities.o \
				$(OBJ_DIR)/mod_writers.o \
				$(OBJ_DIR)/mod_reader.o \
				$(OBJ_DIR)/mod_broydenV2.o \
				$(OBJ_DIR)/mod_local_integrand.o \
				$(OBJ_DIR)/mod_integrate.o \
				$(OBJ_DIR)/mod_types.o \
				$(OBJ_DIR)/mod_self_consistency.o \
				$(OBJ_DIR)/mod_logger.o

###### Plotting dispersion relation, calculating DOS and Chern number
POSTPROCESSING_OBJS = 	$(OBJ_DIR)/main_postprocessing.o \
						$(OBJ_DIR)/mod_postprocessing.o \
						$(OBJ_DIR)/mod_hamiltonians.o \
						$(OBJ_DIR)/mod_parameters.o \
						$(OBJ_DIR)/mod_utilities.o \
						$(OBJ_DIR)/mod_writers.o \
						$(OBJ_DIR)/mod_reader.o \
						$(OBJ_DIR)/mod_local_integrand.o \
						$(OBJ_DIR)/mod_types.o \
						$(OBJ_DIR)/mod_self_consistency.o \
						$(OBJ_DIR)/mod_logger.o

# Superconductivity calculation target
$(TARGET): $(OBJS)
	$(F90) -o $(TARGET) $(F90FLAGS) $^ $(LIBS) $(LIBS_MKL)

# Postprocessing target
$(POSTPROCESSING_TARGET): $(POSTPROCESSING_OBJS)
	$(F90) -o $(POSTPROCESSING_TARGET) $(F90FLAGS) $^ $(LIBS) $(LIBS_MKL)

# Setting where to find .o and .mod files
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(F90FLAGS) -c $< -o $@

$(OBJ_DIR)/%.s : $(SRC_DIR)/%.f90
	mkdir -p $(OBJ_DIR) $(MOD_DIR)
	$(F90) $(F90FLAGS) -S $< -o $@

all: $(TARGET)

ares_all: LIBS = -lscalapack -lflexiblas
ares_all: $(TARGET)

ares_post: LIBS = -lscalapack -lflexiblas
ares_post: $(POSTPROCESSING_TARGET)

gnu: F90 = gfortran
gnu: F90FLAGS = -O3 -Wall -Wextra -ffree-line-length-none $(LIB_OPENMP) $(READ_OLD_FLAG)
gnu: $(TARGET)

debug: F90FLAGS = -O0 -g -fpp $(LIB_OPENMP) -module $(MOD_DIR) $(READ_OLD_FLAG) -check bounds -debug all #-diag-enable sc
debug: $(TARGET)

#To avoid Thread Sanitizer error about bad memory mapping
#echo 0 | sudo tee /proc/sys/kernel/randomize_va_space
#To include suppressions run as
#TSAN_OPTIONS="suppressions=thread_suppressions.txt:history_size=7" bin/lao_sto_qd.x
tsan: F90FLAGS = -O0 -g -fpp -DDEBUG -fsanitize=thread $(LIB_OPENMP) $(READ_OLD_FLAG)
tsan: $(TARGET)

post: $(POSTPROCESSING_TARGET)

post_debug: F90FLAGS = -O0 -g -fpp -DDEBUG -module $(MOD_DIR) -debug all -fpe0 -fstack-protector -traceback -check bounds,pointers $(LIB_OPENMP) $(READ_OLD_FLAG)
post_debug:	$(POSTPROCESSING_TARGET)

post_tsan: F90FLAGS = -O0 -g -fpp -DDEBUG -fsanitize=thread $(LIB_OPENMP) $(READ_OLD_FLAG)
post_tsan: $(POSTPROCESSING_TARGET)

test:
	mkdir -p $(SRC_DIR)/test/$(MOD_DIR)
	cp SRC/*.f90 SRC/test/
	@export FC="$(F90)" && export CC="$(CC)" && export CXX="$(CXX)" && export FSFLAG=-I && export FCFLAGS="$(F90FLAGS)" && cd $(SRC_DIR)/test && funit
	cd ../../

analyze:
	cd Analyzer && python3 mainAnalyzer.py && cd ..

run_slurm:
	cd Runner && python3 mainRunner.py && cd ..

fit:
	cd Analyzer && python3 DosFitter.py && cd ..

clean_plots:
	cd Plots && rm -rf *.png && cd ..

clean:
	rm -rf $(MOD_DIR)
	rm -rf $(OBJ_DIR)
	rm -f $(TARGET)
	rm -f $(POSTPROCESSING_TARGET)
	rm -f $(SRC_DIR)/test/*.f90
	rm -rf $(SRC_DIR)/test/$(MOD_DIR)
	rm -rf $(SRC_DIR)/test/*.o
	rm -f *.mod
	rm -rf $(SRC_DIR)/*.i90
	rm -rf $(SRC_DIR)/*.mod
	rm -rf *.mod
	@export CC="$(CC)" && export CXX="$(CXX)" && cd $(SRC_DIR)/test && funit --clean && cd ../../



########### Dependencies
$(OBJ_DIR)/main.o:	$(OBJ_DIR)/mod_hamiltonians.o \
										$(OBJ_DIR)/mod_parameters.o \
										$(OBJ_DIR)/mod_utilities.o \
										$(OBJ_DIR)/mod_writers.o \
										$(OBJ_DIR)/mod_reader.o \
										$(OBJ_DIR)/mod_broydenV2.o \
										$(OBJ_DIR)/mod_local_integrand.o \
										$(OBJ_DIR)/mod_integrate.o \
										$(OBJ_DIR)/mod_self_consistency.o \
										$(OBJ_DIR)/mod_logger.o

$(OBJ_DIR)/main_postprocessing.o: 	$(OBJ_DIR)/mod_hamiltonians.o \
																		$(OBJ_DIR)/mod_parameters.o \
																		$(OBJ_DIR)/mod_utilities.o \
																		$(OBJ_DIR)/mod_writers.o \
																		$(OBJ_DIR)/mod_reader.o \
																		$(OBJ_DIR)/mod_local_integrand.o \
																		$(OBJ_DIR)/mod_postprocessing.o \
																		$(OBJ_DIR)/mod_logger.o

$(OBJ_DIR)/chern.o: 		$(OBJ_DIR)/mod_hamiltonians.o \
												$(OBJ_DIR)/mod_parameters.o \
												$(OBJ_DIR)/mod_utilities.o \
												$(OBJ_DIR)/mod_writers.o \
												$(OBJ_DIR)/mod_reader.o \
												$(OBJ_DIR)/mod_local_integrand.o

$(OBJ_DIR)/mod_utilities.o: $(OBJ_DIR)/mod_parameters.o \
														$(OBJ_DIR)/mod_reader.o

$(OBJ_DIR)/mod_parameters.o:

$(OBJ_DIR)/mod_hamiltonians.o:	$(OBJ_DIR)/mod_utilities.o \
								$(OBJ_DIR)/mod_parameters.o \
								$(OBJ_DIR)/mod_reader.o

$(OBJ_DIR)/mod_writers.o: $(OBJ_DIR)/mod_parameters.o \
						  $(OBJ_DIR)/mod_reader.o

$(OBJ_DIR)/mod_reader.o: $(OBJ_DIR)/mod_parameters.o \
						 $(OBJ_DIR)/mod_logger.o

$(OBJ_DIR)/mod_local_integrand.o: 	$(OBJ_DIR)/mod_parameters.o \
										$(OBJ_DIR)/mod_utilities.o \
										$(OBJ_DIR)/mod_hamiltonians.o \
										$(OBJ_DIR)/mod_writers.o

$(OBJ_DIR)/mod_integrate.o: $(OBJ_DIR)/mod_parameters.o \
							$(OBJ_DIR)/mod_local_integrand.o \
							$(OBJ_DIR)/mod_logger.o

$(OBJ_DIR)/mod_postprocessing.o: 	$(OBJ_DIR)/mod_hamiltonians.o \
									$(OBJ_DIR)/mod_parameters.o \
									$(OBJ_DIR)/mod_utilities.o \
									$(OBJ_DIR)/mod_writers.o \
									$(OBJ_DIR)/mod_reader.o \
									$(OBJ_DIR)/mod_local_integrand.o \
									$(OBJ_DIR)/mod_self_consistency.o \
									$(OBJ_DIR)/mod_logger.o

$(OBJ_DIR)/mod_self_consistency.o: $(OBJ_DIR)/mod_parameters.o \
																	 $(OBJ_DIR)/mod_reader.o \
																	 $(OBJ_DIR)/mod_logger.o


$(OBJ_DIR)/mod_logger.o:
