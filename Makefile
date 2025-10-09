# This file is part of LAO-STO.
#
# Copyright (C) 2025 Julian Czarnecki
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# If you use this code for scientific research, please cite:
# J. Czarnecki et. al.,
# "Superconducting gap symmetry of 2DEG at (111)-oriented LaAlO3/SrTiO3 interface",
# arXiv:2508.05075 (2025).
# https://arxiv.org/abs/2508.05075

TARGET = ./bin/LAO_STO.x
POSTPROCESSING_TARGET = ./bin/POST_LAO_STO.x
SRC_DIR = SRC
OBJ_DIR = OBJ
MOD_DIR = MOD

ifeq ($(READ_NO_TRIPLET), TRUE)
READ_OLD_FLAG = -DREAD_NO_TRIPLET
endif
ifeq ($(READ_NO_BAND), TRUE)
READ_OLD_FLAG = -DREAD_NO_BAND
endif

# the command shell
SHELL = /bin/sh
F90 = ifx
CC = gcc
CXX = g++

LIB_OPENMP = -qopenmp -qmkl
F90FLAGS = -Ofast -ipo -g -fpp -ipo -I$(SRC_DIR)/input_output $(LIB_OPENMP) -module $(MOD_DIR) $(READ_OLD_FLAG) -diag-disable 5268,7025 -stand f2018
F90_DEBUG_FLAGS = -O0 -g -fpp -I$(SRC_DIR)/input_output -DDEBUG -module $(MOD_DIR) -debug all -fpe0 -fstack-protector -traceback -check all -ftrapuv -heap-arrays $(LIB_OPENMP) $(READ_OLD_FLAG)
LIBS = -llapack -lblas
LIBS_MKL = -I${MKLROOT}/include \
					 -I/opt/intel/mkl/include \
					 -Wl,--start-group \
           ${MKLROOT}/lib/intel64/libmkl_gf_lp64.so \
           ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so \
           ${MKLROOT}/lib/intel64/libmkl_core.so \
           -Wl,--end-group \
           -lpthread -lm -ldl #-lgomp

# --- Automatically find all source files recursively ---
SRC_FILES_ALL := $(shell find $(SRC_DIR) -name '*.f90')

# --- Exclude the two main programs from the common source set ---
SRC_COMMON := $(filter-out $(SRC_DIR)/main/main.f90 $(SRC_DIR)/main_post/main_postprocessing.f90, $(SRC_FILES_ALL))

# --- Define two build sets ---
SRC_FILES_MAIN := $(SRC_COMMON) $(SRC_DIR)/main/main.f90
SRC_FILES_POST := $(SRC_COMMON) $(SRC_DIR)/main_post/main_postprocessing.f90

# --- Define corresponding object files ---
OBJS_MAIN := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES_MAIN))
OBJS_POST := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES_POST))

.PHONY: all  ares_all  ares_post gnu tsan debug clean test post post_debug analyze

# Superconductivity calculation target
$(TARGET): $(OBJS_MAIN)
	$(F90) -o $(TARGET) $(F90FLAGS) $^ $(LIBS) $(LIBS_MKL)

# Postprocessing target
$(POSTPROCESSING_TARGET): $(OBJS_POST)
	$(F90) -o $(POSTPROCESSING_TARGET) $(F90FLAGS) $^ $(LIBS) $(LIBS_MKL)

# Setting where to find .o and .mod files
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(F90) $(F90FLAGS) -c $< -o $@

$(OBJ_DIR)/%.s : $(SRC_DIR)/%.f90
	mkdir -p $(dir $@) $(MOD_DIR)
	$(F90) $(F90FLAGS) -S $< -o $@

all: $(TARGET)

ares_all: LIBS = -lscalapack -lflexiblas
ares_all: $(TARGET)

ares_post: LIBS = -lscalapack -lflexiblas
ares_post: $(POSTPROCESSING_TARGET)

gnu: F90 = gfortran
gnu: LIB_OPENMP = -fopenmp -mkl
gnu: LIBS = -lscalapack -lflexiblas
gnu: F90FLAGS = -O3 -g -cpp -Wall -Wextra -ffree-line-length-none $(LIB_OPENMP) $(READ_OLD_FLAG) -J$(MOD_DIR)
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
	mkdir -p $(SRC_DIR)/physical/test/$(MOD_DIR)
	find $(SRC_DIR) -type f -name "*.f90" ! -path "$(SRC_DIR)/physical/test/*" -exec cp {} $(SRC_DIR)/physical/test/ \;
	@export FC="$(F90)" && export CC="$(CC)" && export CXX="$(CXX)" && export FSFLAG=-I && export FCFLAGS="$(F90_DEBUG_FLAGS)" && cd $(SRC_DIR)/physical/test && funit
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
	rm -f $(SRC_DIR)/**/test/*.f90
	rm -rf $(SRC_DIR)/**/test/$(MOD_DIR)
	rm -rf $(SRC_DIR)/**/test/*.o
	rm -f *.mod
	rm -rf $(SRC_DIR)/**/*.i90
	rm -rf $(SRC_DIR)/**/*.mod
	rm -rf *.mod
	@export CC="$(CC)" && export CXX="$(CXX)" && cd $(SRC_DIR)/physical/test && funit --clean && cd ../../



########### Dependencies
$(OBJ_DIR)/main/main.o: $(OBJ_DIR)/physical/hamiltonians.o \
									 $(OBJ_DIR)/physical/parameters.o \
									 $(OBJ_DIR)/physical/utilities.o \
									 $(OBJ_DIR)/input_output/writers.o \
									 $(OBJ_DIR)/input_output/reader.o \
									 $(OBJ_DIR)/self_consistency/broydenV2.o \
									 $(OBJ_DIR)/integrate/local_integrand.o \
									 $(OBJ_DIR)/integrate/integrate.o \
									 $(OBJ_DIR)/self_consistency/self_consistency.o \
									 $(OBJ_DIR)/input_output/logger.o \
									 $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/main_postprocessing/main_postprocessing.o: $(OBJ_DIR)/physical/hamiltonians.o \
																  $(OBJ_DIR)/physical/parameters.o \
																  $(OBJ_DIR)/physical/utilities.o \
																  $(OBJ_DIR)/input_output/writers.o \
																  $(OBJ_DIR)/input_output/reader.o \
																  $(OBJ_DIR)/integrate/local_integrand.o \
																  $(OBJ_DIR)/postprocessing/postprocessing.o \
																  $(OBJ_DIR)/input_output/logger.o

$(OBJ_DIR)/chern.o: $(OBJ_DIR)/physical/hamiltonians.o \
									  $(OBJ_DIR)/physical/parameters.o \
									  $(OBJ_DIR)/physical/utilities.o \
									  $(OBJ_DIR)/input_output/writers.o \
									  $(OBJ_DIR)/input_output/reader.o \
									  $(OBJ_DIR)/integrate/local_integrand.o

$(OBJ_DIR)/physical/utilities.o: $(OBJ_DIR)/physical/parameters.o \
												$(OBJ_DIR)/input_output/reader.o

$(OBJ_DIR)/physical/parameters.o:

$(OBJ_DIR)/types/types.o: $(OBJ_DIR)/physical/parameters.o

$(OBJ_DIR)/physical/hamiltonians.o: $(OBJ_DIR)/physical/utilities.o \
								           $(OBJ_DIR)/physical/parameters.o \
								           $(OBJ_DIR)/input_output/reader.o \
								           $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/input_output/writers.o: $(OBJ_DIR)/physical/parameters.o \
						          $(OBJ_DIR)/input_output/reader.o \
						          $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/input_output/reader.o: $(OBJ_DIR)/physical/parameters.o \
						         $(OBJ_DIR)/input_output/logger.o \
						         $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/integrate/local_integrand.o: $(OBJ_DIR)/physical/parameters.o \
										          $(OBJ_DIR)/physical/utilities.o \
										          $(OBJ_DIR)/physical/hamiltonians.o \
										          $(OBJ_DIR)/input_output/writers.o \
										          $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/integrate/integrate.o: $(OBJ_DIR)/physical/parameters.o \
							          $(OBJ_DIR)/integrate/local_integrand.o \
							          $(OBJ_DIR)/input_output/logger.o \
							          $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/postprocessing/postprocessing.o: 	$(OBJ_DIR)/physical/hamiltonians.o \
									            $(OBJ_DIR)/physical/parameters.o \
									            $(OBJ_DIR)/physical/utilities.o \
									            $(OBJ_DIR)/input_output/writers.o \
									            $(OBJ_DIR)/input_output/reader.o \
									            $(OBJ_DIR)/integrate/local_integrand.o \
									            $(OBJ_DIR)/self_consistency/self_consistency.o \
									            $(OBJ_DIR)/input_output/logger.o \
									            $(OBJ_DIR)/types/types.o

$(OBJ_DIR)/self_consistency/self_consistency.o: $(OBJ_DIR)/physical/parameters.o \
															 $(OBJ_DIR)/input_output/reader.o \
															 $(OBJ_DIR)/input_output/logger.o \
															 $(OBJ_DIR)/types/types.o \
															 $(OBJ_DIR)/input_output/writers.o

$(OBJ_DIR)/input_output/logger.o:
