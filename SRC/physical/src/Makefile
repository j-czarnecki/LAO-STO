FC = ifx
AR = ar
PFUNIT_DIR ?= /home/czarnecki/pFUnit/build/installed

GIT_ROOT := $(shell git rev-parse --show-toplevel)
SRC_DIR := $(GIT_ROOT)/SRC
OBJ_DIR := OBJ
MOD_DIR := MOD

SRCS := $(SRC_DIR)/physical/src/hamiltonians.f90 \
					$(SRC_DIR)/physical/src/utilities.f90 \
					$(SRC_DIR)/physical/src/parameters.f90 \
					$(SRC_DIR)/physical/src/interaction_factory.f90 \
					$(SRC_DIR)/input_output/src/reader.f90 \
					$(SRC_DIR)/input_output/src/writers.f90 \
					$(SRC_DIR)/input_output/src/logger.f90 \
					$(SRC_DIR)/types/types.f90

OBJS:= $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

LIBS = -qmkl
FFLAGS += -qmkl
FFLAGS += -fpp -I$(SRC_DIR)/input_output/src -module $(MOD_DIR)
FFLAGS += -O0 -g
FFLAGS += -traceback
FFLAGS += -check all
#FFLAGS += -warn all -warn errors
FFLAGS += -fpe0

.PHONY: all clean libsut.a

all: libsut.a

libsut.a: $(OBJS)
	$(AR) -r $@ $?

clean:
	$(RM) -rf *.o *.mod *.a $(OBJ_DIR) $(MOD_DIR)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

include deps.mk
