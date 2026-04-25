FC = ifx

TARGET := ./postprocessing_profiling
GIT_ROOT := $(shell git rev-parse --show-toplevel)
SRC_DIR := $(GIT_ROOT)/SRC
OBJ_DIR := OBJ
MOD_DIR := MOD

SRCS := $(SRC_DIR)/postprocessing/src/energy.f90 \
				$(SRC_DIR)/postprocessing/src/symmetry.f90 \
				$(SRC_DIR)/postprocessing/src/topology.f90 \
				$(SRC_DIR)/physical/src/hamiltonians.f90 \
				$(SRC_DIR)/physical/src/utilities.f90 \
				$(SRC_DIR)/physical/src/parameters.f90 \
				$(SRC_DIR)/input_output/src/reader.f90 \
				$(SRC_DIR)/input_output/src/writers.f90 \
				$(SRC_DIR)/input_output/src/logger.f90 \
				$(SRC_DIR)/integrate/src/local_integrand.f90 \
				$(SRC_DIR)/self_consistency/self_consistency.f90 \
				$(SRC_DIR)/types/types.f90 \
				$(SRC_DIR)/physical/src/interaction_factory.f90 \
				$(SRC_DIR)/postprocessing/test/test_profiling.f90

OBJS:= $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

LIBS = -qmkl
FFLAGS += -qmkl
FFLAGS += -fpp -DBAND_BASIS -I$(SRC_DIR)/input_output/src -module $(MOD_DIR)
OPT_FLAGS := -O3 -ipo -qopt-report
FFLAGS += $(OPT_FLAGS)

.PHONY: profile run clean

clean:
	$(RM) -rf *.o *.mod *.a $(OBJ_DIR) $(MOD_DIR) postprocessing_profiling *.optrpt

profile: $(TARGET)

run:
	./postprocessing_profiling

$(TARGET) : $(OBJS)
	$(FC) -o $(TARGET) $(FFLAGS) $^ $(LIBS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

include $(SRC_DIR)/postprocessing/src/deps.mk
