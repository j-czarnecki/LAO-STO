FC = ifx

TARGET := ./hamiltonians_profiling
GIT_ROOT := $(shell git rev-parse --show-toplevel)
SRC_DIR := $(GIT_ROOT)/SRC
OBJ_DIR := OBJ
MOD_DIR := MOD

SRCS := $(SRC_DIR)/physical/src/hamiltonians.f90 \
				$(SRC_DIR)/physical/src/utilities.f90 \
				$(SRC_DIR)/physical/src/parameters.f90 \
				$(SRC_DIR)/input_output/reader.f90 \
				$(SRC_DIR)/input_output/writers.f90 \
				$(SRC_DIR)/input_output/logger.f90 \
				$(SRC_DIR)/types/types.f90 \
				$(SRC_DIR)/physical/test/test_profiling.f90

OBJS:= $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRCS))

LIBS = -llapack -lblas
FFLAGS += -fpp -I$(SRC_DIR)/input_output -module $(MOD_DIR)
OPT_FLAGS := -O3 -ipo -qopt-report
FFLAGS += $(OPT_FLAGS)

.PHONY: profile run clean

clean:
	$(RM) -rf *.o *.mod *.a $(OBJ_DIR) $(MOD_DIR) hamiltonians_profiling *.optrpt

profile: $(TARGET)

run:
	./hamiltonians_profiling

$(TARGET) : $(OBJS)
	$(FC) -o $(TARGET) $(FFLAGS) $^ $(LIBS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

include $(SRC_DIR)/physical/src/deps.mk