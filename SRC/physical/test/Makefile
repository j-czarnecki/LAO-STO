FC = ifx
PFUNIT_DIR ?= /home/czarnecki/pFUnit/build/installed
GIT_ROOT := $(shell git rev-parse --show-toplevel)
BUILD_DIR := build

.PHONY: all run clean tidy

all: hamiltonians_test_suite utilities_test_suite tidy

run:
	./hamiltonians_test_suite
	./utilities_test_suite

clean:
	$(RM) -rf $(BUILD_DIR) hamiltonians_test_suite utilities_test_suite OutputData/*.dat

tidy:
	@echo "Tidying up build artifacts..."
	@mkdir -p $(BUILD_DIR)
	@mv -f *.o $(BUILD_DIR) 2>/dev/null || true
	@mv -f *.mod $(BUILD_DIR) 2>/dev/null || true
	@mv -f *.F90 $(BUILD_DIR) 2>/dev/null || true
	@mv -f *.inc $(BUILD_DIR) 2>/dev/null || true
	@echo "→ Object, module, and generated files moved to build/."

%.o : %.f90
	$(FC) -c $(FFLAGS) $< $(LIBS)

%.o : %.F90
	$(FC) -c $(FFLAGS) $< $(LIBS)

LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

INC_DIRS := -I../src/OBJ -I../src/MOD
SUT_LIB_DIR := ../src
LIBS := -qmkl
FFLAGS += -qmkl
FFLAGS += -O0 -g
FFLAGS += -traceback
FFLAGS += -check all
#FFLAGS += -warn all -warn errors
FFLAGS += -fpe0

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += $(LIBS)
FFLAGS += $(INC_DIRS)
FFLAGS += $(LIBS)

hamiltonians_test_suite_TESTS := test_hamiltonians.pf
hamiltonians_test_suite_REGISTRY :=
hamiltonians_test_suite_OTHER_SOURCES :=
hamiltonians_test_suite_OTHER_LIBRARIES := -L$(SUT_LIB_DIR) -lsut $(LIBS)
hamiltonians_test_suite_OTHER_INCS := $(INC_DIRS)
$(eval $(call make_pfunit_test,hamiltonians_test_suite))

utilities_test_suite_TESTS := test_utilities.pf
utilities_test_suite_REGISTRY :=
utilities_test_suite_OTHER_SOURCES :=
utilities_test_suite_OTHER_LIBRARIES := -L$(SUT_LIB_DIR) -lsut $(LIBS)
utilities_test_suite_OTHER_INCS := $(INC_DIRS)
$(eval $(call make_pfunit_test,utilities_test_suite))



