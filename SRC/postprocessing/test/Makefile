FC = ifx
PFUNIT_DIR ?= /home/czarnecki/pFUnit/build/installed
GIT_ROOT := $(shell git rev-parse --show-toplevel)
BUILD_DIR := build

.PHONY: all run clean tidy

all: energy_test_suite symmetry_test_suite topology_test_suite tidy

run:
	./energy_test_suite
	./symmetry_test_suite
	./topology_test_suite

clean:
	$(RM) -rf $(BUILD_DIR) energy_test_suite symmetry_test_suite topology_test_suite OutputData/GammaK*

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
FFLAGS += -fpe0
FFLAGS += -fpp -DBAND_BASIS

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += $(LIBS)
FFLAGS += $(INC_DIRS)
FFLAGS += $(LIBS)

energy_test_suite_TESTS := test_energy.pf
energy_test_suite_REGISTRY :=
energy_test_suite_OTHER_SOURCES :=
energy_test_suite_OTHER_LIBRARIES := -L$(SUT_LIB_DIR) -lsut $(LIBS)
energy_test_suite_OTHER_INCS := $(INC_DIRS)
$(eval $(call make_pfunit_test,energy_test_suite))

symmetry_test_suite_TESTS := test_symmetry.pf
symmetry_test_suite_REGISTRY :=
symmetry_test_suite_OTHER_SOURCES :=
symmetry_test_suite_OTHER_LIBRARIES := -L$(SUT_LIB_DIR) -lsut $(LIBS)
symmetry_test_suite_OTHER_INCS := $(INC_DIRS)
$(eval $(call make_pfunit_test,symmetry_test_suite))

topology_test_suite_TESTS := test_topology.pf
topology_test_suite_REGISTRY :=
topology_test_suite_OTHER_SOURCES :=
topology_test_suite_OTHER_LIBRARIES := -L$(SUT_LIB_DIR) -lsut $(LIBS)
topology_test_suite_OTHER_INCS := $(INC_DIRS)
$(eval $(call make_pfunit_test,topology_test_suite))
