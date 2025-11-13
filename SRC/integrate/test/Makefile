FC = ifx
PFUNIT_DIR ?= /home/czarnecki/pFUnit/build/installed
GIT_ROOT := $(shell git rev-parse --show-toplevel)
BUILD_DIR := build

.PHONY: all run clean tidy

all: local_integrand_test_suite tidy

run:
	./local_integrand_test_suite

clean:
	$(RM) -rf $(BUILD_DIR) local_integrand_test_suite OutputData/*.dat

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
LIBS := -llapack -lblas

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += $(LIBS)
FFLAGS += $(INC_DIRS)
FFLAGS += $(LIBS)

local_integrand_test_suite_TESTS := test_local_integrand.pf
local_integrand_test_suite_REGISTRY :=
local_integrand_test_suite_OTHER_SOURCES :=
local_integrand_test_suite_OTHER_LIBRARIES := -L$(SUT_LIB_DIR) -lsut $(LIBS)
local_integrand_test_suite_OTHER_INCS := $(INC_DIRS)

$(eval $(call make_pfunit_test,local_integrand_test_suite))