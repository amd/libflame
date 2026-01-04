###############################################################################
# Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
# 
# Makefile for libflame test suite driver.
#

CC             := gcc
CLINKER        := $(CC)
RM_RF          := rm -rf
OMP_LIB        := -fopenmp
CFLAGS         := -std=c11 -D_GNU_SOURCE -O2 -Wall -Wno-unused-function -Wno-comment ${OMP_LIB}
CXXFLAGS       :=
LDFLAGS        := -lstdc++ -lm -lpthread $(OMP_LIB)

ifndef $(ENABLE_CPP_TEST)
# Flag to enable CPP test
ENABLE_CPP_TEST := 1
endif

ifndef $(ENABLE_AOCL_EXTENSION_APIS)
# Flag to enable testing for extension LAPACK APIs not present in 
# the standard LAPACK library
ENABLE_AOCL_EXTENSION_APIS := 1
endif

#
# --- LIBFLAME \ LAPACK library ------------------------------------------------
#
LIB_PATH          := ../../lib
LIBFLAME          := $(LIB_PATH)/libflame.a
LAPACK_INC_PATH   := ../../include

#
# --- BLAS library -------------------------------------------------------------
#
# A BLAS library must be given in order to run the libflame test suite.
LIBBLAS           :=
BLAS_HEADER_PATH  :=

#
# ----- AOCL-Utils library path
#
LIBAOCLUTILS_LIBRARY_PATH :=

#
# --- Uncomment five lines below and set MKLROOT & LDFLAGS ---------------------
#     appropriately to link to MKL libs
#     To enable MKL counterparts of AOCL extension APIs, set 
#     ENABLE_AOCL_EXTENSION_APIS to 2 and define API prototypes in 
#     validate_src/test_prototype.h
#
#LIBFLAME       :=
#LIBBLAS        :=
#MKLROOT        :=
#LDFLAGS        := -Wl,--start-group ${MKLROOT}/intel64/libmkl_intel_lp64.a \
                   ${MKLROOT}/intel64/libmkl_gnu_thread.a \
                   ${MKLROOT}/intel64/libmkl_core.a -Wl,--end-group \
                   -lpthread -lm -fopenmp -ldl -lstdc++
#ENABLE_AOCL_EXTENSION_APIS := 0

#
# --- ILP64 --- Set to 1 to enable ILP64 ---------------------------------------
#
ifndef $(ILP64)
ILP64 = 0
endif

ifeq ($(ILP64),1)
ifeq ($(ENABLE_CPP_TEST), 1)
    CXXFLAGS    += -DFLA_ENABLE_ILP64=1
endif
    CFLAGS      += -DFLA_ENABLE_ILP64=1
endif

#
# --- Flag to enable uppercase/lowercase API calling conventions --------------
#     Set API_CALLING_CONVENTION = upper  for uppercase without underscore
#                                = upper_ for uppercase with underscore
#                                = lower  for lowercase without underscore
#                                = default for default case (lower_)
#     NOTE: LAPACK default is lowercase with underscore
#
API_CALLING_CONVENTION := default

#
# -- API_CALLING_CONVENTION set to "upper" enables uppercase w/o underscore ----
#
ifeq ($(API_CALLING_CONVENTION),upper)
ifeq ($(ENABLE_CPP_TEST), 1)
    CXXFLAGS    += -DUPPER=1
endif
    CFLAGS      += -DUPPER=1
endif

#
# -- API_CALLING_CONVENTION set to "upper_" enables uppercase w/ underscore ----
#
ifeq ($(API_CALLING_CONVENTION),upper_)
ifeq ($(ENABLE_CPP_TEST), 1)
    CXXFLAGS    += -DUPPER_=1
endif
    CFLAGS      += -DUPPER_=1
endif

#
# -- API_CALLING_CONVENTION set to "lower" enables lowercase w/o underscore ----
#
ifeq ($(API_CALLING_CONVENTION),lower)
ifeq ($(ENABLE_CPP_TEST), 1)
    CXXFLAGS    += -DLOWER=1
endif
    CFLAGS      += -DLOWER=1
endif

# --- enable dynamic memory  unaligned  addres logic --------------------------
ifeq ($(MEM_UNALN),1)
ifeq ($(ENABLE_CPP_TEST), 1)
	CXXFLAGS   += -DFLA_MEM_UNALIGNED=1
endif
    CFLAGS     += -DFLA_MEM_UNALIGNED=1
endif

# --- Enable AOCL Extension APIs  ---------------------------------------------
ifneq ($(ENABLE_AOCL_EXTENSION_APIS), 0)
ifeq ($(ENABLE_CPP_TEST), 1)
    CXXFLAGS   += -DENABLE_AOCL_EXTENSION_APIS=$(ENABLE_AOCL_EXTENSION_APIS)
endif
    CFLAGS     += -DENABLE_AOCL_EXTENSION_APIS=$(ENABLE_AOCL_EXTENSION_APIS)
endif

#
# --- General build definitions ------------------------------------------------
#
TEST_SRC_PATH  := src
TEST_OBJ_PATH  := src/obj
TEST_LAPACKE_INVOKE_PATH := src/lapacke
TEST_LAPACKE_OBJ_PATH := src/obj/lapacke
TEST_VALIDATE_SRC_PATH := validate_src
TEST_VALIDATE_OBJ_PATH := validate_src/obj

ifeq ($(ENABLE_CPP_TEST), 1)
CXX                 := g++
CXXFLAGS            += -std=c++17 -Wall ${OMP_LIB}
TEST_CPP_PATH       := src/cpp
TEST_CPP_OBJ_PATH   := src/obj/cpp

CFLAGS         += -DENABLE_CPP_TEST=1 -I$(BLAS_HEADER_PATH) -I$(LAPACK_INC_PATH) -I$(TEST_SRC_PATH) -I$(TEST_VALIDATE_SRC_PATH) -I$(TEST_CPP_PATH) -I$(TEST_LAPACKE_INVOKE_PATH)
CXXFLAGS       += -DENABLE_CPP_TEST=1 -I$(BLAS_HEADER_PATH) -I$(LAPACK_INC_PATH) -I$(TEST_SRC_PATH) -I$(TEST_VALIDATE_SRC_PATH) -I$(TEST_CPP_PATH) -I$(TEST_LAPACKE_INVOKE_PATH)

else
CFLAGS         += -I$(BLAS_HEADER_PATH) -I$(LAPACK_INC_PATH) -I$(TEST_SRC_PATH) -I$(TEST_VALIDATE_SRC_PATH)
endif

FNAME          := lapack

TEST_VALIDATE_OBJS  := $(patsubst $(TEST_VALIDATE_SRC_PATH)/%.c, \
                       $(TEST_VALIDATE_OBJ_PATH)/%.o, \
                       $(wildcard $(TEST_VALIDATE_SRC_PATH)/*.c))
TEST_LAPACKE_OBJS   := $(patsubst $(TEST_LAPACKE_INVOKE_PATH)/%.c, \
                       $(TEST_LAPACKE_OBJ_PATH)/%.o, \
                       $(wildcard $(TEST_LAPACKE_INVOKE_PATH)/*.c))
ifeq ($(ENABLE_CPP_TEST), 1)
TEST_CPP_OBJS       := $(patsubst $(TEST_CPP_PATH)/%.cc, \
                       $(TEST_CPP_OBJ_PATH)/%.o, \
                       $(wildcard $(TEST_CPP_PATH)/*.cc))
endif
TEST_SRC_OBJS       := $(patsubst $(TEST_SRC_PATH)/%.c, \
                       $(TEST_OBJ_PATH)/%.o, \
                       $(wildcard $(TEST_SRC_PATH)/*.c))

TEST_BIN       := test_$(FNAME).x

$(TEST_VALIDATE_OBJ_PATH)/%.o: $(TEST_VALIDATE_SRC_PATH)/%.c
	mkdir -p $(TEST_VALIDATE_OBJ_PATH)
	$(CC) $(CFLAGS) -c $< -o $@

ifeq ($(ENABLE_CPP_TEST), 1)
$(TEST_CPP_OBJ_PATH)/%.o: $(TEST_CPP_PATH)/%.cc
	mkdir -p $(TEST_CPP_OBJ_PATH)
	$(CXX) $(CXXFLAGS) -c $< -o $@
endif

$(TEST_OBJ_PATH)/%.o: $(TEST_SRC_PATH)/%.c
	mkdir -p $(TEST_OBJ_PATH)
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_LAPACKE_OBJ_PATH)/%.o: $(TEST_LAPACKE_INVOKE_PATH)/%.c
	mkdir -p $(TEST_LAPACKE_OBJ_PATH)
	$(CC) $(CFLAGS) -c $< -o $@

ifeq ($(ENABLE_CPP_TEST), 1)
test_$(FNAME): $(TEST_VALIDATE_OBJS) $(TEST_CPP_OBJS) $(TEST_SRC_OBJS) $(TEST_LAPACKE_OBJS)
	$(CLINKER) $(TEST_VALIDATE_OBJS) $(TEST_CPP_OBJS) $(TEST_SRC_OBJS) $(TEST_LAPACKE_OBJS) $(LIBFLAME) $(LIBBLAS) $(LIBAOCLUTILS_LIBRARY_PATH) \
	$(LDFLAGS) -o $(TEST_BIN)
else
test_$(FNAME): $(TEST_VALIDATE_OBJS) $(TEST_SRC_OBJS) $(TEST_LAPACKE_OBJS)
	$(CLINKER) $(TEST_VALIDATE_OBJS) $(TEST_SRC_OBJS) $(TEST_LAPACKE_OBJS) $(LIBFLAME) $(LIBBLAS) $(LIBAOCLUTILS_LIBRARY_PATH) \
	$(LDFLAGS) -o $(TEST_BIN)
endif

clean:
	$(RM_RF) $(TEST_VALIDATE_OBJ_PATH) $(TEST_OBJ_PATH) $(TEST_BIN)
