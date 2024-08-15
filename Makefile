#
# **************************************************************
# *                Simple C++ Makefile Template                *
# *                                                            *
# * Author: Arash Partow (2003)                                *
# * URL: http://www.partow.net/programming/makefile/index.html *
# *                                                            *
# * Copyright notice:                                          *
# * Free use of this C++ Makefile template is permitted under  *
# * the guidelines and in accordance with the the MIT License  *
# * http://www.opensource.org/licenses/MIT                     *
# *                                                            *
# **************************************************************
#

CXX      := -g++
CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror

GCC_LIB := /usr/local/Compilers/GCC/12.1.0/lib
OpenMPI_LIB := /usr/local/Libraries/OpenMPI/4.1.4/lib
fmt_LIB := /usr/local/Libraries/fmt/9.1.0/lib
LAMMPS_LIB := /usr/local/Software/LAMMPS/OMP_GPU_Kokkos/lib

GCC_INC := /usr/local/Compilers/GCC/12.1.0/include
OpenMPI_INC := /usr/local/Libraries/OpenMPI/4.1.4/include
fmt_INC := /usr/local/Libraries/fmt/9.1.0/include
boost_INC := /usr/local/Libraries/boost/1.82.0/include

BASE_LDFLAGS  := -L${GCC_LIB} -lstdc++ -lm -std=c++17
OPENMPI_LDFLAGS := -L${OpenMPI_LIB} -lmpi -pthread -Wl,-rpath -Wl,${OpenMPI_LIB} -Wl,--enable-new-dtags
FMT_LDFLAGS := -L${fmt_LIB} -lfmt
LAMMPS_LDFLAGS := ${LAMMPS_LIB}/liblammps_OMP_GPU_Kokkos.so
LDFLAGS := ${BASE_LDFLAGS} ${OPENMPI_LDFLAGS} ${FMT_LDFLAGS} ${LAMMPS_LDFLAGS}

BASE_INCLUDE  := -Iinclude/ -I${GCC_INC}
OPENMPI_INCLUDE := -I${OpenMPI_INC}
FMT_INCLUDE := -I${fmt_INC}
BOOST_INCLUDE := -I${boost_INC}
LAMMPS_INCLUDE := -I/home/andrew/Software/Suites/LAMMPS/lammps-27Jun2024/src
INCLUDE := ${BASE_INCLUDE} ${OPENMPI_INCLUDE} ${FMT_INCLUDE} ${BOOST_INCLUDE} ${LAMMPS_INCLUDE}

BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps

TARGET   := btree_chromo

# SRC      :=                      \
#    $(wildcard src/LAMMPS_sys/*.cpp) \
#    $(wildcard src/rep_kinetics/*.cpp) \
#    $(wildcard src/btree/*.cpp) \
#    $(wildcard src/*.cpp)

SRC      := $(wildcard src/*.cpp)

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES \
         := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g -ggdb3
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
	@echo "[*] Include:         ${INCLUDE}     "
