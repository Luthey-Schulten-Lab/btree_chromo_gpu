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

BASE_LDFLAGS  := -L${GCC_LIB} -lstdc++ -lm -std=c++17
OPENMPI_LDFLAGS := -L${OpenMPI_LIB} -lmpi -pthread -Wl,-rpath -Wl,${OpenMPI_LIB} -Wl,--enable-new-dtags
FMT_LDFLAGS := -L${fmt_LIB} -lfmt
LAMMPS_LDFLAGS := ${LAMMPS_LIB}/liblammps_twistable_BD_OMP.so
LDFLAGS := ${BASE_LDFLAGS} ${OPENMPI_LDFLAGS} ${FMT_LDFLAGS} ${LAMMPS_LDFLAGS}

BASE_INCLUDE  := -Iinclude/ -I${GCC_INC}
OPENMPI_INCLUDE := -I${OpenMPI_INC}
FMT_INCLUDE := -I${fmt_INC}
LAMMPS_INCLUDE := -I/home/ben/Software/Suites/LAMMPS/lammps/src
INCLUDE := ${BASE_INCLUDE} ${OPENMPI_INCLUDE} ${FMT_INCLUDE} ${LAMMPS_INCLUDE}

BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps

TARGET   := program

SRC      :=                      \
   $(wildcard src/LAMMPS_sys/*.cpp) \
   $(wildcard src/rep_kinetics/*.cpp) \
   $(wildcard src/btree/*.cpp) \
   $(wildcard src/*.cpp) \

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
