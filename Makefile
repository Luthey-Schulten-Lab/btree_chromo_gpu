CXX      := -g++
CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror

# CUDA toolchain for the fused CG minimizer and fused BD (src/*.cu). CUDA_ARCHS = one or more compute capabilities:
#   make CUDA_ARCHS=100 (B200, default) | 90 (H100) | 80 (A100) | "80 90 100" (one binary for all of them; sm_70 and newer)
NVCC       := nvcc
CUDA_ARCHS ?= 100
NVCC_ARCH  := $(foreach a,$(subst sm_,,$(CUDA_ARCHS)),-gencode arch=compute_$(a),code=sm_$(a))
NVCC_FLAGS := -O3 $(NVCC_ARCH) -rdc=true --std=c++17

CUDA_LIB := /usr/local/cuda/lib64
CUDA_LDFLAGS := -L${CUDA_LIB} -lcuda -lcudart -Wl,-rpath -Wl,${CUDA_LIB}

GCC_LIB := /usr/local/lib64
GCC_INC := /usr/local/include

OpenMPI_LIB := /usr/local/lib
OpenMPI_INC := /usr/local/include

LAMMPS_LIB := /usr/local/lib
LAMMPS_INC := /usr/local/include

# LAMMPS headers: default matches protein_science container layout.
# Override on host, e.g.:
#   LAMMPS_SRC=/path/to/lammps/src make BUILD=./build_fused release
LAMMPS_SRC ?= /Software/LAMMPS/lammps/src

BASE_LDFLAGS  := -L${GCC_LIB} -lstdc++ -lm -std=c++17
OPENMPI_LDFLAGS := -L${OpenMPI_LIB} -lmpi -lmpi_cxx -pthread -Wl,-rpath -Wl,${OpenMPI_LIB} -Wl,--enable-new-dtags
FMT_LDFLAGS := -lfmt
LAMMPS_LDFLAGS := ${LAMMPS_LIB}/liblammps_OMP_GPU_Kokkos.so
# export the executable's LAMMPS_NS::Special::build so liblammps' read_data binds to it
LDFLAGS := -Wl,--export-dynamic-symbol=_ZN9LAMMPS_NS7Special5buildEv -Wl,--export-dynamic-symbol=_ZN9LAMMPS_NS4Atom10data_bondsEiPcPiiiiS2_ -Wl,--export-dynamic-symbol=_ZN9LAMMPS_NS4Atom11data_anglesEiPcPiiiiS2_ -Wl,--export-dynamic-symbol=_ZN9LAMMPS_NS4Atom10data_atomsEiPciiiiPdiPii -Wl,--export-dynamic-symbol=_ZN9LAMMPS_NS10DumpCustom11write_linesEiPd -Wl,--export-dynamic-symbol=_ZN9LAMMPS_NS10DumpCustom14convert_stringEiPd ${BASE_LDFLAGS} ${OPENMPI_LDFLAGS} ${FMT_LDFLAGS} ${LAMMPS_LDFLAGS} ${CUDA_LDFLAGS} -ldl

BASE_INCLUDE  := -Iinclude/ -I${GCC_INC}
OPENMPI_INCLUDE := -I${OpenMPI_INC}
LAMMPS_INCLUDE := -I${LAMMPS_SRC}
INCLUDE := ${BASE_INCLUDE} ${OPENMPI_INCLUDE} ${LAMMPS_INCLUDE}

BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps

TARGET   := btree_chromo

SRC      := $(wildcard src/*.cpp)
CU_SRC   := $(wildcard src/*.cu)

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
CU_OBJECTS := $(CU_SRC:%.cu=$(OBJ_DIR)/%.o)
CU_DLINK := $(OBJ_DIR)/cuda_dlink.o
DEPENDENCIES \
         := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCC_FLAGS) -Iinclude/ -c $< -o $@

$(CU_DLINK): $(CU_OBJECTS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCC_FLAGS) -dlink $^ -o $@ -lcudart

$(APP_DIR)/$(TARGET): $(OBJECTS) $(CU_OBJECTS) $(CU_DLINK)
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
