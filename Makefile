OBJ_DIR      := ./build
APP_DIR      := ./bin
TARGET       := REEF3D
APP          := $(APP_DIR)/$(TARGET)
CXX          := mpicxx
GIT_BRANCH   := $(shell git rev-parse --abbrev-ref HEAD)
GIT_VERSION  := "$(shell git describe --dirty --always --tags)"
HYPRE_DIR    := /usr/local/hypre
CUDA_DIR     := /usr/local/cuda
UMPIRE_DIR   := /usr/local/umpire
CAMP_DIR     := ${UMPIRE_DIR}
EIGEN_DIR    := ThirdParty/eigen-3.3.8

NVCC         := $(CUDA_DIR)/bin/nvcc
NVCCFLAGS    := -std=c++17

CXXFLAGS     := -std=c++17 -DVERSION=\"$(GIT_VERSION)\" -DBRANCH=\"$(GIT_BRANCH)\"

INCLUDE      := -I ${HYPRE_DIR}/include -I ${UMPIRE_DIR}/include -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY

# Shared-library link flags (HYPRE built shared with CUDA/Umpire)
LDFLAGS      := \
    -L${HYPRE_DIR}/lib -lHYPRE \
    -L${UMPIRE_DIR}/lib -lumpire \
    -L${CAMP_DIR}/lib -lcamp \
    -L${CUDA_DIR}/lib64 -L${CUDA_DIR}/lib \
    -lcudart -lcublas -lcusparse -lcurand -lcuda

# Prefer runtime discovery without needing LD_LIBRARY_PATH
RPATH_FLAGS  := \
    -Wl,-rpath,${HYPRE_DIR}/lib \
    -Wl,-rpath,${UMPIRE_DIR}/lib \
    -Wl,-rpath,${CAMP_DIR}/lib \
    -Wl,-rpath,${CUDA_DIR}/lib64 \
    -Wl,-rpath,${CUDA_DIR}/lib

SRC          := $(wildcard src/*.cpp) $(wildcard src/*.cu)
OBJECTS      := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
OBJECTS      := $(OBJECTS:%.cu=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

.PHONY: all clean debug dev info release

.DEFAULT_GOAL := release

all: CXXFLAGS += -O3 -w
all: CXXFLAGS += -DBUILD=\"all\"
all: $(APP)

release: CXXFLAGS += -O3 -DNDEBUG -DEIGEN_NO_DEBUG -march=native -flto -w
release: CXXFLAGS += -DBUILD=\"release\"
release: LDFLAGS += -flto=auto
release: $(APP)

dev: CXXFLAGS += -O3 -Wall -pedantic -Wpedantic -Wextra -Wshadow -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 #-Wold-style-cast 
dev: CXXFLAGS += -DBUILD=\"dev\"
dev: $(APP)

debug: CXXFLAGS += -O0 -g -g3 -Wall
debug: CXXFLAGS += -DBUILD=\"debug\"
debug: $(APP)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -MMD -c $< -o $@

$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -MMD -c $< -o $@

$(APP): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(RPATH_FLAGS)

-include $(DEPENDENCIES)

clean:
	-@rm -rvf $(APP_DIR) $(OBJ_DIR)

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
