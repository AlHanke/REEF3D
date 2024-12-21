BUILD    := ./build
BIN    	 := ./bin
TARGET   := REEF3D
CXX      := mpicxx

OBJ_DIR   := $(BUILD)
APP_DIR   := $(BIN)
HYPRE_DIR := /usr/local/hypre
EIGEN_DIR := ThirdParty/eigen-3.3.8 
CXXFLAGS := -std=c++11 -O3
# CXXFLAGS += -w
CXXFLAGS += -Wall -Wextra
CXXFLAGS += -Wno-unused-parameter
CXXFLAGS += -Wno-deprecated-declarations ## flowfile_in_filename
CXXFLAGS += -Wno-unused-private-field
CXXFLAGS += -Wno-overloaded-virtual ## ptf, vtks
CXXFLAGS += -Wno-return-type ## VOF_PLIC
CXXFLAGS += -Wno-unused-value ## VOF_PLIC
## Done Errors
CXXFLAGS += -Wno-uninitialized ## iowave_*_active-beach
CXXFLAGS += -Wno-tautological-overlap-compare ## bc_ikepsilon, iowave_pressure, rheology_Herschel-Bulkley
# CXXFLAGS += -Wno-sometimes-uninitialized
# CXXFLAGS += -Wno-logical-op-parentheses
# CXXFLAGS += -Wno-unused-variable
# CXXFLAGS += -Wno-unused-but-set-variable
# CXXFLAGS += -Wno-inaccessible-base
# CXXFLAGS += -Wno-sign-compare
# CXXFLAGS += -Wno-reorder-ctor
# CXXFLAGS += -Wno-int-in-bool-context
# CXXFLAGS += -Wno-macro-redefined
# CXXFLAGS += -Wno-parentheses-equality
# CXXFLAGS += -Wno-self-assign-field
# CXXFLAGS += -Wno-extra-tokens
# CXXFLAGS += -Wno-implicit-fallthrough
# CXXFLAGS += -flto
CXXFLAGS += -ffunction-sections
LDFLAGS  := -L ${HYPRE_DIR}/lib/ -lHYPRE 
# LDFLAGS  += -Wl,-dead_strip -Wl,-print_statistics -Wl,-map,mapfile.txt
LDFLAGS  += -Wl,--gc-sections -Wl,--print-gc-sections
INCLUDE  := -I ${HYPRE_DIR}/include -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY 
SRC      := $(wildcard src/*.cpp)
OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

.PHONY: all build clean debug info

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS = -O0 -w -g -g3
debug: all

clean:
	-@rm -rf $(OBJ_DIR)/*
	-@rm -rf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
