# ── Build mode ───────────────────────────────────────────────────────────────
# Objects and the precompiled header live in a per-mode directory. release /
# debug / dev / all use different CXXFLAGS, and make does not track flags, so a
# shared directory silently reuses objects built with the wrong ones. With a
# PCH it is worse than silent: a PCH built with different flags is rejected.
BUILD_MODE   := $(or $(firstword $(filter all release dev debug,$(MAKECMDGOALS))),release)
USE_AMREX    ?= 1
USE_PCH      ?= 1
USE_CCACHE   ?= 0

BUILD_ROOT   := ./build
OBJ_DIR      := $(BUILD_ROOT)/$(BUILD_MODE)$(if $(filter-out 1,$(USE_AMREX)),-noamrex)
APP_DIR      := ./bin
TARGET       := REEF3D
APP          := $(APP_DIR)/$(TARGET)
# mpicxx lives in /opt/homebrew/bin, which is NOT on the PATH macOS gives a
# process launched from a GUI app -- a CodeLite build dies with "mpicxx: command
# not found". Prefer whatever is on PATH (so a non-Homebrew MPI still wins),
# then fall back to the usual install locations. Override with `make CXX=...`.
MPICXX_SEARCH := /opt/homebrew/bin/mpicxx /usr/local/bin/mpicxx /opt/local/bin/mpicxx
CXX          := $(or $(shell command -v mpicxx 2>/dev/null),$(firstword $(wildcard $(MPICXX_SEARCH))),mpicxx)
ifeq ($(USE_CCACHE),1)
# ccache needs to be told that a PCH is in play, or every lookup misses.
# Resolved the same way, and for the same reason.
CCACHE_SEARCH := /opt/homebrew/bin/ccache /usr/local/bin/ccache /opt/local/bin/ccache
CCACHE_BIN   := $(or $(shell command -v ccache 2>/dev/null),$(firstword $(wildcard $(CCACHE_SEARCH))),ccache)
CXX          := CCACHE_SLOPPINESS=pch_defines,time_macros $(CCACHE_BIN) $(CXX)
endif
GIT_BRANCH   := $(shell git rev-parse --abbrev-ref HEAD)
GIT_COMMIT   := $(shell git rev-parse --short=7 HEAD)
GIT_DIRTY    := $(shell git diff --quiet --ignore-submodules HEAD -- || echo -dirty)
GIT_VERSION  := $(GIT_COMMIT)$(GIT_DIRTY)
GIT_DEFS     := -DVERSION=\"$(GIT_VERSION)\" -DBRANCH=\"$(GIT_BRANCH)\"
HYPRE_DIR    := /usr/local/hypre
EIGEN_DIR    := ThirdParty/eigen-5.0.0
# VERSION/BRANCH are deliberately NOT in CXXFLAGS. GIT_DIRTY flips the moment
# the tree is edited, which would change the macro set out from under the PCH
# on every keystroke. They are applied to driver.o / driver_log.o only -- the
# only two TUs that read them -- and those two are built without the PCH.
CXXFLAGS     := -std=c++20
LDFLAGS      := -L ${HYPRE_DIR}/lib/ -lHYPRE
INCLUDE      := -I ${HYPRE_DIR}/include -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY
ifeq ($(USE_AMREX),1)
AMREX_LIBRARY_HOME := ThirdParty/amrex-26.09
CXXFLAGS	 += -DUSE_AMREX=$(USE_AMREX)
LDFLAGS      += -L $(AMREX_LIBRARY_HOME)/lib -lamrex -L /opt/homebrew/lib/gcc/current -lgfortran
INCLUDE      += -I $(AMREX_LIBRARY_HOME)/include
endif
SRC          := $(wildcard src/*.cpp)
OBJECTS      := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

# ── Precompiled header ───────────────────────────────────────────────────────
# NOTE the -c in the rule below: mpicxx expands to `clang++ ... -lmpi`, so
# without it the driver sees two inputs (the header and the library) and fails
# with "cannot specify -o when generating multiple output files".
# ~56% of this project's compile CPU is parsing third-party headers that are
# identical in ~1300 of the 1334 TUs. src/pch.h holds that set. Disable with
# USE_PCH=0 if it ever gets in the way.
PCH_H        := src/pch.h
ifeq ($(USE_PCH),1)
PCH          := $(OBJ_DIR)/pch.pch
PCHFLAGS      = -include-pch $(PCH)
else
PCH          :=
PCHFLAGS      =
endif

# The two TUs that carry -DVERSION/-DBRANCH cannot share the PCH's macro set,
# so they are built without it. ~3 s each.
$(OBJ_DIR)/src/driver.o     : PCHFLAGS :=
$(OBJ_DIR)/src/driver_log.o : PCHFLAGS :=
$(OBJ_DIR)/src/driver.o     : CXXFLAGS += $(GIT_DEFS)
$(OBJ_DIR)/src/driver_log.o : CXXFLAGS += $(GIT_DEFS)

.PHONY: all clean debug dev info release

.DEFAULT_GOAL := release

all: CXXFLAGS += -O3 -g -w
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

ifeq ($(USE_PCH),1)
$(PCH): $(PCH_H)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -x c++-header -c $< -o $@ -MMD -MF $(PCH).d -MT $(PCH)
endif

$(OBJ_DIR)/%.o: %.cpp $(PCH)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(PCHFLAGS) -MMD -c $< -o $@

$(APP): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

-include $(DEPENDENCIES)
ifeq ($(USE_PCH),1)
-include $(PCH).d
endif

clean:
	-@rm -rvf $(APP_DIR) $(BUILD_ROOT)

info:
	@echo "[*] Build mode:      ${BUILD_MODE}     "
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] PCH:             ${PCH}         "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
