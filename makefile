############################################################
# Makefile for Benchmarking Skyline Algorithms             #
#   Darius Sidlauskas (darius.sidlauskas@epfl.ch)          #
#   Sean Chester (schester@cs.au.dk)                       #
#   Copyright (c) 2014 Aarhus University                   #
############################################################

RM = rm -rf
MV = mv
CP = cp -rf
CC = g++

TARGET = $(OUT)/SkyBench

SRC = $(wildcard src/util/*.cpp) \
	  $(wildcard src/common/*.cpp) \
  	  $(wildcard src/bskytree/*.cpp) \
  	  $(wildcard src/pskyline/*.cpp) \
  	  $(wildcard src/qflow/*.cpp) \
  	  $(wildcard src/hybrid/*.cpp) \
      $(wildcard src/*.cpp)

OBJ = $(addprefix $(OUT)/,$(notdir $(SRC:.cpp=.o)))

OUT = bin

LIB_DIR = # used as -L$(LIB_DIR)
INCLUDES = -I ./src/

LIB = 

# Forces make to look these directories
VPATH = src:src/util:src/bskytree:src/pskyline:src/qflow:src/hybrid:src/common

DIMS=6
V=VERBOSE
DT=0
PROFILER=0

# By default compiling for performance (optimal)
CXXFLAGS = -O3 -m64 -DNDEBUG\
		   -DNUM_DIMS=$(DIMS) -D$(V) -DCOUNT_DT=$(DT) -DPROFILER=$(PROFILER)\
	       -Wno-deprecated -Wno-write-strings -nostdlib -Wpointer-arith \
    	   -Wcast-qual -Wcast-align \
       	   -std=c++0x -fopenmp -mavx
           
LDFLAGS=-m64 -lrt -fopenmp

# Target-specific Variable values:
# Compile for debugging (works with valgrind)
dbg : CXXFLAGS = -O0 -g3 -m64\
	  -DNUM_DIMS=$(DIMS) -DVERBOSE -DCOUNT_DT=0 -DPROFILER=1\
	  -Wno-deprecated -Wno-write-strings -nostdlib -Wpointer-arith \
      -Wcast-qual -Wcast-align -std=c++0x
dbg : all

# All Target
all: $(TARGET)

# Tool invocations
$(TARGET): $(OBJ) $(LIB_DIR)$(LIB)
	@echo 'Building target: $@ (GCC C++ Linker)'
	$(CC) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo 'Finished building target: $@'
	@echo ' '

$(OUT)/%.o: %.cpp
	@echo 'Building file: $< (GCC C++ Compiler)'
	$(CC) $(CXXFLAGS) $(INCLUDES) -c -o"$@" "$<" 
	@echo 'Finished building: $<'
	@echo ' '

clean:
	-$(RM) $(OBJ) $(TARGET) $(addprefix $(OUT)/,$(notdir $(SRC:.cpp=.d)))
	-@echo ' '

deepclean:
	-$(RM) bin/*
	-@echo ' '


.PHONY: all clean deepclean dbg tests
