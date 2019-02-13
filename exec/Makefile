# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2015-2016
#
#                    Copyright (C) 2016 Stefano Brambilla
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Stefano Brambilla <s.brambilla93@gmail.com> 
#   DATE        : September 2016
# ====================================================================

include ../config.mk

CPPFLAGS=-I../ -I$(mkGetfemInc) -I$(mkBoostInc) 
CXXFLAGS=-std=c++14

ifeq ($(DEBUG),yes)
  OPTFLAGS=-g -Wall
else
  OPTFLAGS=-O3 -march=native
  CPPFLAGS+=-DNDEBUG
endif

LDFLAGS+=-L../ -L$(mkGetfemLib)  -L$(mkBoostLib) -Wl,-rpath=../ 
LDFLAGS+=-L../muparser/lib -Wl,-rpath=../muparser/lib 
#LDFLAGS += -L/opt/lib/samg/
LIBRARIES+=-lgetfem -l3d3d -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lmuparser
#LIBRARIES += -lamg -liomp5

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
EXEC=M3D1D

OUTDIR=vtk

.PHONY: all clean distclean

all: $(EXEC)
	@echo
	@echo Compilation completed!

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -o $@ -c $<

$(EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@  $(LDFLAGS) $(LIBRARIES)

clean:
	$(RM) $(OBJS) $(EXEC) *~ *.log

distclean: clean
	$(RM) *.txt $(OUTDIR)/*
