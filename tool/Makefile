# (C) Copyright Kirill Lykov 2013.
#
# Distributed under the FreeBSD Software License (See accompanying file license.txt) 
#
# Make file for oni2avi. Check include and lib pathes before runing.

SHELL = /bin/sh
.SUFFIXES: .cpp .u
.IGNORE:

ROOT =  opt-cell
EXE =   $(ROOT)

UNAME := $(shell uname)
CC = g++

CCFLAGS = -O0 -g -std=c++11
LINK = g++
SRC = opt_dmitry_cell.cpp
OBJ = opt_dmitry_cell.o
INCLUDES = -I../

$(EXE):	$(OBJ)
	$(LINK) -o $(EXE) $(OBJ) $(LIBS_PATH) $(LIBS)

%.o:%.cpp
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<

clean:
	rm -fr $(ROOT) *.o *.d
