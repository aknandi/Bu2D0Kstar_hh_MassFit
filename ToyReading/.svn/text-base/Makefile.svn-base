HdrSuf = h
SrcSuf = C
ObjSuf = o
ExecSuf = exe 

SRCS = $(wildcard *.$(SrcSuf))
HDRS = $(wildcard *.$(HdrSuf))
OBJSPRESENT = $(wildcard *.$(ObjSuf))
OBJSTOMAKE = $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(SRCS))
EXECS =$(wildcard *.$(ExecSuf))

CC=c++
CFLAGS=-c -Wall
INCLUDES=-I${ROOTSYS}/include
#INCLUDES+=-I/data/lhcb/sw/lcg/external/Boost/1.44.0_python2.6/x86_64-slc5-gcc43-opt/include/boost-1_44/
LIBS=$(shell root-config --libs)
LIBS+= -lRooFitCore -lRooFit

all: Main_plots.exe

Main_plots.exe: $(OBJSTOMAKE)
	$(CC) $(INCLUDES) $^ $(LIBS) -o $@

%.o : %.$(SrcSuf)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean :
	rm $(OBJSPRESENT) $(EXECS)
