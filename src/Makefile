CSRC   := $(shell ls *.cpp | grep -v ff.cpp)
CC     := g++
FC     := gfortran
CFLAGS := -g -Wno-unused-variable -Wno-write-strings
FFLAGS := 
SFILES := surface
FFFILE := ff.cpp
NLLIB  := ../opt/nlopt-2.4.2/.libs/libnlopt.a

all: ../bin $(patsubst %.cpp, %, $(CSRC)) ConnollyGrid drude-ConnollyGrid clean

%: %.cpp 
	$(CC) $(CFLAGS) $< -o ../bin/$@ $(FFFILE) $(NLLIB)

ConnollyGrid: ConnollyGrid.f
	$(FC) $(FFLAGS) $< -o ../bin/cgrid Surface.f 

drude-ConnollyGrid: drude-ConnollyGrid.f
	$(FC) $(FFLAGS) $< -o ../bin/drude-cgrid drude-Surface.f 

../bin: 
	mkdir -p ../bin

clean:


