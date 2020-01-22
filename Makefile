
HEADERS= struct.h
CFLAGS = $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags)
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)

EDEPSIM = /mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu

EDEPGLIBS = -L$(EDEPSIM)/lib/ -ledepsim -ledepsim_io
EDEPINCDIR = $(EDEPSIM)/include/EDepSim

#all: Digitize Reconstruct Analyze
all: MinDigitize Time Digitize Reconstruct Analyze Ana

struct.cxx: include/struct.h include/Linkdef.h
	cd include && rootcint -f ../src/$@ -c $(CFLAGS) -p $(HEADERS) Linkdef.h && cd ..

libStruct.so: struct.cxx
	g++ -shared -fPIC -o lib/$@ -Iinclude `root-config --ldflags` $(CFLAGS) src/$^ && cp src/struct_rdict.pcm lib/struct_rdict.pcm

MinDigitize: libStruct.so
	g++ src/min_digitization.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom \
	$(EDEPGLIBS) -Llib -lStruct

Time: libStruct.so
	g++ src/time_resolution.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom \
	$(EDEPGLIBS) -Llib -lStruct	
	
Digitize: libStruct.so
	g++ src/digitization.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom \
	$(EDEPGLIBS) -Llib -lStruct

Reconstruct: libStruct.so
	g++ src/reconstruction.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom \
	$(EDEPGLIBS) -Llib -lStruct 

Analyze: libStruct.so
	g++ src/analysis.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom -lEG \
	$(EDEPGLIBS) -Llib -lStruct 
	
Ana: libStruct.so
	g++ src/ana.cpp -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom -lEG \
	$(EDEPGLIBS) -Llib -lStruct 

