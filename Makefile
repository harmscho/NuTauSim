SRC= ./include
CXXFLAGS = -fPIC -Wall -std=c++0x -I$(SRC)

OBJECTS= Earth.o Table.o

all: Earth.o Table.o Simu_elost

Earth.o: $(SRC)/Earth.cc
	$(CXX) -c $(SRC)/Earth.cc -o Earth.o $(CXXFLAGS)
Table.o: $(SRC)/Table.cc
	$(CXX) -c $(SRC)/Table.cc -o Table.o $(CXXFLAGS)
Simu_elost: Simu_elost.cxx $(OBJECTS)
	$(CXX) Simu_elost.cxx -o Simu_elost $(CXXFLAGS) $(OBJECTS)
clean:
	rm *o Simu_elost

