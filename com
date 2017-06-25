#! /bin/csh -f
echo "Cleaning object files and executables"
rm  ./include/Table.o ./include/Earth.o  ./Simu_elost.o ./Simu_elost

echo "Compiling Simu_elost.cxx"
g++ -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual -pthread -Wno-deprecated-declarations -m64 -c -I ./include Simu_elost.cxx
#g++ -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual -pthread -std=c++11 -Wno-deprecated-declarations -m64 -c -I ./include Simu_elost.cxx

echo "Compiling Table.cc"
g++ -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual -pthread -Wno-deprecated-declarations -m64 -c -o ./include/Table.o ./include/Table.cc
#g++ -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual -pthread -std=c++11 -Wno-deprecated-declarations -m64 -c -o ./include/Table.o ./include/Table.cc

echo "Compiling Earth.cc"
g++ -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual -pthread -Wno-deprecated-declarations -m64 -c -o ./include/Earth.o ./include/Earth.cc 
#g++ -g -O2 -Wall -fPIC -Wshadow -Woverloaded-virtual -pthread -std=c++11 -Wno-deprecated-declarations -m64 -c -o ./include/Earth.o ./include/Earth.cc 

echo "Final"
g++ -g -O2  Simu_elost.o include/Table.o include/Earth.o -pthread -lm -ldl -rdynamic -o Simu_elost 
