CPP = g++
CFLAGS = -Wall -Werror -Wextra -O3 -fPIC -std=c++11
INC = -Isrc
LIBS = -lm -lboost_thread -lboost_system

all: film3 contactrecomb

src/pdb.o: src/pdb.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/pdb.cpp -o src/pdb.o

src/ga.o: src/ga.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/ga.cpp -o src/ga.o

src/direct.o: src/direct.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/direct.cpp -o src/direct.o

src/film3.o: src/film3.cpp
	$(CPP) -O3 -fPIC -Wno-write-strings -Wno-unused-result $(INC) -c src/film3.cpp -o src/film3.o

src/contactrecomb.o: src/contactrecomb.cpp
	$(CPP) -O3 -fPIC -Wno-write-strings -Wno-unused-result $(INC) -c src/contactrecomb.cpp -o src/contactrecomb.o

film3: src/pdb.o src/ga.o src/direct.o src/film3.o	
	$(CPP) $(CFLAGS) $(INC) src/pdb.o src/ga.o src/direct.o src/film3.o ${LIBS} -o film3

contactrecomb: src/pdb.o src/ga.o src/direct.o src/contactrecomb.o	
	$(CPP) $(CFLAGS) $(INC) src/pdb.o src/ga.o src/direct.o src/contactrecomb.o ${LIBS} -o contactrecomb

clean:
	rm film3 contactrecomb src/*.o
