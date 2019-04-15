edit : main.o emissionprobabilitycomputer.o copynumber.o
	g++ -o edit main.o emissionprobabilitycomputer.o copynumber.o

main.o : src/main.cpp src/emissionprobabilitycomputer.h src/copynumber.h
	g++ -c src/main.cpp

emissionprobabilitycomputer.o: src/emissionprobabilitycomputer.cpp src/copynumber.h
	g++ -c src/emissionprobabilitycomputer.cpp

copynumber.o : src/copynumber.cpp
	g++ -c src/copynumber.cpp

clean :
	rm edit main.o emissionprobabilitycomputer.o copynumber.o
