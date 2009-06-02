objects = main.o genomic.o candidate.o SetOfCandidates.o LeftAndRight.o

edit: $(objects)
	g++ -o edit $(objects)
main.o: main.cpp GeneralHeader.h genomic.h candidate.h SetOfCandidates.h
	g++ -c main.cpp
genomic.o: genomic.cpp GeneralHeader.h genomic.h
	g++ -c genomic.cpp
candidate.o: candidate.cpp GeneralHeader.h genomic.h candidate.h
	g++ -c candidate.cpp
SetOfCandidates.o: SetOfCandidates.cpp GeneralHeader.h genomic.h candidate.h SetOfCandidates.h LeftAndRight.h
	g++ -c SetOfCandidates.cpp
LeftAndRight.o: LeftAndRight.cpp LeftAndRight.h
	g++ -c LeftAndRight.cpp
.PHONY: clean
clean:
	-rm edit $(objects)