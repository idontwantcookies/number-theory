CC = g++
FLAGS = -lgmp -lgmpxx -Wall -pedantic -g3
all:
	$(CC) -o tests.out tests.cpp algoritmos.cpp $(FLAGS)
	$(CC) -o generator.out test_generator.cpp algoritmos.cpp $(FLAGS)
