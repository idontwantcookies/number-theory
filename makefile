CC = g++
FLAGS = -lgmp -lgmpxx -Wall -pedantic -g3
include:
	$(CC) -o tests.out tests.cpp algoritmos.cpp $(FLAGS)
