CC = gcc
CFLAGS = -g -Wall

all: search cluster

search: search.o lsh.o cube.o io.o math.o tools.o hashtable.o bitmap.o vector.o
	$(CC) -o $@ $^ -lm

cluster: cluster.o lsh.o cube.o io.o math.o tools.o hashtable.o bitmap.o \
	vector.o progbar.o
	$(CC) -o $@ $^ -lm

search.o: lsh.h cube.h io.h math.h tools.h hashtable.h vector.h
cluster.o: cluster.h lsh.h cube.h io.h math.h vector.h tools.h
lsh.o: lsh.h io.h math.h tools.h hashtable.h vector.h
cube.o: cube.h io.h math.h tools.h hashtable.h bitmap.h vector.h
tools.o: tools.h
vector.o: vector.h tools.h
hashtable.o: hashtable.h vector.h math.h
bitmap.o: bitmap.h math.h
math.o: math.h tools.h vector.h
io.o: io.h math.h tools.h vector.h
progbar.o: progbar.h

.PHONY: clean
clean:
	rm -rf search cluster search.o cluster.o lsh.o cube.o cluster.o tools.o \
		vector.o hashtable.o bitmap.o math.o io.o progbar.o
