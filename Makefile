CC=g++
CFLAGS=

DEPS = algorithm.h object.h point.h
OBJ = algorithm.o object.o

Geo: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

algorithm.o: algorithm.cpp algorithm.h object.h point.h
	$(CC) -c -o $@ $< $(CFLAGS)

object.o: object.cpp object.h
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o