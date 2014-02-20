CC = g++
CCFLAGS = -std=c++0x
OBJECTS = mesh.o algorithm.o main.o
TARGET = Geo
Geo: $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(CCFLAGS)
main.o : main.cpp mesh.h algorithm.h
	$(CC) -c $< $(CCFLAGS)
algorithm.o : algorithm.cpp algorithm.h mesh.h base.h point.h
	$(CC) -c $< $(CCFLAGS)
mesh.o : mesh.cpp mesh.h base.h point.h
	$(CC) -c $< $(CCFLAGS)
clean :
	rm $(TARGET) $(OBJECTS)