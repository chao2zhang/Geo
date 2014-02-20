CC = g++
CCFLAGS = -std=c++11
OBJECTS = bitmap.o main.o
TARGET = Geo
GLRenderer: $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS) $(CCFLAGS)
main.o : main.cpp bitmap.hpp
	$(CC) -c $<
bitmap.o : bitmap.cpp bitmap.hpp
	$(CC) -c $<
clean : 
	rm $(TARGET) $(OBJECTS)
