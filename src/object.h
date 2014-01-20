#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "point.h"

using namespace std;

class Triangle {
public:
    unsigned vertexInd[3];
    unsigned textureInd[3];
};

class Object {
public:
    vector<Point3f> vertex;
    vector<Point2f> texture;
    vector<Triangle> face;
    vector<vector<unsigned> > adjacentList;
    vector<vector<unsigned> > facesOfVertex;
    vector<string> text;
private:
    void connect(unsigned u, unsigned v);
public:
    void update();
    void load(istream& in);
    void save(ostream& out);
};

#endif
