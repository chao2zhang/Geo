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
    int vertexInd[3];
    int textureInd[3];
};

class Object {
public:
    vector<Point3f> vertex;
    vector<Point2f> texture;
    vector<Triangle> face;
    vector<vector<int> > adjacentList;
    vector<vector<int> > facesOfVertex;
    vector<string> text;
private:
    void connect(int u, int v);
public:
    void update();
    void load(istream& in);
    void save(ostream& out);
};

#endif
