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
    size_t vertexInd[3];
    size_t textureInd[3];
    bool hasVertex(size_t i) {
        return vertexInd[0] == i || vertexInd[1] == i || vertexInd[2] == i;
    }
    bool hasTexture(size_t i) {
        return textureInd[0] == i || textureInd[1] == i || textureInd[2] == i;
    }
};

class Object {
public:
    vector<Point3f> vertex;
    vector<Point2f> texture;
    vector<Point3f> faceNormal;
    vector<Triangle> face;
    vector<vector<size_t> > adjVertex;
    vector<vector<size_t> > facesOfVertex;
    vector<vector<size_t> > adjFace;
    vector<string> text;
protected:
    void removeUnusedPoints();
    void calculateAdjVertex();
    void calculateAdjFace();
    void calculateFacesOfVertex();
    void calculateFaceNormals();
    void connectVertex(size_t u, size_t v);
    void connectFace(size_t u, size_t v);
public:
    void update();
    void load(istream& in);
    void save(ostream& out);
};

#endif
