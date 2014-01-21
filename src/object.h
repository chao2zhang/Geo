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
    bool hasVertex(int i) const {
        return vertexInd[0] == i || vertexInd[1] == i || vertexInd[2] == i;
    }
    bool hasTexture(int i) const {
        return textureInd[0] == i || textureInd[1] == i || textureInd[2] == i;
    }
};

class Object {
public:
    vector<Point3f> vertex;
    vector<Point2f> texture;
    vector<Point3f> faceNormal;
    vector<Triangle> face;
    vector<vector<int> > adjVertex;
    vector<vector<int> > facesOfVertex;
    vector<vector<int> > adjFace;
    vector<string> text;
protected:
    void removeUnusedPoints();
    void calculateAdjVertex();
    void calculateAdjFace();
    void calculateFacesOfVertex();
    void calculateFaceNormals();
    void connectVertex(int u, int v);
    void connectFace(int u, int v);
public:
    void update();
    void load(istream& in);
    void save(ostream& out);
};

#endif
