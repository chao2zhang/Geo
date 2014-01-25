#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "base.h"
#include "point.h"

using namespace std;

class Triangle {
public:
    int vertex_index[3];
    int texture_index[3];
    bool has_vertex(int i) const {
        return vertex_index[0] == i || vertex_index[1] == i || vertex_index[2] == i;
    }
    bool has_texture(int i) const {
        return texture_index[0] == i || texture_index[1] == i || texture_index[2] == i;
    }
};

class Object {
public:
    vector<Point3f> vertex;
    vector<Point2f> texture;
    vector<Triangle> face;
    vector<string> text;
public:
    vector<Point3f> face_normal;
    vector<Point3f> vertex_normal;
    vector<vector<int> > adj_vertex;
    vector<vector<int> > faces_of_vertex;
    vector<vector<int> > adj_face;
protected:
    void remove_unused_vertex();
    void remove_unused_texture();
    void calculate_adj_vertex();
    void calculate_adj_face();
    void calculate_faces_of_vertex();
    void calculate_face_normal();
    void calculate_vertex_normal();
    void connect_vertex(int u, int v);
    void connect_face(int u, int v);
public:
    void update(bool clean=true);
    void load(istream& in);
    void save(ostream& out);
};

#endif
