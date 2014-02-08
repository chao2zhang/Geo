#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "base.h"
#include "point.h"

class Face {
public:
    int vertex_index[3] {0, 0, 0};
    int texture_index[3] {0, 0, 0};
    Face() {
        vertex_index[0] = 0;
        vertex_index[1] = 0;
        vertex_index[2] = 0;
        texture_index[0] = 0;
        texture_index[1] = 0;
        texture_index[2] = 0;
    }
    Face(int u, int v, int w) {
        vertex_index[0] = u;
        vertex_index[1] = v;
        vertex_index[2] = w;
        texture_index[0] = 0;
        texture_index[1] = 0;
        texture_index[2] = 0;
    }
    Face(int u, int v, int w, int r, int s, int f) {
        vertex_index[0] = u;
        vertex_index[1] = v;
        vertex_index[2] = w;
        texture_index[0] = r;
        texture_index[1] = s;
        texture_index[2] = f;
    }
    bool has_vertex(int i) const {
        return vertex_index[0] == i || vertex_index[1] == i || vertex_index[2] == i;
    }
    bool has_texture(int i) const {
        return texture_index[0] == i || texture_index[1] == i || texture_index[2] == i;
    }
};

class Mesh {
public:
    std::vector<Point3f> vertex;
    std::vector<Point2f> texture;
    std::vector<Face> face;
    std::vector<std::string> text;
public:
    std::vector<Point3f> face_normal;
    std::vector<Point3f> vertex_normal;
    std::vector<std::vector<int> > adj_vertex;
    std::vector<std::vector<int> > faces_of_vertex;
    std::vector<std::vector<int> > adj_face;
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
    void update();
    void clean();
    void load(std::istream& in);
    void save(std::ostream& out);
};

#endif
