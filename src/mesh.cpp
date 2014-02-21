#include "mesh.h"

#include <ctime>
#include <iostream>
#include <algorithm>

using std::istream;
using std::ostream;
using std::cout;
using std::vector;
using std::cerr;
using std::endl;
using std::string;
using std::istringstream;
using std::find;

/**
 * Solve ax = b
 */

inline static bool solve_2(float a[][2], float b[], float x[]) {
    float A = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if (fabs(A) < EPS)
        return false;
    x[0] = (b[0] * a[1][1] - a[0][1] * b[1]) / A;
    x[1] = (a[0][0] * b[1] - b[0] * a[1][0]) / A;
    return true;
}

inline static bool solve_3(float a[][3], float b[], float x[]) {
    float A =
            a[0][0] * a[1][1] * a[2][2] +
            a[0][1] * a[1][2] * a[2][0] +
            a[0][2] * a[1][0] * a[2][1] -
            a[0][0] * a[1][2] * a[2][1] -
            a[0][1] * a[1][0] * a[2][2] -
            a[0][2] * a[1][1] * a[2][0];
    if (fabs(A) < EPS)
        return false;
    x[0] = (b[0] * a[1][1] * a[2][2] +
            a[0][1] * a[1][2] * b[2] +
            a[0][2] * b[1] * a[2][1] -
            b[0] * a[1][2] * a[2][1] -
            a[0][1] * b[1] * a[2][2] -
            a[0][2] * a[1][1] * b[2]) / A;
    x[1] = (a[0][0] * b[1] * a[2][2] +
            b[0] * a[1][2] * a[2][0] +
            a[0][2] * a[1][0] * b[2] -
            a[0][0] * a[1][2] * b[2] -
            b[0] * a[1][0] * a[2][2] -
            a[0][2] * b[1] * a[2][0]) / A;
    x[2] = (a[0][0] * a[1][1] * b[2] +
            a[0][1] * b[1] * a[2][0] +
            b[0] * a[1][0] * a[2][1] -
            a[0][0] * b[1] * a[2][1] -
            a[0][1] * a[1][0] * b[2] -
            b[0] * a[1][1] * a[2][0]) / A;
    return true;
}

void Mesh::connect_vertex(int u, int v) {
    if (u == v)
        return;
    if (find(adj_vertex[u].begin(), adj_vertex[u].end(), v) != adj_vertex[u].end())
        return;
    adj_vertex[u].push_back(v);
}

void Mesh::connect_face(int u, int v) {
    if (u == v)
        return;
    if (find(adj_face[u].begin(), adj_face[u].end(), v) != adj_face[u].end())
        return;
    adj_face[u].push_back(v);
}

void Mesh::remove_unused_vertex() {
    vector<bool> use_vertex(vertex.size(), false);
    vector<int> vertex_index(vertex.size());
    for (int i = 0; i < face.size(); ++i) {
        use_vertex[face[i].vertex_index[0]] = true;
        use_vertex[face[i].vertex_index[1]] = true;
        use_vertex[face[i].vertex_index[2]] = true;
    }
    vector<Point3f> last_vertex;
    last_vertex.swap(vertex);
    int c = 0;
    for (int i = 0; i < last_vertex.size(); ++i) {
        if (use_vertex[i]) {
            vertex_index[i] = c++;
            vertex.push_back(last_vertex[i]);
        }
    }
    for (int i = 0; i < face.size(); ++i) {
        face[i].vertex_index[0] = vertex_index[face[i].vertex_index[0]];
        face[i].vertex_index[1] = vertex_index[face[i].vertex_index[1]];
        face[i].vertex_index[2] = vertex_index[face[i].vertex_index[2]];
    }
}

void Mesh::remove_unused_texture() {
    vector<bool> use_texture(texture.size(), false);
    vector<int> texture_index(texture.size());
    for (int i = 0; i < face.size(); ++i) {
        use_texture[face[i].texture_index[0]] = true;
        use_texture[face[i].texture_index[1]] = true;
        use_texture[face[i].texture_index[2]] = true;
    }
    vector<Point2f> last_texture;
    last_texture.swap(texture);
    int c = 0;
    for (int i = 0; i < last_texture.size(); ++i)
        if (use_texture[i]) {
            texture_index[i] = c++;
            texture.push_back(last_texture[i]);
        }
    for (int i = 0; i < face.size(); ++i) {
        face[i].texture_index[0] = texture_index[face[i].texture_index[0]];
        face[i].texture_index[1] = texture_index[face[i].texture_index[1]];
        face[i].texture_index[2] = texture_index[face[i].texture_index[2]];
    }
}

void Mesh::calculate_adj_face() {
    // Only called after calculate_adj_vertex() and calculate_faces_of_vertex()
    #define conn_adj_face(a, b)\
        for (int j : faces_of_vertex[a])\
            for (int k : faces_of_vertex[b])\
                if (j == k) {\
                    connect_face(i, j);\
                    connect_face(j, i);\
                    break;\
                }
    adj_face.clear();
    adj_face.resize(face.size());
    for (int i = 0; i < face.size(); ++i) {
        conn_adj_face(face[i].vertex_index[0], face[i].vertex_index[1]);
        conn_adj_face(face[i].vertex_index[1], face[i].vertex_index[2]);
        conn_adj_face(face[i].vertex_index[2], face[i].vertex_index[0]);
    }
    #undef conn_adj_face
}

void Mesh::calculate_adj_vertex() {
    // Calculate adjacent list
    adj_vertex.clear();
    adj_vertex.resize(vertex.size());
    for (int i = 0; i < face.size(); ++i) {
        connect_vertex(face[i].vertex_index[0], face[i].vertex_index[1]);
        connect_vertex(face[i].vertex_index[0], face[i].vertex_index[2]);
        connect_vertex(face[i].vertex_index[1], face[i].vertex_index[0]);
        connect_vertex(face[i].vertex_index[1], face[i].vertex_index[2]);
        connect_vertex(face[i].vertex_index[2], face[i].vertex_index[0]);
        connect_vertex(face[i].vertex_index[2], face[i].vertex_index[1]);
    }
}

void Mesh::calculate_faces_of_vertex() {
    // Calculate faces of a specific vertex
    faces_of_vertex.clear();
    faces_of_vertex.resize(vertex.size());
    for (int i = 0; i < face.size(); ++i) {
        faces_of_vertex[face[i].vertex_index[0]].push_back(i);
        faces_of_vertex[face[i].vertex_index[1]].push_back(i);
        faces_of_vertex[face[i].vertex_index[2]].push_back(i);
    }
}
void Mesh::calculate_face_normal() {
    face_normal.clear();
    face_normal.resize(face.size());
    for (int i = 0; i < face.size(); ++i) {
        face_normal[i] = cross_product(vertex[face[i].vertex_index[1]] - vertex[face[i].vertex_index[0]],
                              vertex[face[i].vertex_index[2]] - vertex[face[i].vertex_index[0]]).normalize();
    }
}

void Mesh::calculate_vertex_normal() {
    vertex_normal.clear();
    vertex_normal.resize(vertex.size());
    for (int i = 0; i < vertex.size(); ++i) {
        for (int j : faces_of_vertex[i])
            vertex_normal[i] += face_normal[j];
        vertex_normal[i] /= faces_of_vertex[i].size();
    }
}

void Mesh::clean() {
    DEBUG()
    remove_unused_vertex();
    remove_unused_texture();
}

void Mesh::update() {
    DEBUG()
    calculate_adj_vertex();
    calculate_faces_of_vertex();
    calculate_adj_face();
    calculate_face_normal();
    calculate_vertex_normal();
}

void Mesh::load(istream& in) {
    vertex.clear();
    texture.clear();
    face.clear();
    text.clear();
    string line;
    char ch;
    while (getline(in, line)) {
        istringstream iss(line);
        string first;
        iss >> first;
        if (first == "v") {
            vertex.push_back(Point3f());
            iss >> vertex.back();
        } else if (first == "vt") {
            texture.push_back(Point2f());
            iss >> texture.back();
        } else if (first == "f") {
            face.push_back(Face());
            iss >> face.back().vertex_index[0];
            --face.back().vertex_index[0];
            iss >> ch;
            iss >> face.back().texture_index[0];
            --face.back().texture_index[0];

            iss >> face.back().vertex_index[1];
            --face.back().vertex_index[1];
            iss >> ch;
            iss >> face.back().texture_index[1];
            --face.back().texture_index[1];

            iss >> face.back().vertex_index[2];
            --face.back().vertex_index[2];
            iss >> ch;
            iss >> face.back().texture_index[2];
            --face.back().texture_index[2];
        } else {
            text.push_back(line);
        }
    }
    update();
}

void Mesh::save(ostream& out) {
    out.setf(std::ios::fixed);
    for (int i = 0; i < text.size(); ++i)
        out << text[i] << endl;
    for (int i = 0; i < vertex.size(); ++i)
        out << "v " << vertex[i].x[0] << ' ' << vertex[i].x[1] << ' ' << vertex[i].x[2] << endl;
    for (int i = 0; i < texture.size(); ++i)
        out << "vt " << texture[i].x[0] << ' ' << texture[i].x[1] << endl;
    for (int i = 0; i < face.size(); ++i)
        out << "f " << face[i].vertex_index[0] + 1 << '/'
                    << face[i].texture_index[0] + 1 << ' '
                    << face[i].vertex_index[1] + 1 << '/'
                    << face[i].texture_index[1] + 1 << ' '
                    << face[i].vertex_index[2] + 1 << '/'
                    << face[i].texture_index[2] + 1 << endl;
}
