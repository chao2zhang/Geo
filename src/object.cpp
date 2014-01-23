#include "object.h"

#include <ctime>
#include <fstream>

using namespace std;

void Object::connect_vertex(int u, int v) {
    if (u == v)
        return;
    for (int i = 0; i < adj_vertex[u].size(); i++)
        if (adj_vertex[u][i] == v)
            return;
    adj_vertex[u].push_back(v);
}

void Object::connect_face(int u, int v) {
    if (u == v)
        return;
    for (int i = 0; i < adj_face[u].size(); i++)
        if (adj_face[u][i] == v)
            return;
    adj_face[u].push_back(v);
}

void Object::remove_unused_points() {
    __LOG()
    cout << vertex.size() << ' ' << texture.size() << endl;
    //Remove unused vertexes and textures
    vector<bool> use_vertex(vertex.size(), false);
    vector<bool> use_texture(texture.size(), false);
    for (int i = 0; i < face.size(); i++) {
        use_vertex[face[i].vertex_index[0]] = true;
        use_vertex[face[i].vertex_index[1]] = true;
        use_vertex[face[i].vertex_index[2]] = true;
        use_texture[face[i].texture_index[0]] = true;
        use_texture[face[i].texture_index[1]] = true;
        use_texture[face[i].texture_index[2]] = true;
    }

    __LOG()
    vector<int> vertex_index(vertex.size());
    vector<Point3f> last_vertex;
    last_vertex.swap(vertex);
    int c = 0;
    for (int i = 0; i < last_vertex.size(); i++)
        if (use_vertex[i]) {
            vertex_index[i] = c++;
            vertex.push_back(last_vertex[i]);
        }

    __LOG()
    vector<int> texture_index(texture.size());
    vector<Point2f> lastTexture;
    lastTexture.swap(texture);
    c = 0;
    for (int i = 0; i < lastTexture.size(); i++)
        if (use_texture[i]) {
            texture_index[i] = c++;
            texture.push_back(lastTexture[i]);
        }

    __LOG()
    for (int i = 0; i < face.size(); i++) {
        face[i].vertex_index[0] = vertex_index[face[i].vertex_index[0]];
        face[i].vertex_index[1] = vertex_index[face[i].vertex_index[1]];
        face[i].vertex_index[2] = vertex_index[face[i].vertex_index[2]];
        face[i].texture_index[0] = texture_index[face[i].texture_index[0]];
        face[i].texture_index[1] = texture_index[face[i].texture_index[1]];
        face[i].texture_index[2] = texture_index[face[i].texture_index[2]];
    }
}

void Object::calculate_adj_face() {
    __LOG()
    // Only called after calculate_adj_vertex() and calculate_faces_of_vertex()
    adj_face.clear();
    adj_face.resize(face.size());
    for (int i = 0; i < face.size(); i++) {
        int u = face[i].vertex_index[0];
        int v = face[i].vertex_index[1];
        int w = face[i].vertex_index[2];
        for (int j : faces_of_vertex[u])
            for (int k : faces_of_vertex[v])
                if (j == k) {
                    connect_face(i, j);
                    connect_face(j, i);
                    break;
                }
        for (int j : faces_of_vertex[v])
            for (int k : faces_of_vertex[w])
                if (j == k) {
                    connect_face(i, j);
                    connect_face(j, i);
                    break;
                }
        for (int j : faces_of_vertex[w])
            for (int k : faces_of_vertex[u])
                if (j == k) {
                    connect_face(i, j);
                    connect_face(j, i);
                    break;
                }
    }
}

void Object::calculate_adj_vertex() {
    __LOG()
    // Calculate adjacent list
    adj_vertex.clear();
    adj_vertex.resize(vertex.size());
    for (int i = 0; i < face.size(); i++) {
        connect_vertex(face[i].vertex_index[0], face[i].vertex_index[1]);
        connect_vertex(face[i].vertex_index[0], face[i].vertex_index[2]);
        connect_vertex(face[i].vertex_index[1], face[i].vertex_index[0]);
        connect_vertex(face[i].vertex_index[1], face[i].vertex_index[2]);
        connect_vertex(face[i].vertex_index[2], face[i].vertex_index[0]);
        connect_vertex(face[i].vertex_index[2], face[i].vertex_index[1]);
    }
}

void Object::calculate_faces_of_vertex() {
    __LOG()
    // Calculate faces of a specific vertex
    faces_of_vertex.clear();
    faces_of_vertex.resize(vertex.size());
    for (int i = 0; i < face.size(); i++) {
        faces_of_vertex[face[i].vertex_index[0]].push_back(i);
        faces_of_vertex[face[i].vertex_index[1]].push_back(i);
        faces_of_vertex[face[i].vertex_index[2]].push_back(i);
    }
}
void Object::calculate_face_normals() {
    __LOG()
    face_normal.clear();
    face_normal.resize(face.size());
    for (int i = 0; i < face.size(); i++) {
        face_normal[i] = cross_product(vertex[face[i].vertex_index[1]] - vertex[face[i].vertex_index[0]],
                              vertex[face[i].vertex_index[2]] - vertex[face[i].vertex_index[0]]).normalize();
    }
}
void Object::update() {
    remove_unused_points();
    calculate_adj_vertex();
    calculate_faces_of_vertex();
    calculate_adj_face();
    calculate_face_normals();
}

void Object::load(istream& in) {
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
            iss >> vertex.back().x[0];
            iss >> vertex.back().x[1];
            iss >> vertex.back().x[2];
        } else if (first == "vt") {
            texture.push_back(Point2f());
            iss >> texture.back().x[0];
            iss >> texture.back().x[1];
        } else if (first == "f") {
            face.push_back(Triangle());
            iss >> face.back().vertex_index[0];
            iss >> ch;
            iss >> face.back().texture_index[0];
            iss >> face.back().vertex_index[1];
            iss >> ch;
            iss >> face.back().texture_index[1];
            iss >> face.back().vertex_index[2];
            iss >> ch;
            iss >> face.back().texture_index[2];
            face.back().vertex_index[0]--;
            face.back().vertex_index[1]--;
            face.back().vertex_index[2]--;
            face.back().texture_index[0]--;
            face.back().texture_index[1]--;
            face.back().texture_index[2]--;
        } else {
            text.push_back(line);
        }
    }
    update();
}

void Object::save(ostream& out) {
    out.setf(ios::fixed);
    for (int i = 0; i < text.size(); i++)
        out << text[i] << endl;
    for (int i = 0; i < vertex.size(); i++)
        out << "v " << vertex[i].x[0] << ' ' << vertex[i].x[1] << ' ' << vertex[i].x[2] << endl;
    for (int i = 0; i < texture.size(); i++)
        out << "vt " << texture[i].x[0] << ' ' << texture[i].x[1] << endl;
    for (int i = 0; i < face.size(); i++)
        out << "f " << face[i].vertex_index[0] + 1 << '/'
                    << face[i].texture_index[0] + 1 << ' '
                    << face[i].vertex_index[1] + 1 << '/'
                    << face[i].texture_index[1] + 1 << ' '
                    << face[i].vertex_index[2] + 1 << '/'
                    << face[i].texture_index[2] + 1 << endl;
}
