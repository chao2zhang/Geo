#include "object.h"
#include "algorithm.h"

#include <ctime>
#include <fstream>

using namespace std;

void Object::connect(int u, int v) {
    for (int i = 0; i < adjacentList[u].size(); i++)
        if (adjacentList[u][i] == v)
            return;
    adjacentList[u].push_back(v);
}
void Object::update() {
    //Remove unused vertexes and textures
    vector<bool> useVertex(vertex.size(), false);
    vector<bool> useTexture(texture.size(), false);
    for (int i = 0; i < face.size(); i++) {
        useVertex[face[i].vertexInd[0]] = true;
        useVertex[face[i].vertexInd[1]] = true;
        useVertex[face[i].vertexInd[2]] = true;
        useVertex[face[i].textureInd[0]] = true;
        useVertex[face[i].textureInd[1]] = true;
        useVertex[face[i].textureInd[2]] = true;
    }
    vector<int> vertexIndex(vertex.size());
    vector<int> textureIndex(texture.size());
    for (int i = 0;)

    // Calculate adjacent list
    adjacentList.clear();
    adjacentList.resize(vertex.size());
    for (int i = 0; i < face.size(); i++) {
        connect(face[i].vertexInd[0], face[i].vertexInd[1]);
        connect(face[i].vertexInd[0], face[i].vertexInd[2]);
        connect(face[i].vertexInd[1], face[i].vertexInd[0]);
        connect(face[i].vertexInd[1], face[i].vertexInd[2]);
        connect(face[i].vertexInd[2], face[i].vertexInd[0]);
        connect(face[i].vertexInd[2], face[i].vertexInd[1]);
    }

    // Calculate faces of a specific vertex
    facesOfVertex.clear();
    facesOfVertex.resize(vertex.size());
    for (int i = 0; i < face.size(); i++) {
        facesOfVertex[face[i].vertexInd[0]].push_back(i);
        facesOfVertex[face[i].vertexInd[1]].push_back(i);
        facesOfVertex[face[i].vertexInd[2]].push_back(i);
    }
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
            iss >> face.back().vertexInd[0];
            iss >> ch;
            iss >> face.back().textureInd[0];
            iss >> face.back().vertexInd[1];
            iss >> ch;
            iss >> face.back().textureInd[1];
            iss >> face.back().vertexInd[2];
            iss >> ch;
            iss >> face.back().textureInd[2];
            face.back().vertexInd[0]--;
            face.back().vertexInd[1]--;
            face.back().vertexInd[2]--;
            face.back().textureInd[0]--;
            face.back().textureInd[1]--;
            face.back().textureInd[2]--;
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
        out << "f " << face[i].vertexInd[0] + 1 << '/'
                    << face[i].textureInd[0] + 1 << ' '
                    << face[i].vertexInd[1] + 1 << '/'
                    << face[i].textureInd[1] + 1 << ' '
                    << face[i].vertexInd[2] + 1 << '/'
                    << face[i].textureInd[2] + 1 << endl;
}
