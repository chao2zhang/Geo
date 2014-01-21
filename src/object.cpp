#include "object.h"
#include "algorithm.h"

#include <ctime>
#include <fstream>

using namespace std;

void Object::connectVertex(size_t u, size_t v) {
    if (u == v)
        return;
    for (size_t i = 0; i < adjVertex[u].size(); i++)
        if (adjVertex[u][i] == v)
            return;
    adjVertex[u].push_back(v);
}

void Object::connectFace(size_t u, size_t v) {
    if (u == v)
        return;
    for (size_t i = 0; i < adjFace[u].size(); i++)
        if (adjFace[u][i] == v)
            return;
    adjFace[u].push_back(v);
}

void Object::removeUnusedPoints() {
    //Remove unused vertexes and textures
    vector<bool> useVertex(vertex.size(), false);
    vector<bool> useTexture(texture.size(), false);
    for (size_t i = 0; i < face.size(); i++) {
        useVertex[face[i].vertexInd[0]] = true;
        useVertex[face[i].vertexInd[1]] = true;
        useVertex[face[i].vertexInd[2]] = true;
        useTexture[face[i].textureInd[0]] = true;
        useTexture[face[i].textureInd[1]] = true;
        useTexture[face[i].textureInd[2]] = true;
    }

    vector<size_t> vertexIndex(vertex.size());
    vector<Point3f> lastVertex;
    lastVertex.swap(vertex);
    size_t c = 0;
    for (size_t i = 0; i < lastVertex.size(); i++)
        if (useVertex[i]) {
            vertexIndex[i] = c++;
            vertex.push_back(lastVertex[i]);
        }
    vector<size_t> textureIndex(texture.size());
    vector<Point2f> lastTexture;
    lastTexture.swap(texture);
    c = 0;
    for (size_t i = 0; i < lastTexture.size(); i++)
        if (useTexture[i]) {
            textureIndex[i] = c++;
            texture.push_back(lastTexture[i]);
        }

    for (size_t i = 0; i < face.size(); i++) {
        face[i].vertexInd[0] = vertexIndex[face[i].vertexInd[0]];
        face[i].vertexInd[1] = vertexIndex[face[i].vertexInd[1]];
        face[i].vertexInd[2] = vertexIndex[face[i].vertexInd[2]];
        face[i].textureInd[0] = textureIndex[face[i].textureInd[0]];
        face[i].textureInd[1] = textureIndex[face[i].textureInd[1]];
        face[i].textureInd[2] = textureIndex[face[i].textureInd[2]];
    }
}

void Object::calculateAdjFace() {
    // Only called after calculateAdjVertex() and calculateFacesOfVertex()
    adjFace.clear();
    adjFace.resize(face.size());
    for (size_t i = 0; i < face.size(); i++) {
        size_t u = face[i].vertexInd[0];
        size_t v = face[i].vertexInd[1];
        size_t w = face[i].vertexInd[2];
        for (size_t j = 0; j < facesOfVertex[u].size(); j++)
            for (size_t k = 0; k < facesOfVertex[v].size(); k++)
                if (facesOfVertex[u][j] == facesOfVertex[v][k]) {
                    connectFace(i, facesOfVertex[u][j]);
                    connectFace(facesOfVertex[u][j], i);
                    break;
                }
        for (size_t j = 0; j < facesOfVertex[v].size(); j++)
            for (size_t k = 0; k < facesOfVertex[w].size(); k++)
                if (facesOfVertex[v][j] == facesOfVertex[w][k]) {
                    connectFace(i, facesOfVertex[v][j]);
                    connectFace(facesOfVertex[v][j], i);
                    break;
                }
        for (size_t j = 0; j < facesOfVertex[w].size(); j++)
            for (size_t k = 0; k < facesOfVertex[u].size(); k++)
                if (facesOfVertex[w][j] == facesOfVertex[u][k]) {
                    connectFace(i, facesOfVertex[w][j]);
                    connectFace(facesOfVertex[w][j], i);
                    break;
                }
    }
}

void Object::calculateAdjVertex() {
    // Calculate adjacent list
    adjVertex.clear();
    adjVertex.resize(vertex.size());
    for (size_t i = 0; i < face.size(); i++) {
        connectVertex(face[i].vertexInd[0], face[i].vertexInd[1]);
        connectVertex(face[i].vertexInd[0], face[i].vertexInd[2]);
        connectVertex(face[i].vertexInd[1], face[i].vertexInd[0]);
        connectVertex(face[i].vertexInd[1], face[i].vertexInd[2]);
        connectVertex(face[i].vertexInd[2], face[i].vertexInd[0]);
        connectVertex(face[i].vertexInd[2], face[i].vertexInd[1]);
    }
}

void Object::calculateFacesOfVertex() {
    // Calculate faces of a specific vertex
    facesOfVertex.clear();
    facesOfVertex.resize(vertex.size());
    for (size_t i = 0; i < face.size(); i++) {
        facesOfVertex[face[i].vertexInd[0]].push_back(i);
        facesOfVertex[face[i].vertexInd[1]].push_back(i);
        facesOfVertex[face[i].vertexInd[2]].push_back(i);
    }
}
void Object::calculateFaceNormals() {
    faceNormal.clear();
    faceNormal.resize(face.size());
    for (size_t i = 0; i < face.size(); i++)
        faceNormal[i] = cross(vertex[face[i].vertexInd[1]] - vertex[face[i].vertexInd[0]],
                              vertex[face[i].vertexInd[2]] - vertex[face[i].vertexInd[0]]).normalize();
}
void Object::update() {
    removeUnusedPoints();
    calculateAdjVertex();
    calculateFacesOfVertex();
    calculateAdjFace();
    calculateFaceNormals();
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
    for (size_t i = 0; i < text.size(); i++)
        out << text[i] << endl;
    for (size_t i = 0; i < vertex.size(); i++)
        out << "v " << vertex[i].x[0] << ' ' << vertex[i].x[1] << ' ' << vertex[i].x[2] << endl;
    for (size_t i = 0; i < texture.size(); i++)
        out << "vt " << texture[i].x[0] << ' ' << texture[i].x[1] << endl;
    for (size_t i = 0; i < face.size(); i++)
        out << "f " << face[i].vertexInd[0] + 1 << '/'
                    << face[i].textureInd[0] + 1 << ' '
                    << face[i].vertexInd[1] + 1 << '/'
                    << face[i].textureInd[1] + 1 << ' '
                    << face[i].vertexInd[2] + 1 << '/'
                    << face[i].textureInd[2] + 1 << endl;
}
