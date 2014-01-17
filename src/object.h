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
    vector<Triangle> surface;
    vector<vector<int> > adjacentList;
    vector<string> text;
private:
    void connect(int u, int v) {
        for (int i = 0; i < adjacentList[u].size(); i++)
            if (adjacentList[u][i] == v)
                return;
        adjacentList[u].push_back(v);
    }
public:
    void update() {
        adjacentList.clear();
        adjacentList.resize(vertex.size());
        for (int i = 0; i < surface.size(); i++) {
            connect(surface[i].vertexInd[0], surface[i].vertexInd[1]);
            connect(surface[i].vertexInd[0], surface[i].vertexInd[2]);
            connect(surface[i].vertexInd[1], surface[i].vertexInd[0]);
            connect(surface[i].vertexInd[1], surface[i].vertexInd[2]);
            connect(surface[i].vertexInd[2], surface[i].vertexInd[0]);
            connect(surface[i].vertexInd[2], surface[i].vertexInd[1]);
        }
    }
    void load(istream& in) {
        vertex.clear();
        texture.clear();
        surface.clear();
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
                surface.push_back(Triangle());
                iss >> surface.back().vertexInd[0];
                iss >> ch;
                iss >> surface.back().textureInd[0];
                iss >> surface.back().vertexInd[1];
                iss >> ch;
                iss >> surface.back().textureInd[1];
                iss >> surface.back().vertexInd[2];
                iss >> ch;
                iss >> surface.back().textureInd[2];
                surface.back().vertexInd[0]--;
                surface.back().vertexInd[1]--;
                surface.back().vertexInd[2]--;
                surface.back().textureInd[0]--;
                surface.back().textureInd[1]--;
                surface.back().textureInd[2]--;
            } else {
                text.push_back(line);
            }
        }
        update();
    }
    void save(ostream& out) const {
        for (int i = 0; i < text.size(); i++)
            out << text[i] << endl;
        for (int i = 0; i < vertex.size(); i++)
            out << "v " << vertex[i].x[0] << ' ' << vertex[i].x[1] << ' ' << vertex[i].x[2] << endl;
        for (int i = 0; i < texture.size(); i++)
            out << "vt " << texture[i].x[0] << ' ' << texture[i].x[1] << endl;
        for (int i = 0; i < surface.size(); i++)
            out << "f " << surface[i].vertexInd[0] + 1 << '/'
                        << surface[i].textureInd[0] + 1 << ' '
                        << surface[i].vertexInd[1] + 1 << '/'
                        << surface[i].textureInd[1] + 1 << ' '
                        << surface[i].vertexInd[2] + 1 << '/'
                        << surface[i].textureInd[2] + 1 << endl;
    }
};

#endif
