#include "algorithm.h"

#include <vector>

using namespace std;

void laplacian_hc_smooth(Object& obj, int times, float alpha, float beta) {
    vector<Point3f> o(obj.vertex);
    vector<Point3f> p(o);
    for (int i = 0; i < times; i++) {
        vector<Point3f> b(o.size());
        vector<Point3f> q(p);
        for (unsigned i = 0; i < o.size(); i++) {
            if (obj.adjacentList[i].size()) {
                Point3f s(0, 0, 0);
                for (unsigned j = 0; j < obj.adjacentList[i].size(); j++)
                    s += q[obj.adjacentList[i][j]];
                p[i] = s / obj.adjacentList[i].size();
            }
            b[i] = p[i] - (o[i] * alpha + q[i] * (1 - alpha));
        }
        for (unsigned i = 0; i < o.size(); i++) {
            if (obj.adjacentList[i].size()) {
                Point3f s(0, 0, 0);
                for (unsigned j = 0; j < obj.adjacentList[i].size(); j++)
                    s += b[obj.adjacentList[i][j]];
                p[i] = p[i] - (b[i] * beta + s * (1 - beta) / obj.adjacentList[i].size());
            }
        }
    }
    obj.vertex.assign(p.begin(), p.end());
    obj.update();
}

void center_positioning(Object& o) {
    Point3f u(o.vertex[0]);
    Point3f d(o.vertex[0]);
    for (unsigned i = 0; i < o.vertex.size(); i++) {
        if (o.vertex[0].x[0] > u.x[0]) u.x[0] = o.vertex[0].x[0];
        if (o.vertex[0].x[0] < d.x[0]) d.x[0] = o.vertex[0].x[0];
        if (o.vertex[0].x[1] > u.x[1]) u.x[1] = o.vertex[0].x[1];
        if (o.vertex[0].x[1] < d.x[1]) d.x[1] = o.vertex[0].x[1];
        if (o.vertex[0].x[2] > u.x[2]) u.x[2] = o.vertex[0].x[2];
        if (o.vertex[0].x[2] < d.x[2]) d.x[2] = o.vertex[0].x[2];
    }
    Point3f mid = (u + d) / 2;
    for (unsigned i = 0; i < o.vertex.size(); i++)
        o.vertex[i] -= mid;
}

inline static Point3f cross_point(const Point3f& x, const Point3f& y, float a, float b, float c, float d) {
    float t = (a * x.x[0] + b * x.x[1] + c * x.x[2] + d) / (a * (x.x[0] - y.x[0]) + b * (x.x[1] - y.x[1]) + c * (x.x[2] - y.x[2]));
    Point3f ret;
    ret.x[0] = x.x[0] + (y.x[0] - x.x[0]) * t;
    ret.x[1] = x.x[1] + (y.x[1] - x.x[1]) * t;
    ret.x[2] = x.x[2] + (y.x[2] - x.x[2]) * t;
    return ret;
}

void partition_by_plane(Object& o, float a, float b, float c, float d) {
    vector<bool> vertexFlag(o.vertex.size());
    int flagCount = 0;
    bool remainFlag;
    for (unsigned i = 0; i < o.vertex.size(); i++) {
        vertexFlag[i] = (a * o.vertex[i].x[0] + b * o.vertex[i].x[1] + c * o.vertex[i].x[2] + d) >= 0;
        if (vertexFlag[i])
            flagCount++;
        else
            flagCount--;
    }
    remainFlag = flagCount >= 0;
    vector<bool> faceFlag(o.face.size());
    for (unsigned i = 0; i < o.face.size(); i++) {
        unsigned remainCount = 0;
        if (vertexFlag[o.face[i].vertexInd[0]] == remainFlag) remainCount++;
        if (vertexFlag[o.face[i].vertexInd[1]] == remainFlag) remainCount++;
        if (vertexFlag[o.face[i].vertexInd[2]] == remainFlag) remainCount++;
        if (remainCount == 3) {
            faceFlag[i] = true;
        } else if (remainCount == 2) {
            unsigned j = 0;
            if (vertexFlag[o.face[i].vertexInd[j]] == remainFlag) j = 1;
            if (vertexFlag[o.face[i].vertexInd[j]] == remainFlag) j = 2;
            Point3f midpoint = (o.vertex[o.face[i].vertexInd[0]] +
                                o.vertex[o.face[i].vertexInd[1]] +
                                o.vertex[o.face[i].vertexInd[2]] -
                                o.vertex[o.face[i].vertexInd[j]]) / 2;
            o.vertex[o.face[i].vertexInd[j]] = cross_point(
                midpoint,
                o.vertex[o.face[i].vertexInd[j]],
                a, b, c, d);
            vertexFlag[o.face[i].vertexInd[j]] = remainFlag;
            faceFlag[i] = true;
        } else if (remainCount == 1) {
            unsigned j = 0;
            if (vertexFlag[o.face[i].vertexInd[j]] != remainFlag) j = 1;
            if (vertexFlag[o.face[i].vertexInd[j]] != remainFlag) j = 2;
            for (unsigned k = 0; k < 3; k++)
                if (k != j) {
                    o.vertex[o.face[i].vertexInd[k]] = cross_point(
                        o.vertex[o.face[i].vertexInd[j]],
                        o.vertex[o.face[i].vertexInd[k]],
                        a, b, c, d
                    );
                    vertexFlag[o.face[i].vertexInd[k]] = remainFlag;
                }
            faceFlag[i] = true;
        } else if (remainCount == 0) {
            faceFlag[i] = false;
        }
    }
    vector<Triangle> lastFace;
    lastFace.swap(o.face);
    for (unsigned i = 0; i < lastFace.size(); i++)
        if (faceFlag[i])
            o.face.push_back(lastFace[i]);
    o.update();
}
