#include "algorithm.h"

#include <vector>
#include <queue>
#include <utility>

using namespace std;

void laplacian_hc_smooth(Object& obj, int times, float alpha, float beta) {
    vector<Point3f> o(obj.vertex);
    vector<Point3f> p(o);
    for (int i = 0; i < times; i++) {
        vector<Point3f> b(o.size());
        vector<Point3f> q(p);
        for (unsigned i = 0; i < o.size(); i++) {
            if (obj.adjVertex[i].size()) {
                Point3f s(0, 0, 0);
                for (unsigned j = 0; j < obj.adjVertex[i].size(); j++)
                    s += q[obj.adjVertex[i][j]];
                p[i] = s / obj.adjVertex[i].size();
            }
            b[i] = p[i] - (o[i] * alpha + q[i] * (1 - alpha));
        }
        for (unsigned i = 0; i < o.size(); i++) {
            if (obj.adjVertex[i].size()) {
                Point3f s(0, 0, 0);
                for (unsigned j = 0; j < obj.adjVertex[i].size(); j++)
                    s += b[obj.adjVertex[i][j]];
                p[i] = p[i] - (b[i] * beta + s * (1 - beta) / obj.adjVertex[i].size());
            }
        }
    }
    obj.vertex.swap(p);
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

inline static Point3f cross_point(const Point3f& x, const Point3f& y, const Plane& p) {
    float t = (p.a * x.x[0] + p.b * x.x[1] + p.c * x.x[2] + p.d) /
              (p.a * (x.x[0] - y.x[0]) + p.b * (x.x[1] - y.x[1]) + p.c * (x.x[2] - y.x[2]));
    Point3f ret;
    ret.x[0] = x.x[0] + (y.x[0] - x.x[0]) * t;
    ret.x[1] = x.x[1] + (y.x[1] - x.x[1]) * t;
    ret.x[2] = x.x[2] + (y.x[2] - x.x[2]) * t;
    return ret;
}

inline static Point3f vertical_point(const Point3f& x, const Point3f& y, const Point3f& p) {
    // Find ret on (x, y) satisfying that y - x is vertical to p - ret
    float t = (p.x[0] - x.x[0]) * (y.x[0] - x.x[0]) + (p.x[1] - x.x[1]) * (y.x[1] - x.x[1]) + (p.x[2] - x.x[2]) * (y.x[2] - x.x[2]);
    t /= (y.x[0] - x.x[0]) * (y.x[0] - x.x[0]) + (y.x[1] - x.x[1]) * (y.x[1] - x.x[1]) + (y.x[2] - x.x[2]) * (y.x[2] - x.x[2]);
    Point3f ret;
    ret.x[0] = x.x[0] + (y.x[0] - x.x[0]) * t;
    ret.x[1] = x.x[1] + (y.x[1] - x.x[1]) * t;
    ret.x[2] = x.x[2] + (y.x[2] - x.x[2]) * t;
    return ret;
}

inline static float point_value_on_line(const Point3f& x, const Plane& p) {
    return p.a * x.x[0] + p.b * x.x[1] + p.c * x.x[2] + p.d;
}

inline static float dihedral_angle_cosin(const Point3f& x, const Point3f& y, const Point3f& a, const Point3f& b) {
    // (x, y) is their shared edge, calculate cosin of dihedral angle between triangle (x, y, a) and triangle (x, y, b)
    Point3f va = a - vertical_point(x, y, a);
    Point3f vb = b - vertical_point(x, y, b);
    return va * vb / va.length() / vb.length();
}

inline static void reverse_face(Object& o, size_t i) {
    swap(o.face[i].vertexInd[0], o.face[i].vertexInd[1]);
    swap(o.face[i].textureInd[0], o.face[i].textureInd[1]);
    o.faceNormal[i] *= -1;
}

void dfs_face_normals(vector<bool>& visited, unsigned& reverse_count, Object& o, size_t i) {
    if (visited[i])
        return;
    visited[i] = true;
    for (size_t k = 0; k < o.adjFace[i].size(); k++) {
        size_t j = o.adjFace[i][k];
        size_t x = o.face[i].vertexInd[0];
        if (!o.face[j].hasVertex(x)) x = o.face[i].vertexInd[1];
        size_t y = o.face[i].vertexInd[1];
        if (!o.face[j].hasVertex(y) || x == y) y = o.face[i].vertexInd[2];
        // (x, y) are the shared edge of face[i] and face[adjFace[i][j]] now
        size_t a = o.face[i].vertexInd[0];
        if (a == x || a == y) a = o.face[i].vertexInd[1];
        if (a == x || a == y) a = o.face[i].vertexInd[2];
        size_t b = o.face[j].vertexInd[0];
        if (b == x || b == y) b = o.face[j].vertexInd[1];
        if (b == x || b == y) b = o.face[j].vertexInd[2];
        float cos = dihedral_angle_cosin(o.vertex[x], o.vertex[y], o.vertex[a], o.vertex[b]);
        float prod = o.faceNormal[i] * o.faceNormal[j];
        // This mean face[adjFace[i][j]] normal is not correct
        if (prod * cos > 0) {
            reverse_face(o, j);
            reverse_count++;
        }
        dfs_face_normals(visited, reverse_count, o, j);
    }
}

void unify_face_normals(Object& o) {
    // Suppose the orientation of first face in a connected component is correct
    // Fix the orientation of adjacent faces
    // If more than half faces is fixed, reverse all face normals
    // only called after calculateAdjFace()
    if (o.face.empty())
        return;
    vector<bool> visited(o.face.size(), false);
    unsigned reverse_count = 0;
    for (size_t i = 0; i < o.face.size(); i++) {
        dfs_face_normals(visited, reverse_count, o, i);
    }
    if (reverse_count > o.face.size() / 2)
        for (size_t i = 0; i < o.face.size(); i++)
            reverse_face(o, i);
    cout << reverse_count << " faces reverted" << endl;
}

void partition_by_plane(Object& o, const Plane& p) {
    vector<bool> vertexFlag(o.vertex.size());
    int flagCount = 0;
    bool remainFlag;
    for (unsigned i = 0; i < o.vertex.size(); i++) {
        vertexFlag[i] = point_value_on_line(o.vertex[i], p) >= 0;
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
            o.vertex[o.face[i].vertexInd[j]] = cross_point(midpoint, o.vertex[o.face[i].vertexInd[j]], p);
            vertexFlag[o.face[i].vertexInd[j]] = remainFlag;
            faceFlag[i] = true;
        } else if (remainCount == 1) {
            unsigned j = 0;
            if (vertexFlag[o.face[i].vertexInd[j]] != remainFlag) j = 1;
            if (vertexFlag[o.face[i].vertexInd[j]] != remainFlag) j = 2;
            for (unsigned k = 0; k < 3; k++)
                if (k != j) {
                    o.vertex[o.face[i].vertexInd[k]] = cross_point(o.vertex[o.face[i].vertexInd[j]],o.vertex[o.face[i].vertexInd[k]],p);
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
