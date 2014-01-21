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
        for (int i = 0; i < o.size(); i++) {
            if (obj.adjVertex[i].size()) {
                Point3f s(0, 0, 0);
                for (int j = 0; j < obj.adjVertex[i].size(); j++)
                    s += q[obj.adjVertex[i][j]];
                p[i] = s / obj.adjVertex[i].size();
            }
            b[i] = p[i] - (o[i] * alpha + q[i] * (1 - alpha));
        }
        for (int i = 0; i < o.size(); i++) {
            if (obj.adjVertex[i].size()) {
                Point3f s(0, 0, 0);
                for (int j = 0; j < obj.adjVertex[i].size(); j++)
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
    for (int i = 0; i < o.vertex.size(); i++) {
        if (o.vertex[0].x[0] > u.x[0]) u.x[0] = o.vertex[0].x[0];
        if (o.vertex[0].x[0] < d.x[0]) d.x[0] = o.vertex[0].x[0];
        if (o.vertex[0].x[1] > u.x[1]) u.x[1] = o.vertex[0].x[1];
        if (o.vertex[0].x[1] < d.x[1]) d.x[1] = o.vertex[0].x[1];
        if (o.vertex[0].x[2] > u.x[2]) u.x[2] = o.vertex[0].x[2];
        if (o.vertex[0].x[2] < d.x[2]) d.x[2] = o.vertex[0].x[2];
    }
    Point3f mid = (u + d) / 2;
    for (int i = 0; i < o.vertex.size(); i++)
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

inline static void reverse_face(Object& o, int i) {
    swap(o.face[i].vertexInd[0], o.face[i].vertexInd[1]);
    swap(o.face[i].textureInd[0], o.face[i].textureInd[1]);
    o.faceNormal[i] *= -1;
}

inline static void get_shared_edge(const Object& o, int i, int j, int& x, int& y, int& a, int& b) {
    // (x, y) is the shared edge
    x = o.face[i].vertexInd[0];
    if (!o.face[j].hasVertex(x))
        x = o.face[i].vertexInd[1];

    y = o.face[i].vertexInd[1];
    if (!o.face[j].hasVertex(y) || x == y)
        y = o.face[i].vertexInd[2];

    a = o.face[i].vertexInd[0];
    if (a == x || a == y)
        a = o.face[i].vertexInd[1];
    if (a == x || a == y)
        a = o.face[i].vertexInd[2];

    b = o.face[j].vertexInd[0];
    if (b == x || b == y)
        b = o.face[j].vertexInd[1];
    if (b == x || b == y)
        b = o.face[j].vertexInd[2];
}

void dfs_face_normals(vector<bool>& visited, unsigned& reverse_count, Object& o, int i) {
    if (visited[i])
        return;
    visited[i] = true;
    for (int k = 0; k < o.adjFace[i].size(); k++) {
        int j = o.adjFace[i][k];
        int x, y, a, b;
        get_shared_edge(o, i, j, x, y, a, b);
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
    for (int i = 0; i < o.face.size(); i++) {
        dfs_face_normals(visited, reverse_count, o, i);
    }
    if (reverse_count > o.face.size() / 2)
        for (int i = 0; i < o.face.size(); i++)
            reverse_face(o, i);
    cout << reverse_count << " faces reverted" << endl;
}

void partition_by_plane(Object& o, const Plane& p) {
    vector<bool> vertexFlag(o.vertex.size());
    int flagCount = 0;
    bool remainFlag; // vertex i satisfying that vertexFlag[i] is equal to remainFlag would remain
    for (unsigned i = 0; i < o.vertex.size(); i++) {
        vertexFlag[i] = point_value_on_line(o.vertex[i], p) >= 0;
        if (vertexFlag[i])
            flagCount++;
        else
            flagCount--;
    }
    remainFlag = flagCount >= 0;
    vector<bool> faceFlag(o.face.size()); // tells whether a face should remain

    // First remove faces with all vertices outside
    vector<Triangle> lastFace;
    lastFace.swap(o.face);
    for (int i = 0; i < lastFace.size(); i++) {
        if (vertexFlag[lastFace[i].vertexInd[0]] == remainFlag ||
            vertexFlag[lastFace[i].vertexInd[1]] == remainFlag ||
            vertexFlag[lastFace[i].vertexInd[2]] == remainFlag)
        o.face.push_back(lastFace[i]);
    }

    //Second round
    for (int i = 0; i < o.face.size(); i++) {
        unsigned remainCount = 0;
        if (vertexFlag[o.face[i].vertexInd[0]] == remainFlag) remainCount++;
        if (vertexFlag[o.face[i].vertexInd[1]] == remainFlag) remainCount++;
        if (vertexFlag[o.face[i].vertexInd[2]] == remainFlag) remainCount++;
        if (remainCount == 2) { // Move the point outside to the border
            int j = 0;
            if (vertexFlag[o.face[i].vertexInd[j]] == remainFlag) j = 1;
            if (vertexFlag[o.face[i].vertexInd[j]] == remainFlag) j = 2;
            Point3f midpoint = (o.vertex[o.face[i].vertexInd[0]] +
                                o.vertex[o.face[i].vertexInd[1]] +
                                o.vertex[o.face[i].vertexInd[2]] -
                                o.vertex[o.face[i].vertexInd[j]]) / 2;
            o.vertex[o.face[i].vertexInd[j]] = cross_point(midpoint, o.vertex[o.face[i].vertexInd[j]], p);
            vertexFlag[o.face[i].vertexInd[j]] = remainFlag;
        } else if (remainCount == 1) { // Move the points outside to the border
            int j = 0;
            if (vertexFlag[o.face[i].vertexInd[j]] != remainFlag) j = 1;
            if (vertexFlag[o.face[i].vertexInd[j]] != remainFlag) j = 2;
            for (int k = 0; k < 3; k++)
                if (k != j) {
                    o.vertex[o.face[i].vertexInd[k]] = cross_point(o.vertex[o.face[i].vertexInd[j]],o.vertex[o.face[i].vertexInd[k]],p);
                    vertexFlag[o.face[i].vertexInd[k]] = remainFlag;
                }
        }
    }
    o.update();
}

void count_spikes(Object& o) {
    // A pair of adjacent triangles with dihedral angle smaller than 90 deg is regarded as a spike
    unsigned spike = 0, total_pair = 0;
    for (int i = 0; i < o.face.size(); i++)
        for (int k = 0; k < o.adjFace[i].size(); k++) {
            int j = o.adjFace[i][k];
            int x, y, a, b;
            get_shared_edge(o, i, j, x, y, a, b);
            float cos = dihedral_angle_cosin(o.vertex[x], o.vertex[y], o.vertex[a], o.vertex[b]);
            if (cos > 0)
                spike++;
            total_pair++;
        }
    cout << "spikes: " << spike << '/' << total_pair << endl;
}

void euclian_circuit(vector<pair<int, int> >& path, vector<vector<int> >& adjacentGraph, int i) {
    while (!adjacentGraph[i].empty()) {
        int j = adjacentGraph[i].back();
        adjacentGraph[i].pop_back();
        path.push_back(pair<int, int>(i, j));
        euclian_circuit(path, adjacentGraph, j);
    }
}

void find_texture_with_border_edge(const Object& o, int u, int v, int& tu, int& tv) {
    for (int i = 0; i < o.facesOfVertex[u].size(); i++) {
        if (o.face[o.facesOfVertex[u][i]].hasVertex(v)) {
            Triangle f = o.face[o.facesOfVertex[u][i]];
            if (f.vertexInd[0] == u) tu = f.textureInd[0];
            if (f.vertexInd[1] == u) tu = f.textureInd[1];
            if (f.vertexInd[2] == u) tu = f.textureInd[2];
            if (f.vertexInd[0] == v) tv = f.textureInd[0];
            if (f.vertexInd[1] == v) tv = f.textureInd[1];
            if (f.vertexInd[2] == v) tv = f.textureInd[2];
            return;
        }
    }
}

void shell(Object& o, float offset) {
    int u = o.vertex.size();
    Point3f n;
    for (int i = 0; i < o.face.size(); i++)
        n += o.faceNormal[i];
    n.normalize();
    for (int i = 0; i < u; i++)
        o.vertex.push_back(o.vertex[i] + n * offset);

    int f = o.face.size();
    for (int i = 0; i < f; i++) {
        Triangle t;
        t.vertexInd[0] = o.face[i].vertexInd[1] + u;
        t.vertexInd[1] = o.face[i].vertexInd[0] + u;
        t.vertexInd[2] = o.face[i].vertexInd[2] + u;
        t.textureInd[0] = o.face[i].textureInd[1];
        t.textureInd[1] = o.face[i].textureInd[0];
        t.textureInd[2] = o.face[i].textureInd[2];
        o.face.push_back(t);
    }

    // Tell which points are at the border
    vector<bool> border(u, false);
    for (int i = 0; i < u; i++)
        if (o.adjVertex[i].size() > o.facesOfVertex[i].size())
            border[i] = true;

    // Store the graph of vertices at the border in adj
    vector<vector<int> > adj(o.adjVertex);
    for (int i = 0; i < u; i++)
        if (border[i]) {
            vector<int> v;
            for (int j = 0; j < adj[i].size(); j++)
                if (border[adj[i][j]])
                    v.push_back(adj[i][j]);
            v.swap(adj[i]);
        } else
            adj[i].clear();
    vector<pair<int, int> > path;
    for (int i = 0; i < u; i++)
        euclian_circuit(path, adj, i);

    // Add faces according to path
    for (int i = 0; i < path.size(); i++) {
        int v = path[i].first;
        int w = path[i].second;
        Triangle t;
        t.vertexInd[0] = w;
        t.vertexInd[1] = v;
        t.vertexInd[2] = v + u;
        find_texture_with_border_edge(o, v, w, t.textureInd[1], t.textureInd[0]);
        t.textureInd[2] = t.textureInd[1];
        o.face.push_back(t);
        t.vertexInd[0] = w;
        t.vertexInd[1] = v + u;
        t.vertexInd[2] = w + u;
        find_texture_with_border_edge(o, v, w, t.textureInd[1], t.textureInd[2]);
        t.textureInd[0] = t.textureInd[2];
        o.face.push_back(t);
    }
    o.update();
}
