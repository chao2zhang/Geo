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
            if (obj.adj_vertex[i].size()) {
                Point3f s(0, 0, 0);
                for (int j : obj.adj_vertex[i])
                    s += q[j];
                p[i] = s / obj.adj_vertex[i].size();
            }
            b[i] = p[i] - (o[i] * alpha + q[i] * (1 - alpha));
        }
        for (int i = 0; i < o.size(); i++) {
            if (obj.adj_vertex[i].size()) {
                Point3f s(0, 0, 0);
                for (int j : obj.adj_vertex[i])
                    s += b[j];
                p[i] = p[i] - (b[i] * beta + s * (1 - beta) / obj.adj_vertex[i].size());
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

// Find ret on Plane p and Line xy
inline static Point3f intersect_point(const Point3f& x, const Point3f& y, const Plane& p) {
    float t = (p.a * x.x[0] + p.b * x.x[1] + p.c * x.x[2] + p.d) /
              (p.a * (x.x[0] - y.x[0]) + p.b * (x.x[1] - y.x[1]) + p.c * (x.x[2] - y.x[2]));
    Point3f ret;
    ret.x[0] = x.x[0] + (y.x[0] - x.x[0]) * t;
    ret.x[1] = x.x[1] + (y.x[1] - x.x[1]) * t;
    ret.x[2] = x.x[2] + (y.x[2] - x.x[2]) * t;
    return ret;
}

inline static float point_value_on_line(const Point3f& x, const Plane& p) {
    return p.a * x.x[0] + p.b * x.x[1] + p.c * x.x[2] + p.d;
}

/*
 * Find Point ret on Line xy satisfying that y - x is vertical to p - ret
 */
inline static Point3f vertical_point(const Point3f& x, const Point3f& y, const Point3f& p) {
    float t = (p.x[0] - x.x[0]) * (y.x[0] - x.x[0]) + (p.x[1] - x.x[1]) * (y.x[1] - x.x[1]) + (p.x[2] - x.x[2]) * (y.x[2] - x.x[2]);
    t /= (y.x[0] - x.x[0]) * (y.x[0] - x.x[0]) + (y.x[1] - x.x[1]) * (y.x[1] - x.x[1]) + (y.x[2] - x.x[2]) * (y.x[2] - x.x[2]);
    Point3f ret;
    ret.x[0] = x.x[0] + (y.x[0] - x.x[0]) * t;
    ret.x[1] = x.x[1] + (y.x[1] - x.x[1]) * t;
    ret.x[2] = x.x[2] + (y.x[2] - x.x[2]) * t;
    return ret;
}

/*
 * (x, y) is their shared edge, calculate cosin value of dihedral angle between triangle (x, y, a) and triangle (x, y, b)
 */
inline static float dihedral_angle_cosin(const Point3f& x, const Point3f& y, const Point3f& a, const Point3f& b) {
    Point3f va = a - vertical_point(x, y, a);
    Point3f vb = b - vertical_point(x, y, b);
    return va * vb / va.length() / vb.length();
}

/*
 * Reverse face normal orientation
 */
inline static void reverse_face_normal_orientation(Object& o, int i) {
    swap(o.face[i].vertex_index[0], o.face[i].vertex_index[1]);
    swap(o.face[i].texture_index[0], o.face[i].texture_index[1]);
    o.face_normal[i] *= -1;
}

/*
 * Tell if Line xy intersects with Face i
 */
inline static bool is_intersect(const Object& o, const Point3f& x, const Point3f& y, int i) {
    #define smaller_than_bounding_box(p,i,axis) \
       (p.x[axis] < o.vertex[o.face[i].vertex_index[0]].x[axis] &&\
        p.x[axis] < o.vertex[o.face[i].vertex_index[1]].x[axis] &&\
        p.x[axis] < o.vertex[o.face[i].vertex_index[2]].x[axis])
    #define larger_than_bounding_box(p,i,axis) \
       (p.x[axis] > o.vertex[o.face[i].vertex_index[0]].x[axis] &&\
        p.x[axis] > o.vertex[o.face[i].vertex_index[1]].x[axis] &&\
        p.x[axis] > o.vertex[o.face[i].vertex_index[2]].x[axis])
    if (smaller_than_bounding_box(x,i,0) && smaller_than_bounding_box(y,i,0))
        return false;
    if (smaller_than_bounding_box(x,i,1) && smaller_than_bounding_box(y,i,1))
        return false;
    if (smaller_than_bounding_box(x,i,2) && smaller_than_bounding_box(y,i,2))
        return false;
    if (larger_than_bounding_box(x,i,0) && larger_than_bounding_box(y,i,0))
        return false;
    if (larger_than_bounding_box(x,i,1) && larger_than_bounding_box(y,i,1))
        return false;
    if (larger_than_bounding_box(x,i,2) && larger_than_bounding_box(y,i,2))
        return false;
    #undef smaller_than_bounding_box
    #undef larger_than_bounding_box
    const Point3f& n = o.face_normal[i];
    Plane p(n.x[0], n.x[1], n.x[2], o.vertex[o.face[i].vertex_index[0]] * n * -1);
    const Point3f&& t = intersect_point(x, y, p);
    if ((t - x) * (t - y) > -eps)
        return false;
    const Point3f&& d1 = o.vertex[o.face[i].vertex_index[0]] - t;
    const Point3f&& d2 = o.vertex[o.face[i].vertex_index[1]] - t;
    const Point3f&& d3 = o.vertex[o.face[i].vertex_index[2]] - t;
    const Point3f&& n1 = cross_product(d1, d2);
    const Point3f&& n2 = cross_product(d2, d3);
    const Point3f&& n3 = cross_product(d3, d1);
    return !(n1 * n2 > -eps && n1 * n3 > -eps && n2 * n3 > -eps);
}

/*
 * (x, y) is the shared edge
 */
inline static void get_shared_edge(const Object& o, int i, int j, int& x, int& y, int& ri, int& rj) {
    x = o.face[i].vertex_index[0];
    if (!o.face[j].has_vertex(x))
        x = o.face[i].vertex_index[1];

    y = o.face[i].vertex_index[1];
    if (!o.face[j].has_vertex(y) || x == y)
        y = o.face[i].vertex_index[2];

    ri = o.face[i].vertex_index[0];
    if (ri == x || ri == y)
        ri = o.face[i].vertex_index[1];
    if (ri == x || ri == y)
        ri = o.face[i].vertex_index[2];

    rj = o.face[j].vertex_index[0];
    if (rj == x || rj == y)
        rj = o.face[j].vertex_index[1];
    if (rj == x || rj == y)
        rj = o.face[j].vertex_index[2];
}

/*
 * i is the face with correct normal orientation,
 * while j is the face with possibly wrong normal orientation.
 */
inline static bool correct_winding(Object& o, int i, int j) {
    int x, y, ri, rj;
    get_shared_edge(o, i, j, x, y, ri, rj);
    int ix = 0, iy = 0, jx = 0, jy = 0;
    while (o.face[i].vertex_index[ix] != x) ix++;
    while (o.face[i].vertex_index[iy] != y) iy++;
    while (o.face[j].vertex_index[jx] != x) jx++;
    while (o.face[j].vertex_index[jy] != y) jy++;

    if ((ix + 1 == iy || ix - 2 == iy) == (jx + 1 == jy || jx - 2 == jy)) {
        reverse_face_normal_orientation(o, j);
        return true;
    }
    return false;
}

void dfs_face_normals(vector<bool>& visited, unsigned& reverse_count, Object& o, int i) {
    if (visited[i])
        return;
    visited[i] = true;
    for (int j : o.adj_face[i]) {
        reverse_count += correct_winding(o, i, j);
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
            reverse_face_normal_orientation(o, i);
    cout << reverse_count << " faces reverted" << endl;
}

void partition_by_plane(Object& o, const Plane& p) {
    vector<bool> vertex_flag(o.vertex.size());
    int flag_count = 0;
    bool remain_flag; // vertex i satisfying that vertex_flag[i] is equal to remain_flag would remain
    for (unsigned i = 0; i < o.vertex.size(); i++) {
        vertex_flag[i] = point_value_on_line(o.vertex[i], p) >= 0;
        if (vertex_flag[i])
            flag_count++;
        else
            flag_count--;
    }
    remain_flag = flag_count >= 0;
    vector<bool> face_flag(o.face.size()); // tells whether a face should remain

    // First remove faces with all vertices outside
    vector<Triangle> last_face;
    last_face.swap(o.face);
    for (int i = 0; i < last_face.size(); i++) {
        if (vertex_flag[last_face[i].vertex_index[0]] == remain_flag ||
            vertex_flag[last_face[i].vertex_index[1]] == remain_flag ||
            vertex_flag[last_face[i].vertex_index[2]] == remain_flag)
        o.face.push_back(last_face[i]);
    }

    // Second round
    for (int i = 0; i < o.face.size(); i++) {
        unsigned remainCount = 0;
        if (vertex_flag[o.face[i].vertex_index[0]] == remain_flag) remainCount++;
        if (vertex_flag[o.face[i].vertex_index[1]] == remain_flag) remainCount++;
        if (vertex_flag[o.face[i].vertex_index[2]] == remain_flag) remainCount++;
        if (remainCount == 2) { // Move the point outside to the border
            int j = 0;
            if (vertex_flag[o.face[i].vertex_index[j]] == remain_flag) j = 1;
            if (vertex_flag[o.face[i].vertex_index[j]] == remain_flag) j = 2;
            Point3f midpoint = (o.vertex[o.face[i].vertex_index[0]] +
                                o.vertex[o.face[i].vertex_index[1]] +
                                o.vertex[o.face[i].vertex_index[2]] -
                                o.vertex[o.face[i].vertex_index[j]]) / 2;
            o.vertex[o.face[i].vertex_index[j]] = intersect_point(midpoint, o.vertex[o.face[i].vertex_index[j]], p);
            vertex_flag[o.face[i].vertex_index[j]] = remain_flag;
        } else if (remainCount == 1) { // Move the points outside to the border
            int j = 0;
            if (vertex_flag[o.face[i].vertex_index[j]] != remain_flag) j = 1;
            if (vertex_flag[o.face[i].vertex_index[j]] != remain_flag) j = 2;
            for (int k = 0; k < 3; k++)
                if (k != j) {
                    o.vertex[o.face[i].vertex_index[k]] = intersect_point(o.vertex[o.face[i].vertex_index[j]],o.vertex[o.face[i].vertex_index[k]],p);
                    vertex_flag[o.face[i].vertex_index[k]] = remain_flag;
                }
        }
    }
    o.update();
}

void count_spikes(Object& o) {
    // A pair of adjacent triangles with dihedral angle smaller than 90 deg is regarded as a spike
    unsigned spike = 0, total_pair = 0;
    for (int i = 0; i < o.face.size(); i++)
        for (int j : o.adj_face[i]) {
            int x, y, ri, rj;
            get_shared_edge(o, i, j, x, y, ri, rj);
            float cos = dihedral_angle_cosin(o.vertex[x], o.vertex[y], o.vertex[ri], o.vertex[rj]);
            if (cos > eps)
                spike++;
            total_pair++;
        }
    cout << "Spikes: " << spike << '/' << total_pair << endl;
}

void euclian_circuit(vector<pair<int, int> >& path, vector<vector<int> >& adjacent_graph, int i) {
    while (!adjacent_graph[i].empty()) {
        int j = adjacent_graph[i].back();
        adjacent_graph[i].pop_back();
        path.push_back(pair<int, int>(i, j));
        euclian_circuit(path, adjacent_graph, j);
    }
}

inline static bool is_border_edge(const Object& o, int u, int v) {
    unsigned count = 0;
    for (int i = 0; i < o.faces_of_vertex[u].size(); i++) {
        count += o.face[o.faces_of_vertex[u][i]].has_vertex(v);
    }
    return count == 1;
}

inline static void get_texture_with_border_edge(const Object& o, int u, int v, int& tu, int& tv) {
    for (int i = 0; i < o.faces_of_vertex[u].size(); i++) {
        if (o.face[o.faces_of_vertex[u][i]].has_vertex(v)) {
            Triangle f = o.face[o.faces_of_vertex[u][i]];
            if (f.vertex_index[0] == u) tu = f.texture_index[0];
            if (f.vertex_index[1] == u) tu = f.texture_index[1];
            if (f.vertex_index[2] == u) tu = f.texture_index[2];
            if (f.vertex_index[0] == v) tv = f.texture_index[0];
            if (f.vertex_index[1] == v) tv = f.texture_index[1];
            if (f.vertex_index[2] == v) tv = f.texture_index[2];
            return;
        }
    }
}

void shell(Object& o, float offset) {
    int u = o.vertex.size();
    int f = o.face.size();
    for (int i = 0; i < u; i++) {
        Point3f n;
        for (int j : o.faces_of_vertex[i]) {
            n += o.face_normal[j];
        }
        n.normalize();
        Point3f m;
        bool flag = true;
        do {
            flag = false;
            m = o.vertex[i] + n * offset;
            for (int j = 0; j < f; j++) {
                if (is_intersect(o, m, o.vertex[i], j) && !o.face[j].has_vertex(i)) {
                    flag = true;
                    cout << "HAHA" << endl;
                    n /= 4;
                    break;
                }
            }
        } while (flag);
        o.vertex.push_back(m);
    }
    for (int i = 0; i < f; i++) {
        Triangle t;
        t.vertex_index[0] = o.face[i].vertex_index[1] + u;
        t.vertex_index[1] = o.face[i].vertex_index[0] + u;
        t.vertex_index[2] = o.face[i].vertex_index[2] + u;
        t.texture_index[0] = o.face[i].texture_index[1];
        t.texture_index[1] = o.face[i].texture_index[0];
        t.texture_index[2] = o.face[i].texture_index[2];
        o.face.push_back(t);
    }

    // Tell which points are at the border
    // border is only for acceleration
    vector<bool> border(u, false);
    for (int i = 0; i < u; i++)
        if (o.adj_vertex[i].size() > o.faces_of_vertex[i].size())
            border[i] = true;

    // Store the graph of border edges
    for (int i = 0; i < u; i++)
        if (border[i]) {
            for (int j : o.adj_vertex[i])
                if (j > i && border[j] && is_border_edge(o, i, j)) {
                    Triangle t;
                    t.vertex_index[0] = j;
                    t.vertex_index[1] = i;
                    t.vertex_index[2] = i + u;
                    get_texture_with_border_edge(o, i, j, t.texture_index[1], t.texture_index[0]);
                    t.texture_index[2] = t.texture_index[1];
                    o.face.push_back(t);
                    Triangle r;
                    r.vertex_index[0] = j;
                    r.vertex_index[1] = i + u;
                    r.vertex_index[2] = j + u;
                    get_texture_with_border_edge(o, i, j, r.texture_index[1], r.texture_index[2]);
                    r.texture_index[0] = r.texture_index[2];
                    o.face.push_back(r);
                }
        }

    // Seems euclian circuit is redundant
    //vector<pair<int, int> > path;
    //for (int i = 0; i < u; i++)
    //    euclian_circuit(path, adj, i);

    o.update();
}
