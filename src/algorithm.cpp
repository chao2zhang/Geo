#include "algorithm.h"

#include <vector>
#include <queue>
#include <utility>

using namespace std;

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
    const Point3f&& va = a - vertical_point(x, y, a);
    const Point3f&& vb = b - vertical_point(x, y, b);
    return va * vb / va.length() / vb.length();
}

/*
 * Reverse face normal orientation
 */
inline static void reverse_face_orientation(Object& o, int f) {
    swap(o.face[f].vertex_index[0], o.face[f].vertex_index[1]);
    swap(o.face[f].texture_index[0], o.face[f].texture_index[1]);
    o.face_normal[f] *= -1;
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
    if ((t - x) * (t - y) > 0)
        return false;

    // If t is inside the triangle
    // http://blog.sina.com.cn/s/blog_61feffe10100n65y.html
    const Point3f&& a = o.vertex[o.face[i].vertex_index[2]] - t;
    const Point3f&& b = o.vertex[o.face[i].vertex_index[1]] - t;
    const Point3f&& c = o.vertex[o.face[i].vertex_index[0]] - t;
    const Point3f&& u = cross_product(b, c);
    const Point3f&& v = cross_product(c, a);
    if (u * v < 0) return false;
    const Point3f&& w = cross_product(a, b);
    if (u * w < 0) return false;
    return true;
}

/*
 * (u, v) is the shared edge
 */
inline static void get_shared_edge(const Object& o, int f, int g, int& u, int& v, int& rf, int& rg) {
    u = o.face[f].vertex_index[0];
    if (!o.face[g].has_vertex(u))
        u = o.face[f].vertex_index[1];

    v = o.face[f].vertex_index[1];
    if (!o.face[g].has_vertex(v) || u == v)
        v = o.face[f].vertex_index[2];

    rf = o.face[f].vertex_index[0];
    if (rf == u || rf == v)
        rf = o.face[f].vertex_index[1];
    if (rf == u || rf == v)
        rf = o.face[f].vertex_index[2];

    rg = o.face[g].vertex_index[0];
    if (rg == u || rg == v)
        rg = o.face[g].vertex_index[1];
    if (rg == u || rg == v)
        rg = o.face[g].vertex_index[2];
}

inline static bool is_correct_winding(Object& o, int f, int g) {
    int x, y, rf, rg;
    get_shared_edge(o, f, g, x, y, rf, rg);
    int fx = 0, fy = 0, gx = 0, gy = 0;
    while (o.face[f].vertex_index[fx] != x) fx++;
    while (o.face[f].vertex_index[fy] != y) fy++;
    while (o.face[g].vertex_index[gx] != x) gx++;
    while (o.face[g].vertex_index[gy] != y) gy++;

    return (fx + 1 == fy || fx - 2 == fy) == (gx + 1 == gy || gx - 2 == gy);
}
/*
 * f is the face with correct normal orientation,
 * while g is the face with possibly wrong normal orientation.
 */
inline static bool correct_winding(Object& o, int f, int g) {
    bool answer = is_correct_winding(o, f, g);
    if (answer)
        reverse_face_orientation(o, g);
    return answer;
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

/*
 * If the orientation of a vertex normal is inconsistent with all face normals, the normal is not valid.
 */
inline static bool is_vertex_normal_valid(const Object& o, int v) {
    for (int f : o.faces_of_vertex[v])
        if (o.vertex_normal[v] * o.face_normal[f] < eps)
            return false;
    return true;
}

/*
 * If any vertex appears at least twice in the triangle, it is actually a degenerate triangle
 */
inline static bool is_degenerate_triangle(const Triangle& t) {
    return (t.vertex_index[0] == t.vertex_index[1]
         || t.vertex_index[1] == t.vertex_index[2]
         || t.vertex_index[2] == t.vertex_index[0]);
}

/*
 * If two triangle share a vertex, they should not be considered in self-intersection
 */
inline static bool is_triangle_overlap(const Triangle& t, const Triangle& s) {
    return (t.vertex_index[0] == s.vertex_index[0]
         || t.vertex_index[0] == s.vertex_index[1]
         || t.vertex_index[0] == s.vertex_index[2]
         || t.vertex_index[1] == s.vertex_index[0]
         || t.vertex_index[1] == s.vertex_index[1]
         || t.vertex_index[1] == s.vertex_index[2]
         || t.vertex_index[2] == s.vertex_index[0]
         || t.vertex_index[2] == s.vertex_index[1]
         || t.vertex_index[2] == s.vertex_index[2]
            );
}

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
            reverse_face_orientation(o, i);
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
    for (const Triangle& f : o.face) {
        unsigned remainCount = 0;
        remainCount += vertex_flag[f.vertex_index[0]] == remain_flag;
        remainCount += vertex_flag[f.vertex_index[1]] == remain_flag;
        remainCount += vertex_flag[f.vertex_index[2]] == remain_flag;
        if (remainCount == 2) { // Move the point outside to the border
            int j = 0;
            while (vertex_flag[f.vertex_index[j]] == remain_flag) j++;
            Point3f midpoint = (o.vertex[f.vertex_index[0]] +
                                o.vertex[f.vertex_index[1]] +
                                o.vertex[f.vertex_index[2]] -
                                o.vertex[f.vertex_index[j]]) / 2;
            o.vertex[f.vertex_index[j]] = intersect_point(midpoint, o.vertex[f.vertex_index[j]], p);
            vertex_flag[f.vertex_index[j]] = remain_flag;
        } else if (remainCount == 1) { // Move the points outside to the border
            int j = 0;
            while (vertex_flag[f.vertex_index[j]] != remain_flag) j++;
            for (int k = 0; k < 3; k++)
                if (k != j) {
                    o.vertex[f.vertex_index[k]] = intersect_point(
                        o.vertex[f.vertex_index[j]],
                        o.vertex[f.vertex_index[k]],
                        p);
                    vertex_flag[f.vertex_index[k]] = remain_flag;
                }
        }
    }
    o.update();
}
/*
 * A pair of adjacent triangles with dot product of their normals
 * smaller than threshold is regarded as spike.
 */
void count_spike(Object& o, float dot_product) {
    unsigned spike = 0, total = 0;
    for (int i = 0; i < o.face.size(); i++)
        for (int j : o.adj_face[i]) {
            if (j <= i)
                continue;
            int x, y, ri, rj;
            get_shared_edge(o, i, j, x, y, ri, rj);
            float cos = o.face_normal[i] * o.face_normal[j];
            if (is_correct_winding(o, i, j))
                cos *= -1;
            if (cos < dot_product - eps)
                spike++;
            total++;
        }
    cout << "Detect " << spike << '/' << total << " spikes" << endl;
}

void clear_spike(Object& o, float dot_product, vector<bool>& invalid) {
    unsigned spike = 0, total = 0;
    for (int i = 0; i < o.face.size(); i++)
        for (int j : o.adj_face[i]) {
            if (j <= i)
                continue;
            int x, y, ri, rj;
            get_shared_edge(o, i, j, x, y, ri, rj);
            float cos = o.face_normal[i] * o.face_normal[j];
            if (is_correct_winding(o, i, j))
                cos *= -1;
            if (cos < dot_product - eps) {
                spike++;
                invalid[i] = true;
                invalid[j] = true;
            }
            total++;
        }
    cout << "Clear " << spike << '/' << total << " spikes" << endl;
}

void euclian_circuit(vector<pair<int, int> >& path, vector<vector<int> >& adjacent_graph, int i) {
    while (!adjacent_graph[i].empty()) {
        int j = adjacent_graph[i].back();
        adjacent_graph[i].pop_back();
        path.push_back(pair<int, int>(i, j));
        euclian_circuit(path, adjacent_graph, j);
    }
}

void gen_offset_vertex(const Object& o, float offset, Object& m) {
    m.vertex.clear();
    m.vertex.resize(o.vertex.size());
    for (int i = 0; i < o.vertex.size(); i++)
        m.vertex[i] = o.vertex[i] + o.vertex_normal[i] * offset;
}

void gen_offset_invalid_vertex(const Object& o, const Object& m, vector<bool>& invalid) {
    // Mark Vertex i_m as invalid if Line (i_m, i_o) intersects with any face.
    for (int i = 0; i < o.vertex.size(); i++) {
        if (!is_vertex_normal_valid(o, i)) {
            cout << "Invalid: Vertex " << i << endl;
            invalid[i] = true;
            continue;
        }
        for (int j = 0; j < o.face.size(); j++)
            if (is_intersect(o, m.vertex[i], o.vertex[i], j) && !o.face[j].has_vertex(i)) {
                cout << "Intersect: Vertex " << i << " and Face " << j << endl;
                invalid[i] = true;
                break;
            }
    }
}

void gen_offset_face(const Object& o, Object& m, const vector<bool>& invalid) {
    for (int i = 0; i < o.face.size(); i++) {
        const Triangle& face = o.face[i];
        if (invalid[face.vertex_index[0]] ||
            invalid[face.vertex_index[1]] ||
            invalid[face.vertex_index[2]])
            continue;
        Triangle t;
        t.vertex_index[0] = face.vertex_index[1];
        t.vertex_index[1] = face.vertex_index[0];
        t.vertex_index[2] = face.vertex_index[2];
        // If the offset face is flipped, they should be rejected.
        Point3f n = cross_product(o.vertex[t.vertex_index[1]] - o.vertex[t.vertex_index[0]],
                                  o.vertex[t.vertex_index[2]] - o.vertex[t.vertex_index[0]]).normalize();
        if (n * o.face_normal[i] > -eps)
            continue;
        t.texture_index[0] = face.texture_index[1];
        t.texture_index[1] = face.texture_index[0];
        t.texture_index[2] = face.texture_index[2];
        m.face.push_back(t);
    }
}

template <typename T>
void append_vector(vector<T>& v, const vector<T>& w) {
    v.reserve(v.size() + w.size());
    v.insert(v.end(), w.begin(), w.end());
}

void enclose_offset_mesh(Object& o, Object& m, const vector<bool>& invalid) {
    // Tell which points are at the border,
    // and link new vertices along the border
    int u = o.vertex.size();
    append_vector(o.vertex, m.vertex);
    for (Triangle& f : m.face) {
        f.vertex_index[0] += u;
        f.vertex_index[1] += u;
        f.vertex_index[2] += u;
    }
    append_vector(o.face, m.face);
    vector<bool> border(u, false);
    for (int i = 0; i < u; i++)
        if (o.adj_vertex[i].size() > o.faces_of_vertex[i].size())
            border[i] = true;
    for (int i = 0; i < u; i++)
        if (border[i]) {
            for (int j : o.adj_vertex[i])
                if (j > i && border[j] && is_border_edge(o, i, j)) {
                    if (invalid[i] || invalid[j])
                        continue;
                    Triangle t;
                    t.vertex_index[0] = j;
                    t.vertex_index[1] = i;
                    t.vertex_index[2] = i + u;
                    get_texture_with_border_edge(o, i, j, t.texture_index[1], t.texture_index[0]);
                    t.texture_index[2] = t.texture_index[1];
                    o.face.push_back(t);

                    t.vertex_index[0] = j;
                    t.vertex_index[1] = i + u;
                    t.vertex_index[2] = j + u;
                    get_texture_with_border_edge(o, i, j, t.texture_index[1], t.texture_index[2]);
                    t.texture_index[0] = t.texture_index[2];
                    o.face.push_back(t);
                }
        }
}

void detect_self_intersect(const Object& o, vector<bool>& intersect) {
    for (int i = 0; i < o.face.size(); i++)
        for (int j = 0; j < o.face.size(); j++)
            if (!intersect[j] && !is_triangle_overlap(o.face[i], o.face[j])) {
                if (is_intersect(o, o.vertex[o.face[j].vertex_index[0]], o.vertex[o.face[j].vertex_index[1]], i)
                 || is_intersect(o, o.vertex[o.face[j].vertex_index[1]], o.vertex[o.face[j].vertex_index[2]], i)
                 || is_intersect(o, o.vertex[o.face[j].vertex_index[2]], o.vertex[o.face[j].vertex_index[0]], i)) {
                    intersect[i] = true;
                    intersect[j] = true;
                    break;
                 }
            }
}

int flood_fill_face_group(const Object& o, const vector<bool>& invalid, vector<int>& face_group) {
    int cur_group = 0;
    int max_group = 0;
    int max_group_count = 0;
    for (int i = 0; i < o.face.size(); i++) {
        if (face_group[i] || invalid[i])
            continue;
        ++cur_group;
        face_group[i] = cur_group;
        vector<int> q;
        q.push_back(i);
        int cur_pos = 0;
        while (cur_pos < q.size()) {
            int cur_face = q[cur_pos];
            for (int f : o.adj_face[cur_face])
                if (face_group[f] == 0 && !invalid[f]) {
                    q.push_back(f);
                    face_group[f] = cur_group;
                }
            cur_pos++;
        }
        cout << "Group " << cur_group << " Number " << q.size() << endl;
        if (q.size() > max_group_count) {
            max_group_count = q.size();
            max_group = cur_group;
        }
    }
    cout << "Total groups " << cur_group - 1 << endl;
    return max_group;
}

void slice_group(Object& o, const vector<int>& face_group, int group) {
    vector<Triangle> face;
    face.swap(o.face);
    for (int i = 0; i < face.size(); i++)
        if (face_group[i] == group)
            o.face.push_back(face[i]);
}

void fill_trivial_hole(Object& o, float threshold=0.01) {
}

void mesh_offset(Object& o, float offset) {
    Object m;
    vector<bool> vertex_invalid(o.vertex.size(), false);
    gen_offset_vertex(o, offset, m);
    m.texture = o.texture;
    gen_offset_invalid_vertex(o, m, vertex_invalid);
    gen_offset_face(o, m, vertex_invalid);
    m.update(false);

    vector<bool> face_invalid(m.face.size(), false);
    detect_self_intersect(m, face_invalid);

    vector<int> face_group(m.face.size(), 0);
    clear_spike(m, -0.2588, face_invalid);
    int max_group = flood_fill_face_group(m, face_invalid, face_group);
    cout << "Choose Group " << max_group << endl;
    slice_group(m, face_group, max_group);
    o.vertex = m.vertex;
    o.texture = m.texture;
    o.face = m.face;
//
//    fill_trivial_hole(m);
//
//    enclose_offset_mesh(o, m, invalid);
    o.update();
    laplacian_hc_smooth(o,3);
}
