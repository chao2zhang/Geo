#include "algorithm.h"

#include <vector>
#include <list>
#include <queue>
#include <utility>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <algorithm>

using std::vector;
using std::list;
using std::queue;
using std::pair;
using std::swap;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::iterator_traits;
using std::next;

inline static float safe_acos(float x) {
    if (x > ONE_EPS)
        return 0;
    if (x < MINUS_ONE_EPS)
        return PI;
    return acos(x);
}

inline static void get_texture_with_edge(const Mesh& m, int u, int v, int& tu, int& tv) {
    for (int i = 0; i < m.faces_of_vertex[u].size(); ++i) {
        if (m.face[m.faces_of_vertex[u][i]].has_vertex(v)) {
            Face f = m.face[m.faces_of_vertex[u][i]];
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

inline static void get_texture_with_vertex(const Mesh& m, int u, int &tu) {
    for (int i = 0; i < m.faces_of_vertex[u].size(); ++i) {
        const Face& f = m.face[i];
        if (f.vertex_index[0] == u) tu = f.texture_index[0];
        if (f.vertex_index[1] == u) tu = f.texture_index[1];
        if (f.vertex_index[2] == u) tu = f.texture_index[2];
        return;
    }
}


/**
 * Find ret on Plane p and Line xy
 */
inline static Point3f get_intersect_point(const Point3f& x, const Point3f& y, const Plane& p) {
    float f = (p.a * x.x[0] + p.b * x.x[1] + p.c * x.x[2] + p.d) /
              (p.a * (x.x[0] - y.x[0]) + p.b * (x.x[1] - y.x[1]) + p.c * (x.x[2] - y.x[2]));
    Point3f ret;
    ret.x[0] = x.x[0] + (y.x[0] - x.x[0]) * f;
    ret.x[1] = x.x[1] + (y.x[1] - x.x[1]) * f;
    ret.x[2] = x.x[2] + (y.x[2] - x.x[2]) * f;
    return ret;
}

inline static float get_value_on_plane(const Point3f& x, const Plane& p) {
    return p.a * x.x[0] + p.b * x.x[1] + p.c * x.x[2] + p.d;
}

/**
 * Find Point r on Line xy satisfying that Line xy is vertical to Line rp
 */
inline static Point3f get_vertical_point(const Point3f& x, const Point3f& y, const Point3f& p) {
    float f = (p.x[0] - x.x[0]) * (y.x[0] - x.x[0]) + (p.x[1] - x.x[1]) * (y.x[1] - x.x[1]) + (p.x[2] - x.x[2]) * (y.x[2] - x.x[2]);
    f /= (y.x[0] - x.x[0]) * (y.x[0] - x.x[0]) + (y.x[1] - x.x[1]) * (y.x[1] - x.x[1]) + (y.x[2] - x.x[2]) * (y.x[2] - x.x[2]);
    Point3f r;
    r.x[0] = x.x[0] + (y.x[0] - x.x[0]) * f;
    r.x[1] = x.x[1] + (y.x[1] - x.x[1]) * f;
    r.x[2] = x.x[2] + (y.x[2] - x.x[2]) * f;
    return r;
}

/**
 * Find Point r on Plane p that Line xr is vertical to Plane p
 */
inline static Point3f get_vertical_point(const Point3f& x, const Plane& p) {
    float f = (p.a * x.x[0] + p.b * x.x[1] + p.c * x.x[2] + p.d) / (p.a * p.a + p.b * p.b + p.c * p.c);
    Point3f r;
    r.x[0] = x.x[0] - p.a * f;
    r.x[1] = x.x[1] - p.b * f;
    r.x[2] = x.x[2] - p.c * f;
    return r;
}

/**
 * Calculate cosin value of dihedral angle between (x, y, a) and (x, y, b)
 */
inline static float get_dihedral_angle_cosin(const Point3f& x, const Point3f& y, const Point3f& a, const Point3f& b) {
    const Point3f&& va = a - get_vertical_point(x, y, a);
    const Point3f&& vb = b - get_vertical_point(x, y, b);
    return va * vb / va.length() / vb.length();
}

/**
 * (u, v) is the shared edge of Face f and Face g
 * (u, v, rf) constitute Face f
 * (u, v, rg) constitute Face g
 */
inline static void get_shared_edge(const Mesh& m, int f, int g, int& u, int& v, int& rf, int& rg) {
    u = m.face[f].vertex_index[0];
    if (!m.face[g].has_vertex(u))
        u = m.face[f].vertex_index[1];

    v = m.face[f].vertex_index[1];
    if (!m.face[g].has_vertex(v) || u == v)
        v = m.face[f].vertex_index[2];

    rf = m.face[f].vertex_index[0];
    if (rf == u || rf == v)
        rf = m.face[f].vertex_index[1];
    if (rf == u || rf == v)
        rf = m.face[f].vertex_index[2];

    rg = m.face[g].vertex_index[0];
    if (rg == u || rg == v)
        rg = m.face[g].vertex_index[1];
    if (rg == u || rg == v)
        rg = m.face[g].vertex_index[2];
}

/**
 * Tell if Point p lies inside triangle xyz
 */
inline static bool is_inside_triangle(const Point3f& x, const Point3f& y, const Point3f& z, const Point3f& p) {
    const Point3f&& a = x - p;
    const Point3f&& b = y - p;
    const Point3f&& c = z - p;
    const Point3f&& u = cross_product(b, c);
    const Point3f&& v = cross_product(c, a);
    if (u * v < 0)
        return false;
    const Point3f&& w = cross_product(a, b);
    if (u * w < 0)
        return false;
    return true;
}

inline static bool is_less_than_bounding_box_along_axis(const Mesh& m, const Point3f& p, int f, int axis) {
    return
        p.x[axis] < m.vertex[m.face[f].vertex_index[0]].x[axis] &&
        p.x[axis] < m.vertex[m.face[f].vertex_index[1]].x[axis] &&
        p.x[axis] < m.vertex[m.face[f].vertex_index[2]].x[axis];
}

inline static bool is_less_than_bounding_box_along_axis(const Point3f& x, const Point3f& y, const Point3f& p, int axis) {
    return p.x[axis] < x.x[axis] && p.x[axis] < y.x[axis];
}

inline static bool is_more_than_bounding_box_along_axis(const Mesh& m, const Point3f& p, int f, int axis) {
    return
        p.x[axis] > m.vertex[m.face[f].vertex_index[0]].x[axis] &&
        p.x[axis] > m.vertex[m.face[f].vertex_index[1]].x[axis] &&
        p.x[axis] > m.vertex[m.face[f].vertex_index[2]].x[axis];
}

inline static bool is_more_than_bounding_box_along_axis(const Point3f& x, const Point3f& y, const Point3f& p, int axis) {
    return p.x[axis] > x.x[axis] && p.x[axis] > y.x[axis];
}

/**
 * Tell if Line xy intersects with Line pq
 */
inline static bool is_intersected(const Point3f& x, const Point3f& y, const Point3f& p, const Point3f& q) {
    if (is_less_than_bounding_box_along_axis(x, y, p, 0) &&
        is_less_than_bounding_box_along_axis(x, y, q, 0))
        return false;
    if (is_less_than_bounding_box_along_axis(x, y, p, 1) &&
        is_less_than_bounding_box_along_axis(x, y, q, 1))
        return false;
    if (is_less_than_bounding_box_along_axis(x, y, p, 2) &&
        is_less_than_bounding_box_along_axis(x, y, q, 2))
        return false;
    if (is_more_than_bounding_box_along_axis(x, y, p, 0) &&
        is_more_than_bounding_box_along_axis(x, y, q, 0))
        return false;
    if (is_more_than_bounding_box_along_axis(x, y, p, 1) &&
        is_more_than_bounding_box_along_axis(x, y, q, 1))
        return false;
    if (is_more_than_bounding_box_along_axis(x, y, p, 2) &&
        is_more_than_bounding_box_along_axis(x, y, q, 2))
        return false;
    if (cross_product(y - x, p - x).normalize() * cross_product(y - x, q - x).normalize() > 0)
        return false;
    if (cross_product(q - p, x - p).normalize() * cross_product(q - p, y - p).normalize() > 0)
        return false;
    return true;
}
/**
 * Tell if Line xy intersects with Face f
 */
inline static bool is_intersected(const Mesh& m, const Point3f& x, const Point3f& y, int f) {
    // Acceleration:
    // If Line xy is outside the bounding box of Face f,
    // it cannot intersect with the face.
    if (is_less_than_bounding_box_along_axis(m,x,f,0) &&
        is_less_than_bounding_box_along_axis(m,y,f,0))
        return false;
    if (is_less_than_bounding_box_along_axis(m,x,f,1) &&
        is_less_than_bounding_box_along_axis(m,y,f,1))
        return false;
    if (is_less_than_bounding_box_along_axis(m,x,f,2) &&
        is_less_than_bounding_box_along_axis(m,y,f,2))
        return false;
    if (is_more_than_bounding_box_along_axis(m,x,f,0) &&
        is_more_than_bounding_box_along_axis(m,y,f,0))
        return false;
    if (is_more_than_bounding_box_along_axis(m,x,f,1) &&
        is_more_than_bounding_box_along_axis(m,y,f,1))
        return false;
    if (is_more_than_bounding_box_along_axis(m,x,f,2) &&
        is_more_than_bounding_box_along_axis(m,y,f,2))
        return false;

    // If the intersect point of Line xy and Plane p lies between x and y.
    const Point3f& n = m.face_normal[f];
    Plane p(n.x[0], n.x[1], n.x[2], m.vertex[m.face[f].vertex_index[0]] * n * -1);
    const Point3f&& i = get_intersect_point(x, y, p);
    if ((i - x) * (i - y) > 0)
        return false;

    return is_inside_triangle(m.vertex[m.face[f].vertex_index[0]],
                              m.vertex[m.face[f].vertex_index[1]],
                              m.vertex[m.face[f].vertex_index[2]],
                              i);
}

/**
 * If the winding of Face f is correct, tell if the winding of Face g is correct.
 * CONDITION: Face f should be adjacent to Face g.
 */
inline static bool is_correct_winding(Mesh& m, int f, int g) {
    int x, y, rf, rg;
    get_shared_edge(m, f, g, x, y, rf, rg);
    int fx = 0, fy = 0, gx = 0, gy = 0;
    while (m.face[f].vertex_index[fx] != x) fx++;
    while (m.face[f].vertex_index[fy] != y) fy++;
    while (m.face[g].vertex_index[gx] != x) gx++;
    while (m.face[g].vertex_index[gy] != y) gy++;

    return (fx + 1 == fy || fx - 2 == fy) == (gx + 1 == gy || gx - 2 == gy);
}


/**
 * Tell if Vertex u is adjacent to Vertex v
 */
inline static bool is_adjacent(const Mesh& m, int u, int v) {
    return find(m.adj_vertex[u].begin(), m.adj_vertex[u].end(), v) != m.adj_vertex[u].end();
}

/**
 * Tell if Vertex v is a vertex on border of Mesh m
 */
inline static bool is_border_vertex(const Mesh& m, int v) {
    return m.faces_of_vertex[v].size() < m.adj_vertex[v].size();
}

inline static bool is_border_edge(const Mesh& m, int u, int v) {
    unsigned count = 0;
    for (int i = 0; i < m.faces_of_vertex[u].size(); ++i)
        count += m.face[m.faces_of_vertex[u][i]].has_vertex(v);
    return count == 1;
}

/**
 * If the orientation of a vertex normal is inconsistent with all face normals, the normal is not valid.
 */
inline static bool is_vertex_normal_valid(const Mesh& m, int v) {
    for (int f : m.faces_of_vertex[v])
        if (m.vertex_normal[v] * m.face_normal[f] < 0)
            return false;
    return true;
}

/**
 * If any vertex appears twice or more in the face, the face is degenerated.
 */
inline static bool is_degenerate_face(const Face& f) {
    return (f.vertex_index[0] == f.vertex_index[1]
         || f.vertex_index[1] == f.vertex_index[2]
         || f.vertex_index[2] == f.vertex_index[0]);
}

/**
 * If two triangle share a vertex, they should not be considered in self-intersection
 */
inline static bool is_face_overlap(const Face& f, const Face& g) {
    return f.has_vertex(g.vertex_index[0])
        || f.has_vertex(g.vertex_index[1])
        || f.has_vertex(g.vertex_index[2]);
}

/**
 * Reverse face normal orientation
 */
void reverse_face_orientation(Mesh& m, int f) {
    swap(m.face[f].vertex_index[0], m.face[f].vertex_index[1]);
    swap(m.face[f].texture_index[0], m.face[f].texture_index[1]);
    m.face_normal[f] *= -1;
}
/**
 * Face f is the face with correct normal orientation,
 * while Face g is the face with possibly wrong normal orientation.
 */
bool correct_winding(Mesh& m, int f, int g) {
    bool answer = is_correct_winding(m, f, g);
    if (answer)
        reverse_face_orientation(m, g);
    return answer;
}

/**
 * Merge Vertex u with Vertex v(replace all occurrences of v by u)
 */
void merge_vertex(Mesh& m, int u, int v) {
    for (int f : m.faces_of_vertex[v]) {
        if (m.face[f].vertex_index[0] == v)
            m.face[f].vertex_index[0] = u;
        if (m.face[f].vertex_index[1] == v)
            m.face[f].vertex_index[1] = u;
        if (m.face[f].vertex_index[2] == v)
            m.face[f].vertex_index[2] = u;
        // Partial update, remember to call update()
        if (find(m.faces_of_vertex[u].begin(), m.faces_of_vertex[u].end(), f) ==
                 m.faces_of_vertex[u].end())
            m.faces_of_vertex[u].push_back(f);
    }
}

/**
 * Generate vertices at the border.
 */
void gen_border_vertex(const Mesh& m, vector<bool>& is_border) {
    for (int i = 0; i < m.vertex.size(); ++i)
        if (is_border_vertex(m, i))
            is_border[i] = true;
}

/**
 * Generate the graph of only edges at border.
 */
void gen_border_graph(const Mesh& m, vector<vector<int> >& adj_vertex) {
    vector<bool> is_border(m.vertex.size(), false);
    gen_border_vertex(m, is_border);
    for (int i = 0; i < m.vertex.size(); ++i)
        if (is_border[i])
            for (int j : m.adj_vertex[i])
                if (is_border[j] && is_border_edge(m, i, j)) {
                    adj_vertex[i].push_back(j);
                }
}

/**
 * Modified Euler circuit algorithm.
 *
 * Instead of finding a euler circuit, find its sub-circuit where no vertex is visited twice.
 *
 */
bool euler_simple_circuit(vector<int>& path, vector<bool>& visited, vector<vector<int> >& adj_vertex, int s, int u) {
    if (visited[s] && u == s)
        return true;
    if (visited[u])
        return false;
    visited[u] = true;
    path.push_back(u);
    for (int i = 0; i < adj_vertex[u].size(); ++i)
        if (adj_vertex[u][i] != -1) {
            int v = adj_vertex[u][i];
            adj_vertex[u][i] = -1;
            int j = 0;
            for (j = 0; j < adj_vertex[v].size(); ++j)
                if (adj_vertex[v][j] == u) {
                    adj_vertex[v][j] = -1;
                    break;
                }
            if (euler_simple_circuit(path, visited, adj_vertex, s, v))
                return true;
            adj_vertex[v][j] = u;
            adj_vertex[u][i] = v;
        }
    visited[u] = false;
    path.pop_back();
    return false;
}

void euler_simple_circuit(vector<int>& path, vector<vector<int> >& adj_vertex, int v) {
    vector<bool> visited(adj_vertex.size(), false);
    euler_simple_circuit(path, visited, adj_vertex, v, v);
}

void laplacian_smooth(Mesh& m, int times, bool reserve_border, float lambda) {
    vector<Point3f> p(m.vertex);
    for (int f = 0; f < times; ++f) {
        vector<Point3f> q(p);
        for (int i = 0; i < m.vertex.size(); ++i) {
            if (m.adj_vertex[i].size() && m.adj_vertex[i].size() == m.faces_of_vertex[i].size()) {
                Point3f s(0, 0, 0);
                for (int j : m.adj_vertex[i])
                    s += q[j];
                p[i] = lambda * s / m.adj_vertex[i].size() + (1 - lambda) * p[i];
            }
        }
    }
    m.vertex.swap(p);
}

void laplacian_hc_smooth(Mesh& m, int times, bool reserve_border, float alpha, float beta) {
    vector<Point3f> o(m.vertex);
    vector<Point3f> p(o);
    for (int f = 0; f < times; ++f) {
        vector<Point3f> b(o.size());
        vector<Point3f> q(p);
        for (int i = 0; i < o.size(); ++i) {
            if (m.adj_vertex[i].size() && m.adj_vertex[i].size() == m.faces_of_vertex[i].size()) {
                Point3f s(0, 0, 0);
                for (int j : m.adj_vertex[i])
                    s += q[j];
                p[i] = s / m.adj_vertex[i].size();
            }
            b[i] = p[i] - (o[i] * alpha + q[i] * (1 - alpha));
        }
        for (int i = 0; i < o.size(); ++i) {
            if (m.adj_vertex[i].size() && m.adj_vertex[i].size() == m.faces_of_vertex[i].size()) {
                Point3f s(0, 0, 0);
                for (int j : m.adj_vertex[i])
                    s += b[j];
                p[i] = p[i] - (b[i] * beta + s * (1 - beta) / m.adj_vertex[i].size());
            }
        }
    }
    m.vertex.swap(p);
}

void center_positioning_by_bounding_box(Mesh& m) {
    Point3f u(m.vertex[0]);
    Point3f d(m.vertex[0]);
    for (const Point3f& p : m.vertex) {
        if (p.x[0] > u.x[0]) u.x[0] = p.x[0];
        if (p.x[0] < d.x[0]) d.x[0] = p.x[0];
        if (p.x[1] > u.x[1]) u.x[1] = p.x[1];
        if (p.x[1] < d.x[1]) d.x[1] = p.x[1];
        if (p.x[2] > u.x[2]) u.x[2] = p.x[2];
        if (p.x[2] < d.x[2]) d.x[2] = p.x[2];
    }
    Point3f mid = (u + d) / 2;
    for (Point3f& p : m.vertex)
        p -= mid;
}

void center_positioning_by_averaging_vertex(Mesh& m) {
    Point3f s;
    for (const Point3f& p : m.vertex)
        s += p;
    s /= m.vertex.size();
    cout << "Center: " << s.x[0] << ' ' << s.x[1] << ' ' << s.x[2] << endl;
    for (Point3f& p : m.vertex)
        p -= s;
}

void unify_face_normals(Mesh& m) {
    // Suppose the orientation of first face in a connected component is correct
    // Fix the orientation of adjacent faces
    // If more than half faces is fixed, reverse all face normals
    // only called after calculateAdjFace()
    if (m.face.empty())
        return;
    vector<bool> visited(m.face.size(), false);
    unsigned reverse_count = 0;
    for (int i = 0; i < m.face.size(); ++i) {
        if (visited[i])
            continue;
        queue<int> q;
        q.push(i);
        visited[i] = true;
        while (!q.empty()) {
            int i = q.front();
            for (int j : m.adj_face[i]) {
                reverse_count += correct_winding(m, i, j);
                if (!visited[j]) {
                    visited[j] = true;
                    q.push(j);
                }
            }
            q.pop();
        }
    }
    if (reverse_count > m.face.size() / 2)
        for (int i = 0; i < m.face.size(); ++i)
            reverse_face_orientation(m, i);
    cout << reverse_count << " faces reverted" << endl;
}

void partition_by_plane(Mesh& m, const Plane& p) {
    vector<bool> vertex_flag(m.vertex.size());
    int flag_count = 0;
    bool remain_flag; // vertex i satisfying that vertex_flag[i] is equal to remain_flag would remain
    for (unsigned i = 0; i < m.vertex.size(); ++i) {
        vertex_flag[i] = get_value_on_plane(m.vertex[i], p) >= 0;
        if (vertex_flag[i])
            flag_count++;
        else
            flag_count--;
    }
    remain_flag = flag_count >= 0;
    vector<bool> face_flag(m.face.size()); // tells whether a face should remain

    // First remove faces with all vertices outside
    vector<Face> last_face;
    last_face.swap(m.face);
    for (int i = 0; i < last_face.size(); ++i) {
        if (vertex_flag[last_face[i].vertex_index[0]] == remain_flag ||
            vertex_flag[last_face[i].vertex_index[1]] == remain_flag ||
            vertex_flag[last_face[i].vertex_index[2]] == remain_flag)
        m.face.push_back(last_face[i]);
    }

    // Second round
    for (const Face& f : m.face) {
        unsigned remainCount = 0;
        remainCount += vertex_flag[f.vertex_index[0]] == remain_flag;
        remainCount += vertex_flag[f.vertex_index[1]] == remain_flag;
        remainCount += vertex_flag[f.vertex_index[2]] == remain_flag;
        if (remainCount == 2) {
            // Move the point outside to the border
            // Texture might be wrong
            int j = 0;
            while (vertex_flag[f.vertex_index[j]] == remain_flag) ++j;
            Point3f&& midpoint = (m.vertex[f.vertex_index[0]] +
                                m.vertex[f.vertex_index[1]] +
                                m.vertex[f.vertex_index[2]] -
                                m.vertex[f.vertex_index[j]]) / 2;
            m.vertex[f.vertex_index[j]] = get_intersect_point(midpoint, m.vertex[f.vertex_index[j]], p);
            vertex_flag[f.vertex_index[j]] = remain_flag;
        } else if (remainCount == 1) {
            // Move the points outside to the border
            // Texture might be wrong
            int j = 0;
            while (vertex_flag[f.vertex_index[j]] != remain_flag) ++j;
            for (int k = 0; k < 3; ++k)
                if (k != j) {
                    m.vertex[f.vertex_index[k]] = get_intersect_point(
                        m.vertex[f.vertex_index[j]],
                        m.vertex[f.vertex_index[k]],
                        p);
                    vertex_flag[f.vertex_index[k]] = remain_flag;
                }
        }
    }
    m.clean();
    m.update();
}

void project_by_plane(Mesh& m, const Plane& p) {
    int v = m.vertex.size();
    vector<vector<int> > adj_vertex(v);
    gen_border_graph(m, adj_vertex);
    vector<int> max_path;
    DEBUG()
    for (int i = 0; i < v; ++i)
        for (int j = 0; j < adj_vertex[i].size(); ++j)
            if (adj_vertex[i][j] != -1) {
                vector<int> path;
                euler_simple_circuit(path, adj_vertex, i);
                if (path.size() > max_path.size())
                    max_path = path;
            }
    DEBUG()
    for (int i = 0; i < max_path.size(); ++i) {
        Point3f&& r = get_vertical_point(m.vertex[max_path[i]], p);
        m.vertex.push_back(r);
    }
    for (int k = 0; k < max_path.size(); ++k) {
        int i = max_path[k];
        int j = max_path[(k + 1) % max_path.size()];
        int iu = v + k;
        int ju = v + (k + 1) % max_path.size();
        Face f;
        f.vertex_index[0] = j;
        f.vertex_index[1] = i;
        f.vertex_index[2] = iu;
        get_texture_with_edge(m, i, j, f.texture_index[1], f.texture_index[0]);
        f.texture_index[2] = f.texture_index[1];
        m.face.push_back(f);

        f.vertex_index[0] = j;
        f.vertex_index[1] = iu;
        f.vertex_index[2] = ju;
        get_texture_with_edge(m, i, j, f.texture_index[1], f.texture_index[2]);
        f.texture_index[0] = f.texture_index[2];
        m.face.push_back(f);
    }
    m.update();
}

/*
 * A pair of adjacent triangles with dot product of their normals
 * smaller than threshold is regarded as spike.
 */
void count_spike(Mesh& m, float dot_product) {
    unsigned spike = 0, total = 0;
    for (int i = 0; i < m.face.size(); ++i)
        for (int j : m.adj_face[i]) {
            if (j <= i)
                continue;
            int x, y, ri, rj;
            get_shared_edge(m, i, j, x, y, ri, rj);
            float cos = m.face_normal[i] * m.face_normal[j];
            if (is_correct_winding(m, i, j))
                cos *= -1;
            if (cos < dot_product)
                spike++;
            total++;
        }
    cout << "Detect " << spike << '/' << total << " spikes" << endl;
}

void mark_spike(Mesh& m, float dot_product, vector<bool>& is_invalid) {
    unsigned spike = 0, total = 0;
    for (int i = 0; i < m.face.size(); ++i)
        for (int j : m.adj_face[i]) {
            if (j <= i)
                continue;
            int x, y, ri, rj;
            get_shared_edge(m, i, j, x, y, ri, rj);
            float cos = m.face_normal[i] * m.face_normal[j];
            if (is_correct_winding(m, i, j))
                cos *= -1;
            if (cos < dot_product) {
                spike++;
                is_invalid[i] = true;
                is_invalid[j] = true;
            }
            total++;
        }
    cout << "Clear " << spike << '/' << total << " spikes" << endl;
}

void gen_offset_vertex(const Mesh& m, float offset, Mesh& n) {
    n.vertex.clear();
    n.vertex.resize(m.vertex.size());
    for (int i = 0; i < m.vertex.size(); ++i)
        n.vertex[i] = m.vertex[i] + m.vertex_normal[i] * offset;
}

void gen_offset_invalid_vertex(const Mesh& m, const Mesh& n, float offset, vector<bool>& is_invalid) {
    // Mark Vertex i_m as invalid if Line (i_m, i_n) intersects with any face.
    for (int i = 0; i < m.vertex.size(); ++i) {
        if (!is_vertex_normal_valid(m, i)) {
            cout << "Invalid Normal: Vertex " << i << endl;
            is_invalid[i] = true;
            continue;
        }
        if (is_invalid[i])
            continue;
        for (int j = 0; j < m.face.size(); ++j)
            if (is_intersected(m, m.vertex[i], n.vertex[i], j) && !m.face[j].has_vertex(i)) {
                cout << "Intersect: Vertex " << i << " and Face " << j << endl;
                is_invalid[i] = true;
                break;
            }
    }
}

void gen_offset_face(const Mesh& m, Mesh& n, const vector<bool>& is_invalid) {
    for (int i = 0; i < m.face.size(); ++i) {
        const Face& face = m.face[i];
        if (is_invalid[face.vertex_index[0]] ||
            is_invalid[face.vertex_index[1]] ||
            is_invalid[face.vertex_index[2]])
            continue;
        Face f;
        f.vertex_index[0] = face.vertex_index[1];
        f.vertex_index[1] = face.vertex_index[0];
        f.vertex_index[2] = face.vertex_index[2];
        // If the offset face is flipped, they should be rejected.
        Point3f c = cross_product(m.vertex[f.vertex_index[1]] - m.vertex[f.vertex_index[0]],
                                  m.vertex[f.vertex_index[2]] - m.vertex[f.vertex_index[0]]).normalize();
        if (c * m.face_normal[i] > 0)
            continue;
        f.texture_index[0] = face.texture_index[1];
        f.texture_index[1] = face.texture_index[0];
        f.texture_index[2] = face.texture_index[2];
        n.face.push_back(f);
    }
}

void enclose_offset_mesh(Mesh& m, Mesh& n, const vector<bool>& is_invalid) {
    // Tell which points are at the border,
    // and link new vertices along the border
    int v = m.vertex.size();
    vector<vector<int> > adj_vertex(v);
    gen_border_graph(m, adj_vertex);
    for (const Point3f& p : n.vertex)
        m.vertex.push_back(p);
    for (Face& f : n.face) {
        f.vertex_index[0] += v;
        f.vertex_index[1] += v;
        f.vertex_index[2] += v;
    }
    for (const Face& f : n.face)
        m.face.push_back(f);
    for (int i = 0; i < v; ++i)
        for (int j : adj_vertex[i])
            if (j > i && !is_invalid[i] && !is_invalid[j]) {
                Face f;
                f.vertex_index[0] = j;
                f.vertex_index[1] = i;
                f.vertex_index[2] = i + v;
                get_texture_with_edge(m, i, j, f.texture_index[1], f.texture_index[0]);
                f.texture_index[2] = f.texture_index[1];
                m.face.push_back(f);

                f.vertex_index[0] = j;
                f.vertex_index[1] = i + v;
                f.vertex_index[2] = j + v;
                get_texture_with_edge(m, i, j, f.texture_index[1], f.texture_index[2]);
                f.texture_index[0] = f.texture_index[2];
                m.face.push_back(f);
            }
}

void detect_self_intersect(const Mesh& m, vector<bool>& is_int) {
    for (int i = 0; i < m.face.size(); ++i)
        for (int j = 0; j < m.face.size(); ++j)
            if (!is_int[j] && !is_face_overlap(m.face[i], m.face[j])) {
                if (is_intersected(m, m.vertex[m.face[j].vertex_index[0]], m.vertex[m.face[j].vertex_index[1]], i)
                 || is_intersected(m, m.vertex[m.face[j].vertex_index[1]], m.vertex[m.face[j].vertex_index[2]], i)
                 || is_intersected(m, m.vertex[m.face[j].vertex_index[2]], m.vertex[m.face[j].vertex_index[0]], i)) {
                    is_int[i] = true;
                    is_int[j] = true;
                    break;
                 }
            }
}

/*
 * All faces in one cluster shares at least one edge of another face in the same cluster.
 */
int flood_fill_face_group(const Mesh& m, vector<int>& face_group) {
    int cur_group = 0;
    for (int i = 0; i < m.face.size(); ++i) {
        if (face_group[i])
            continue;
        ++cur_group;
        face_group[i] = cur_group;
        vector<int> q;
        q.push_back(i);
        int cur_pos = 0;
        while (cur_pos < q.size()) {
            int cur_face = q[cur_pos];
            for (int f : m.adj_face[cur_face])
                if (face_group[f] == 0) {
                    q.push_back(f);
                    face_group[f] = cur_group;
                }
            ++cur_pos;
        }
    }
    return cur_group;
}

void fill_hole(Mesh& m, const vector<int>& path) {
    if (path.size() <= 2)
        return;
    for (int i = 1 ; i < path.size() - 1; ++i) {
        Face f;
        f.vertex_index[0] = path[0];
        f.vertex_index[1] = path[i];
        f.vertex_index[2] = path[i+1];
        get_texture_with_edge(m, path[0], path[1], f.texture_index[0], f.texture_index[1]);
        get_texture_with_edge(m, path[i], path[i+1], f.texture_index[1], f.texture_index[2]);
        m.face.push_back(f);
    }
}

void fill_trivial_hole(Mesh& m, float threshold=0.02) {
    int u = m.vertex.size();
    vector<vector<int> > adj_vertex(u);
    gen_border_graph(m, adj_vertex);
    for (int i = 0; i < u; ++i)
        for (int j = 0; j < adj_vertex[i].size(); ++j)
            if (adj_vertex[i][j] != -1) {
                vector<int> path;
                euler_simple_circuit(path, adj_vertex, i);
                if (path.size() < threshold * u) {
                    fill_hole(m, path);
                }
            }
}

void mesh_offset(Mesh& m, float offset) {
    Mesh n;
    vector<bool> is_vertex_invalid(m.vertex.size(), false);
    gen_offset_vertex(m, offset, n);
    n.texture = m.texture;
    n.update();
    gen_offset_invalid_vertex(m, n, offset, is_vertex_invalid);
    gen_offset_face(m, n, is_vertex_invalid);
    n.update();
//    vector<bool> face_invalid(n.face.size(), false);
//    detect_self_intersect(n, face_invalid);
//    vector<int> face_group(n.face.size(), 0);
//    n.update();
//    //laplacian_smooth(n, 50);
//    flood_fill_face_group(n, face_invalid, face_group, 10);
    enclose_offset_mesh(m, n, is_vertex_invalid);
    m.clean();
    m.update();
}

void remove_face_by_plane(Mesh& m, const Plane& p) {
    vector<bool> is_vertex_on_plane(m.vertex.size(), false);
    for (int i = 0; i < m.vertex.size(); ++i) {
        if (fabs(get_value_on_plane(m.vertex[i], p)) < EPS)
            is_vertex_on_plane[i] = true;
    }
    vector<Face> face;
    for (int i = 0; i < m.face.size(); ++i) {
        if (!is_vertex_on_plane[m.face[i].vertex_index[0]] ||
            !is_vertex_on_plane[m.face[i].vertex_index[1]] ||
            !is_vertex_on_plane[m.face[i].vertex_index[2]])
            face.push_back(m.face[i]);
    }
    m.face.swap(face);
    m.clean();
    m.update();
}

void remove_face_by_largest_component(Mesh& m) {
    vector<int> face_group(m.face.size(), 0);
    flood_fill_face_group(m, face_group);
    vector<int> component_size(m.face.size() + 1, 0);
    int max_group = 0;
    int max_group_size = 0;
    for (int g : face_group) {
        ++component_size[g];
        if (component_size[g] > max_group_size) {
            max_group = g;
            max_group_size = component_size[g];
        }
    }
    vector<Face> face;
    for (int i = 0; i < face_group.size(); ++i)
        if (face_group[i] == max_group)
            face.push_back(m.face[i]);
    m.face.swap(face);
    m.clean();
    m.update();
}

void remove_face_by_long_edge(Mesh& m, float threshold) {

}

void get_max_border_path_by_plane(Mesh& m, const Plane& p, list<int>& max_path) {
    vector<vector<int> > adj_vertex(m.vertex.size());
    gen_border_graph(m, adj_vertex);
    for (int i = 0; i < m.vertex.size(); ++i)
        for (int j = 0; j < adj_vertex[i].size(); ++j)
            if (adj_vertex[i][j] != -1) {
                vector<int> path;
                euler_simple_circuit(path, adj_vertex, i);
                if (path.size() > max_path.size())
                    max_path.assign(path.begin(), path.end());
            }
    for (list<int>::iterator i = max_path.begin(); i != max_path.end();)
        if (fabs(get_value_on_plane(m.vertex[*i], p)) > EPS)
            i = max_path.erase(i);
        else
            ++i;
}

/*
 * An angle in a polygon may not be convex.
 * This function is to return a normal vector R(orientation vector)
 * If uv X vw is in the same direction of R, then the angle uvw is smaller than 180 degree.
 */
Point3f get_orientation_vector_of_path_on_plane(const Mesh& m, const list<int>& path) {
    Point3f inner_orientation;
    list<int>::const_iterator i = path.begin();
    do {
        list<int>::const_iterator j = prev<list<int> >(path, i);
        list<int>::const_iterator k = next<list<int> >(path, i);
        inner_orientation = cross_product(m.vertex[*i] - m.vertex[*j],
                                          m.vertex[*k] - m.vertex[*i]);
        ++i;
    } while (inner_orientation.length() < EPS);
    i = path.begin();
    float inner_sum = 0;
    while (i != path.end()) {
        list<int>::const_iterator j = prev<list<int> >(path, i);
        list<int>::const_iterator k = next<list<int> >(path, i);
        Point3f&& a = m.vertex[*i] - m.vertex[*j];
        Point3f&& b = m.vertex[*k] - m.vertex[*i];
        float angle = safe_acos(cos(a, b));
        if (cross_product(a, b) * inner_orientation > 0)
            inner_sum += angle;
        else
            inner_sum -= angle;
        ++i;
    }
    // inner_sum must be +2PI or -2PI
    cout << "Inner sum = " << inner_sum << endl;
    if (inner_sum > 0)
        return inner_orientation.normalize();
    else
        return -inner_orientation.normalize();
}

bool is_intersected_with_path(const Mesh& m, const list<int>& path, const Face& f) {
    for (list<int>::const_iterator i = path.begin(); i != path.end(); ++i) {
        list<int>::const_iterator j = next(path, i);
        if (is_intersected(m.vertex[f.vertex_index[0]],
                         m.vertex[f.vertex_index[1]],
                         m.vertex[*i],
                         m.vertex[*j]) &&
            !is_set_intersected(f.vertex_index[0], f.vertex_index[1], *i, *j))
            return true;
        if (is_intersected(m.vertex[f.vertex_index[1]],
                         m.vertex[f.vertex_index[2]],
                         m.vertex[*i], m.vertex[*j]) &&
            !is_set_intersected(f.vertex_index[1], f.vertex_index[2], *i, *j))
            return true;
        if (is_intersected(m.vertex[f.vertex_index[2]],
                         m.vertex[f.vertex_index[0]],
                         m.vertex[*i], m.vertex[*j]) &&
            !is_set_intersected(f.vertex_index[2], f.vertex_index[0], *i, *j))
            return true;
    }
    return false;
}

void advance_front_on_plane(Mesh& m, list<int>& path) {
    if (path.size() < 3)
        return;
    Point3f&& inner_orientation = get_orientation_vector_of_path_on_plane(m, path);
    int last_size = path.size();
    bool success = true;
    while (success) {
        success = false;
        vector<pair<float, list<int>::iterator> > v;
        for (list<int>::iterator i = path.begin(); i != path.end(); ++i) {
            list<int>::iterator j = prev(path, i);
            list<int>::iterator k = next(path, i);
            Point3f a = (m.vertex[*j] - m.vertex[*i]).normalize();
            Point3f b = (m.vertex[*k] - m.vertex[*i]).normalize();
            float cos_value = a * b;
            bool convex = cross_product(-a, b) * inner_orientation > EPS;
            if (convex)
                v.push_back(pair<float, list<int>::iterator>(cos_value, i));
        }
        sort(v.begin(), v.end(), associative_compare<float, list<int>::iterator>);
        // First try
        // Link when angle in (0, 75)
        // Add a new point when angle in (75, 135)
        // Add two new points when angle in (135, 180)
        // If failed(intersected, try to add fewer points.
        for (auto it = v.rbegin(); it != v.rend(); ++it) {
            list<int>::iterator i = it->second;
            list<int>::iterator j = prev(path, i);
            list<int>::iterator k = next(path, i);
            Point3f& u = m.vertex[*j];
            Point3f& v = m.vertex[*i];
            Point3f& w = m.vertex[*k];
            Point3f&& a = u - v;
            Point3f&& b = w - v;
            float la = a.length();
            float lb = b.length();
            a.normalize();
            b.normalize();
            if (it->first < -COS_45) {
                Point3f&& p = v + ((a + a + b).normalize() * (la + la + lb) / 2);
                Point3f&& q = v + ((a + b + b).normalize() * (la + lb + lb) / 2);
                m.vertex.push_back(p);
                m.vertex.push_back(q);
                Face f(*j, *i, m.vertex.size() - 2);
                Face g(m.vertex.size() - 2, *i, m.vertex.size() - 1);
                Face h(m.vertex.size() - 1, *i, *k);
                if (!is_intersected_with_path(m, path, f) &&
                    !is_intersected_with_path(m, path, g) &&
                    !is_intersected_with_path(m, path, h)) {
                    m.face.push_back(f);
                    m.face.push_back(g);
                    m.face.push_back(h);
                    *i = m.vertex.size() - 1;
                    path.insert(i, m.vertex.size() - 2);
                    success = true;
                    break;
                }
            }
            if (it->first < COS_75) {
                Point3f&& p = v + ((a + b).normalize() * (la + lb) / 2);
                m.vertex.push_back(p);
                Face f(*j, *i, m.vertex.size() - 1);
                Face g(m.vertex.size() - 1, *i, *k);
                if (!is_intersected_with_path(m, path, f) &&
                    !is_intersected_with_path(m, path, g)) {
                    m.face.push_back(f);
                    m.face.push_back(g);
                    *i = m.vertex.size() - 1;
                    success = true;
                    break;
                }
            }
            if (it->first < COS_0) {
                Face f(*j, *i, *k);
                if (!is_intersected_with_path(m, path, f)) {
                    m.face.push_back(f);
                    path.erase(i);
                    success = true;
                    break;
                }
            }
        }
    }
}

int parent(vector<int>& p, int x) {
    if (p[x] == x)
        return x;
    else
        return p[x] = parent(p, p[x]);
}

void join(vector<int>& p, int x, int y) {
    p[parent(p, y)] = parent(p, x);
}

void shrink_edge(Mesh& m) {
    vector<float> length;
    for (int i = 0; i < m.vertex.size(); ++i)
        for (int j : m.adj_vertex[i])
            length.push_back((m.vertex[i] - m.vertex[j]).length());
    float avg = average(length);
    float threshold = avg / 2;
    vector<int> p(m.vertex.size());
    for (int i = 0; i < p.size(); ++i)
        p[i] = i;
    vector<bool> border(m.vertex.size(), false);
    gen_border_vertex(m, border);
    for (int i = 0; i < m.vertex.size(); ++i)
        if (!border[i] && parent(p, i) == i)
            for (int j : m.adj_vertex[i])
                if (!border[j] && parent(p, j) == j) {
                    if ((m.vertex[i] - m.vertex[j]).length() < threshold) {
                        merge_vertex(m, i, j);
                        m.vertex[i] = (m.vertex[i] + m.vertex[j]) / 2;
                        join(p, i, j);
                    }
                }
    vector<Face> face;
    for (const Face& f : m.face) {
        if (parent(p, f.vertex_index[0]) != parent(p, f.vertex_index[1]) &&
            parent(p, f.vertex_index[1]) != parent(p, f.vertex_index[2]) &&
            parent(p, f.vertex_index[2]) != parent(p, f.vertex_index[0]))
            face.push_back(f);
    }
    m.face.swap(face);
    m.update();
}

void fill_max_border_face_by_plane(Mesh& m, const Plane& p) {
    Point3f bbox[2];
    memcpy(bbox, m.bbox, sizeof(bbox));
    list<int> path;
    get_max_border_path_by_plane(m, p, path);

    DEBUG()
    // Split edge
    list<float> length;
    for (list<int>::const_iterator i = path.begin(); i != path.end(); ++i) {
        list<int>::const_iterator j = next(path, i);
        length.push_back((m.vertex[*i] - m.vertex[*j]).length());
    }
    float avg = filtered_average(length);
    float sdv = filtered_standard_deviation(length);
    cout << "Avg length = " << avg << endl;
    cout << "Std dev length = " << sdv << endl;
    list<int>::iterator i = path.begin();
    list<float>::iterator k = length.begin();
    for (;i != path.end(); ++i, ++k) {
        list<int>::iterator j = next(path, i);
        list<float>::iterator l = next(length, k);
        if (!is_adjacent(m, *i, *j) && *k > avg + 2 * sdv) {
            int add_vertices = floorf(*k / (avg + sdv));
            Point3f&& unit = (m.vertex[*j] - m.vertex[*i]) / (add_vertices + 1);
            *k = unit.length();
            for (int n = 1; n <= add_vertices; ++n) {
                m.vertex.push_back(m.vertex[*i] + unit * n);
                path.insert(j, m.vertex.size() - 1);
                length.insert(l, unit.length());
            }
        }
    }
    int f = m.face.size();
    // Main AFM
    advance_front_on_plane(m, path);

    Mesh n;
    n.vertex = m.vertex;
    n.texture = m.texture;
    n.face.assign(m.face.begin() + f, m.face.end());
    n.update();
    // Shrink
    shrink_edge(n);
    // Smooth
    laplacian_smooth(n, 20);

    m.vertex = n.vertex;
    m.face.erase(m.face.begin() + f, m.face.end());
    m.face.insert(m.face.end(), n.face.begin(), n.face.end());
    m.update();
    if (m.bbox[0].x[0] < bbox[0].x[0] - EPS ||
        m.bbox[0].x[1] < bbox[0].x[1] - EPS ||
        m.bbox[0].x[2] < bbox[0].x[2] - EPS ||
        m.bbox[1].x[0] > bbox[1].x[0] + EPS ||
        m.bbox[1].x[1] > bbox[1].x[1] + EPS ||
        m.bbox[1].x[2] > bbox[1].x[2] + EPS) {
        cout << "Failed fill face" << endl;
        exit(EXIT_FAILURE);
    }
}


void brute_force_fill_max_border_face_by_plane(Mesh& m, const Plane& p) {
    Point3f bbox[2];
    memcpy(bbox, m.bbox, sizeof(bbox));
    list<int> path;
    get_max_border_path_by_plane(m, p, path);
    Point3f&& inner_orientation = get_orientation_vector_of_path_on_plane(m, path);
    int last_size = path.size();
    while (path.size() > 2) {
        for (list<int>::iterator i = path.begin(); i != path.end(); ) {
            list<int>::iterator j = prev<list<int> >(path, i);
            list<int>::iterator k = next<list<int> >(path, i);
            Point3f& u = m.vertex[*j];
            Point3f& v = m.vertex[*i];
            Point3f& w = m.vertex[*k];
            // uvw is a valid triangle
            // iff uvw < 180 degree
            // and no other vertex on path lies inside uvw
            // and uvw is a triangle
            bool convex = cross_product(v - u, w - v) * inner_orientation > EPS;
            if (!convex) {
                ++i;
                continue;
            }
            bool triangle = cross_product(u - v, w - v).length() > EPS;
            if (!triangle) {
                ++i;
                continue;
            }
            bool inside = false;
            for (list<int>::iterator p = path.begin(); p != path.end(); ++p)
                if (p != i && p != j && p != k && is_inside_triangle(u, v, w, m.vertex[*p])) {
                    inside = true;
                    break;
                }
            if (inside) {
                ++i;
                continue;
            }
            m.face.push_back(Face(*i, *j, *k));
            i = path.erase(i);
        }
        if (path.size() == last_size) {
            cout << "Brute-force fill failed" << endl;
            break;
        }
        last_size = path.size();
    }
    m.update();
    if (m.bbox[0].x[0] < bbox[0].x[0] - EPS ||
        m.bbox[0].x[1] < bbox[0].x[1] - EPS ||
        m.bbox[0].x[2] < bbox[0].x[2] - EPS ||
        m.bbox[1].x[0] > bbox[0].x[0] + EPS ||
        m.bbox[1].x[1] > bbox[0].x[1] + EPS ||
        m.bbox[1].x[2] > bbox[0].x[2] + EPS) {
        cout << "Failed fill face" << endl;
        exit(EXIT_FAILURE);
    }
}

void rotate_point(float rot[3][3], Point3f& p) {
    Point3f np;
    np.x[0] = rot[0][0] * p.x[0] + rot[0][1] * p.x[1] + rot[0][2] * p.x[2];
    np.x[1] = rot[1][0] * p.x[0] + rot[1][1] * p.x[1] + rot[1][2] * p.x[2];
    np.x[2] = rot[2][0] * p.x[0] + rot[2][1] * p.x[1] + rot[2][2] * p.x[2];
    p = np;
}

void rotate_mesh(Mesh& m, const Point3f& axis, float angle) {
    float cosa = cos(angle);
    float sina = sin(angle);
    float onec = 1 - cosa;
    float rot[3][3];
    rot[0][0] = axis.x[0] * axis.x[0] * onec + cosa;
    rot[0][1] = axis.x[0] * axis.x[1] * onec - axis.x[2] * sina;
    rot[0][2] = axis.x[0] * axis.x[2] * onec + axis.x[1] * sina;
    rot[1][0] = axis.x[1] * axis.x[0] * onec + axis.x[2] * sina;
    rot[1][1] = axis.x[1] * axis.x[1] * onec + cosa;
    rot[1][2] = axis.x[1] * axis.x[2] * onec - axis.x[0] * sina;
    rot[2][0] = axis.x[2] * axis.x[0] * onec - axis.x[1] * sina;
    rot[2][1] = axis.x[2] * axis.x[1] * onec + axis.x[0] * sina;
    rot[2][2] = axis.x[2] * axis.x[2] * onec + cosa;
    for (Point3f& p : m.vertex)
        rotate_point(rot, p);
    for (Point3f& p : m.vertex_normal)
        rotate_point(rot, p);
    for (Point3f& p : m.face_normal)
        rotate_point(rot, p);
}

void rotate_mesh(Mesh& m, const Point3f& from, const Point3f& to) {
    // No need to rotate if the total orientation is the same with target destination.
    if (from * to > ONE_EPS)
        return;

    // Reverse
    if (from * to < MINUS_ONE_EPS) {
        for (Point3f& p : m.vertex)
            p = -p;
        return;
    }

    // Rotate around u
    Point3f&& u = cross_product(from, to);
    u.normalize();
    float angle = acos(from * to);
    rotate_mesh(m, u, angle);
}

void auto_rotate_mesh(Mesh& m) {
    Point3f t;
    Point3f z(0, 0, 1);
    for (int i = 0; i < m.face.size(); ++i)
        t += m.face_normal[i];
    t.normalize();
    t.x[0] = 0;
    rotate_mesh(m, t, z);
    rotate_mesh(m,
                Point3f(0, 1, 0),
                Point3f(1, 0, 0));
}

void analyze_z(const Mesh& m) {
    float maxz = -FLT_MAX;
    float minz = FLT_MAX;
    for (const Point3f& p : m.vertex) {
        if (p.x[2] > maxz) maxz = p.x[2];
        if (p.x[2] < minz) minz = p.x[2];
    }
}
