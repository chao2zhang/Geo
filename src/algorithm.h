#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "base.h"
#include "point.h"
#include "plane.h"
#include "mesh.h"
#include <vector>
#include <list>

void laplacian_hc_smooth(Mesh& m, int times=1, bool reserver_border=true, float alpha=0.1, float beta=0.6);
void laplacian_smooth(Mesh& m, int times=1, bool reserve_border=true, float lambda=0.5);
void center_positioning_by_bounding_box(Mesh& m);
void center_positioning_by_averaging_vertex(Mesh& m);
void partition_by_plane(Mesh& m, const Plane& p);
void project_by_z_plane(Mesh& m, float z);
void unify_face_normals(Mesh& m);
void count_spike(Mesh& m, float dot_product=0);
void mark_spike(Mesh& m, float dot_product=0);
void mesh_offset(Mesh& m, float offset);
void remove_face_by_y_plane(Mesh& m, float y);
void remove_face_by_largest_component(Mesh& m);
void remove_face_by_long_edge(Mesh& m, float rate=5);
void fill_max_border_face_by_plane(Mesh& m, const Plane& p);
void brute_force_fill_max_border_face_by_plane(Mesh& m, const Plane& p);
void auto_rotate_mesh(Mesh& m);
void rotate_mesh(Mesh& m, const Point3f& axis, float angle);
void rotate_mesh(Mesh& m, const Point3f& from, const Point3f& to);

void reverse_face_orientation(Mesh& m, int f);
bool correct_winding(Mesh& m, int f, int g);
void merge_vertex(Mesh& m, int u, int v);
void mark_face_of_edge(const Mesh& m, int u, int v, std::vector<bool>& mark);
void gen_border_vertex(const Mesh& m, std::vector<bool>& is_border);
void gen_border_graph(const Mesh& m, std::vector<std::vector<int> >& adj_vertex);
void brute_force_fill_path(Mesh& m, const std::list<int>& p);
bool euler_simple_circuit(std::vector<int>& path, std::vector<bool>& visited, std::vector<std::vector<int> >& adj_vertex, int s, int u);
void euler_simple_circuit(std::vector<int>& path, std::vector<std::vector<int> >& adj_vertex, int v);
void gen_offset_vertex(const Mesh& m, float offset, Mesh& n);
void gen_offset_invalid_vertex(const Mesh& m, const Mesh& n, float offset, std::vector<bool>& is_invalid);
void gen_offset_face(const Mesh& m, Mesh& n, const std::vector<bool>& is_invalid);
void enclose_offset_mesh(Mesh& m, Mesh& n, const std::vector<bool>& is_invalid);
void detect_self_intersect(const Mesh& m, std::vector<bool>& is_int);
int flood_fill_face_group(const Mesh& m, std::vector<int>& face_group);
Point3f get_orientation_vector_of_path_on_plane(const Mesh& m, const std::list<int>& path);
bool is_intersected_with_path(const Mesh& m, const std::list<int>& path, const Face& f);
void advance_front_on_plane(Mesh& m, std::list<int>& path);
int parent(std::vector<int>& p, int x);
void join(std::vector<int>& p, int x, int y);
void shrink_edge(Mesh& m);
void rotate_point(float rot[3][3], Point3f& p);

#endif
