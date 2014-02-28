#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "base.h"
#include "point.h"
#include "plane.h"
#include "mesh.h"

void laplacian_hc_smooth(Mesh& m, int times=1, bool reserver_border=true, float alpha=0.1, float beta=0.6);
void laplacian_smooth(Mesh& m, int times=1, bool reserve_border=true, float lambda=0.5);
void center_positioning_by_bounding_box(Mesh& m);
void center_positioning_by_averaging_vertex(Mesh& m);
void partition_by_plane(Mesh& m, const Plane& p);
void project_by_plane(Mesh& m, const Plane& p);
void unify_face_normals(Mesh& m);
void count_spike(Mesh& m, float dot_product=0);
void mark_spike(Mesh& m, float dot_product=0);
void mesh_offset(Mesh& m, float offset);
void remove_face_by_plane(Mesh& m, const Plane& p);
void remove_face_by_largest_component(Mesh& m);
void remove_face_by_long_edge(Mesh& m, float threshold);
void fill_max_border_face_by_plane(Mesh& m, const Plane& p);
void brute_force_fill_max_border_face_by_plane(Mesh& m, const Plane& p);
void auto_rotate_mesh(Mesh& m);
void rotate_mesh(Mesh& m, const Point3f& axis, float angle);
void rotate_mesh(Mesh& m, const Point3f& from, const Point3f& to);
void analyze_z(const Mesh& m);

#endif
