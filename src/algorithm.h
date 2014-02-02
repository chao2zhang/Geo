#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "base.h"
#include "point.h"
#include "plane.h"
#include "mesh.h"

void laplacian_hc_smooth(Mesh& f, int times=1, float alpha=0.1, float beta=0.6);
void laplacian_smooth(Mesh& f, int times=1, float lambda=0.5);
void center_positioning_by_bounding_box(Mesh& f);
void center_positioning_by_averaging_vertex(Mesh& f);
void partition_by_plane(Mesh& f, const Plane& p);
void project_by_plane(Mesh& f, const Plane& p);
void unify_face_normals(Mesh& f);
void count_spike(Mesh& f, float dot_product=0);
void mark_spike(Mesh& f, float dot_product=0);
void mesh_offset(Mesh& f, float offset);
void remove_face_by_plane(Mesh& f, const Plane& p);
void fill_max_border_face_by_plane(Mesh& f, const Plane& p);

#endif
