#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "base.h"
#include "point.h"
#include "plane.h"
#include "object.h"

void laplacian_hc_smooth(Object& obj, int times=1, float alpha=0, float beta=0.5);
void center_positioning(Object& o);
void partition_by_plane(Object& o, const Plane& p);
void unify_face_normals(Object& o);
void count_spike(Object& o, float dot_product=0);
void clear_spike(Object& o, float dot_product=0);
void mesh_offset(Object& o, float offset);

#endif
