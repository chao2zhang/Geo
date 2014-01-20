#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "object.h"

void laplacian_hc_smooth(Object& obj, int times=1, float alpha=0, float beta=0.5);
void center_positioning(Object& o);
void partition_by_plane(Object& o, float a, float b, float c, float d);

#endif
