#ifndef PLANE_H
#define PLANE_H

class Plane {
public:
    float a, b, c, d;
    Plane(float a_, float b_, float c_, float d_): a(a_), b(b_), c(c_), d(d_){}
    Plane(): a(0), b(0), c(0), d(0){}
};

#endif // PLANE_H
