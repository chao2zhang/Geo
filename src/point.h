#ifndef POINT_H
#define POINT_H

#include <cstring>
#include <cmath>

using namespace std;

class Point3f;

const float eps = 0.000001;

Point3f operator+(const Point3f& self, const Point3f& other);
Point3f operator+(const Point3f& self, float other);
Point3f operator-(const Point3f& self, const Point3f& other);
Point3f operator-(const Point3f& self, float other);
float operator*(const Point3f& self, const Point3f& other);
Point3f operator*(const Point3f& self, float other);
Point3f operator/(const Point3f& self, float other);
Point3f cross(const Point3f& self, const Point3f& other);

class Point3f {
public:
    float x[3];
    Point3f(float a, float b, float c) {
        x[0] = a;
        x[1] = b;
        x[2] = c;
    }
    Point3f():Point3f(0, 0, 0) {}
    Point3f(const Point3f& other) {
        memcpy(x, other.x, sizeof(float) * 3);
    }
    Point3f& operator=(const Point3f& other) {
        if (&other == this)
            return *this;
        memcpy(x, other.x, sizeof(float) * 3);
        return *this;
    }
    Point3f& operator+=(const Point3f& other) {
        x[0] += other.x[0];
        x[1] += other.x[1];
        x[2] += other.x[2];
        return *this;
    }
    Point3f& operator+=(float other) {
        x[0] += other;
        x[1] += other;
        x[2] += other;
        return *this;
    }
    Point3f& operator-=(const Point3f& other) {
        x[0] -= other.x[0];
        x[1] -= other.x[1];
        x[2] -= other.x[2];
        return *this;
    }
    Point3f& operator-=(float other) {
        x[0] -= other;
        x[1] -= other;
        x[2] -= other;
        return *this;
    }
    Point3f& operator*=(float other) {
        x[0] *= other;
        x[1] *= other;
        x[2] *= other;
        return *this;
    }
    Point3f& operator/=(float other) {
        x[0] /= other;
        x[1] /= other;
        x[2] /= other;
        return *this;
    }
    float length() const {
        return sqrt(*this * *this);
    }
    Point3f& normalize() {
        float l = length();
        x[0] /= l;
        x[1] /= l;
        x[2] /= l;
        return *this;
    }
    ~Point3f() {
    }
};

class Point2f;

Point2f operator+(const Point2f& self, const Point2f& other);
Point2f operator+(const Point2f& self, float other);
Point2f operator-(const Point2f& self, const Point2f& other);
Point2f operator-(const Point2f& self, float other);
float operator*(const Point2f& self, const Point2f& other);
Point2f operator*(const Point2f& self, float other);
Point2f operator/(const Point2f& self, float other);

class Point2f {
public:
    float x[2];
    Point2f(float a, float b) {
        x[0] = a;
        x[1] = b;
    }
    Point2f():Point2f(0, 0) {}
    Point2f(const Point2f& other) {
        memcpy(x, other.x, sizeof(float) * 2);
    }
    Point2f& operator=(const Point2f& other) {
        if (&other == this)
            return *this;
        memcpy(x, other.x, sizeof(float) * 2);
        return *this;
    }
    Point2f& operator+=(const Point2f& other) {
        x[0] += other.x[0];
        x[1] += other.x[1];
        return *this;
    }
    Point2f& operator+=(float other) {
        x[0] += other;
        x[1] += other;
        return *this;
    }
    Point2f& operator-=(const Point2f& other) {
        x[0] -= other.x[0];
        x[1] -= other.x[1];
        return *this;
    }
    Point2f& operator-=(float other) {
        x[0] -= other;
        x[1] -= other;
        return *this;
    }
    Point2f& operator*=(float other) {
        x[0] *= other;
        x[1] *= other;
        return *this;
    }
    Point2f& operator/=(float other) {
        x[0] /= other;
        x[1] /= other;
        return *this;
    }
    float length() const {
        return sqrt(*this * *this);
    }
    Point2f& normalize() {
        float l = length();
        x[0] /= l;
        x[1] /= l;
        return *this;
    }
    ~Point2f() {
    }
};


class Plane {
public:
    float a, b, c, d;
    Plane(float a_, float b_, float c_, float d_) {
        a = a_;
        b = b_;
        c = c_;
        d = d_;
    }
    Plane():Plane(0, 0, 0, 0) {}
};


#endif
