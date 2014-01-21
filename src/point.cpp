#include "point.h"

Point3f operator+(const Point3f& self, const Point3f& other) {
    Point3f ret(self);
    ret += other;
    return ret;
}
Point3f operator+(const Point3f& self, float other) {
    Point3f ret(self);
    ret += other;
    return ret;
}
Point3f operator-(const Point3f& self, const Point3f& other) {
    Point3f ret(self);
    ret -= other;
    return ret;
}
Point3f operator-(const Point3f& self, float other) {
    Point3f ret(self);
    ret -= other;
    return ret;
}
float operator*(const Point3f& self, const Point3f& other) {
    return self.x[0] * other.x[0] + self.x[1] * other.x[1] + self.x[2] * other.x[2];
}
Point3f operator*(const Point3f& self, float other) {
    Point3f ret(self);
    ret *= other;
    return ret;
}
Point3f operator/(const Point3f& self, float other) {
    Point3f ret(self);
    ret /= other;
    return ret;
}
Point3f cross(const Point3f& self, const Point3f& other) {
    Point3f ret;
    ret.x[0] = self.x[1] * other.x[2] - self.x[2] * other.x[1];
    ret.x[1] = self.x[2] * other.x[0] - self.x[0] * other.x[2];
    ret.x[2] = self.x[0] * other.x[1] - self.x[1] * other.x[0];
    return ret;
}
Point2f operator+(const Point2f& self, const Point2f& other) {
    Point2f ret(self);
    ret += other;
    return ret;
}
Point2f operator+(const Point2f& self, float other) {
    Point2f ret(self);
    ret += other;
    return ret;
}
Point2f operator-(const Point2f& self, const Point2f& other) {
    Point2f ret(self);
    ret -= other;
    return ret;
}
Point2f operator-(const Point2f& self, float other) {
    Point2f ret(self);
    ret -= other;
    return ret;
}
float operator*(const Point2f& self, const Point2f& other) {
    return self.x[0] * other.x[0] + self.x[1] * other.x[1];
}
Point2f operator*(const Point2f& self, float other) {
    Point2f ret(self);
    ret *= other;
    return ret;
}
Point2f operator/(const Point2f& self, float other) {
    Point2f ret(self);
    ret /= other;
    return ret;
}
