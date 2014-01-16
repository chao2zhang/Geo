#ifndef POINT_H
#define POINT_H

#include <cstdio>
#include <cstdlib>

using namespace std;

class Point3f {
public:
    float x[3];
    Point3f(float a, float b, float c) {
        x[0] = a;
        x[1] = b;
        x[2] = c;
    }
    Point3f() {
        Point3f(0, 0, 0);
    }
    Point3f(const Point3f& other) {
        memcpy(x, other.x, sizeof(float) * 3);
    }
    Point3f& operator=(const Point3f& other) {
        if (&other == this)
            return;
        memcpy(x, other.x, sizeof(float) * 3);
    }
    Point3f& operator+=(const Point3f& other) {
        x[0] += other.x[0];
        x[1] += other.x[1];
        x[2] += other.x[2];
    }
    Point3f& operator+=(float other) {
        x[0] += other;
        x[1] += other;
        x[2] += other;
    }
    Point3f& operator-=(const Point3f& other) {
        x[0] -= other.x[0];
        x[1] -= other.x[1];
        x[2] -= other.x[2];
    }
    Point3f& operator-=(float other) {
        x[0] -= other;
        x[1] -= other;
        x[2] -= other;
    }
    Point3f& operator*=(float other) {
        x[0] *= other;
        x[1] *= other;
        x[2] *= other;
    }
    Point3f& operator/=(float other) {
        x[0] /= other;
        x[1] /= other;
        x[2] /= other;
    }
    ~Point3f() {
    }
};

inline Point3f operator+(const Point3f& self, const Point3f& other) {
    Point3f ret(self);
    self += other;
    return ret;
}
inline Point3f operator+(const Point3f& self, float other) {
    Point3f ret(self);
    self += other;
    return ret;
}
inline Point3f operator-(const Point3f& self, const Point3f& other) {
    Point3f ret(self);
    self -= other;
    return ret;
}
inline Point3f operator-(const Point3f& self, float other) {
    Point3f ret(self);
    self -= other;
    return ret;
}
inline float operator*(const Point3f& self, const Point3f& other) {
    return self.x[0] * other.x[0] + self.x[1] * other.x[1] + self.x[2] * other.x[2];
}
inline Point3f operator/(const Point3f& self, float other) {
    Point3f ret(self);
    self /= other;
    return ret;
}

class Point2f {
public:
    float x[2];
    Point2f(float a, float b) {
        x[0] = a;
        x[1] = b;
    }
    Point2f() {
        Point2f(0, 0);
    }
    Point2f(const Point2f& other) {
        memcpy(x, other.x, sizeof(float) * 2);
    }
    Point2f& operator=(const Point2f& other) {
        if (&other == this)
            return;
        memcpy(x, other.x, sizeof(float) * 2);
    }
    Point2f& operator+=(const Point2f& other) {
        x[0] += other.x[0];
        x[1] += other.x[1];
    }
    Point2f& operator+=(float other) {
        x[0] += other;
        x[1] += other;
    }
    Point2f& operator-=(const Point2f& other) {
        x[0] -= other.x[0];
        x[1] -= other.x[1];
    }
    Point2f& operator-=(float other) {
        x[0] -= other;
        x[1] -= other;
    }
    Point2f& operator*=(float other) {
        x[0] *= other;
        x[1] *= other;
    }
    Point2f& operator/=(float other) {
        x[0] /= other;
        x[1] /= other;
    }
    ~Point2f() {
    }
};

inline Point2f operator+(const Point2f& self, const Point2f& other) {
    Point2f ret(self);
    self += other;
    return ret;
}
inline Point2f operator+(const Point2f& self, float other) {
    Point2f ret(self);
    self += other;
    return ret;
}
inline Point2f operator-(const Point2f& self, const Point2f& other) {
    Point2f ret(self);
    self -= other;
    return ret;
}
inline Point2f operator-(const Point2f& self, float other) {
    Point2f ret(self);
    self -= other;
    return ret;
}
inline float operator*(const Point2f& self, const Point2f& other) {
    return self.x[0] * other.x[0] + self.x[1] * other.x[1];
}
inline Point2f operator*(const Point2f& self, float other) {
    Point2f ret(self);
    self *= other;
    return ret;
}
inline Point2f operator/(const Point2f& self, float other) {
    Point2f ret(self);
    self /= other;
    return ret;
}

#endif