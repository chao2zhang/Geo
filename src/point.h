#ifndef POINT_H
#define POINT_H

#include <cstring>
#include <cmath>

#include "base.h"

class Point3f;

inline Point3f operator+(const Point3f& self, const Point3f& other);
inline Point3f operator+(float other, const Point3f& self);
inline Point3f operator+(const Point3f& self, float other);
inline Point3f operator-(const Point3f& self, const Point3f& other);
inline Point3f operator-(const Point3f& self, float other);
inline Point3f operator-(const Point3f& self);
inline float operator*(const Point3f& self, const Point3f& other);
inline Point3f operator*(const Point3f& self, float other);
inline Point3f operator*(float other, const Point3f& self);
inline Point3f operator/(const Point3f& self, float other);
inline Point3f cross_product(const Point3f& self, const Point3f& other);

class Point3f {
public:
    float x[3];
    Point3f() {
        x[0] = 0;
        x[1] = 0;
        x[2] = 0;
    }
    Point3f(float a, float b, float c) {
        x[0] = a;
        x[1] = b;
        x[2] = c;
    }
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
    float length_square() const {
        return *this * *this;
    }
    Point3f& normalize() {
        float l = length();
        if (l == 0.0f)
            return *this;
        x[0] /= l;
        x[1] /= l;
        x[2] /= l;
        return *this;
    }
    ~Point3f() {
    }
};

inline Point3f operator+(const Point3f& self, const Point3f& other) {
    Point3f ret(self);
    ret += other;
    return ret;
}
inline Point3f operator+(const Point3f& self, float other) {
    Point3f ret(self);
    ret += other;
    return ret;
}
inline Point3f operator+(float other, const Point3f& self) {
    return self + other;
}
inline Point3f operator-(const Point3f& self, const Point3f& other) {
    Point3f ret(self);
    ret -= other;
    return ret;
}
inline Point3f operator-(const Point3f& self, float other) {
    Point3f ret(self);
    ret -= other;
    return ret;
}
inline Point3f operator-(const Point3f& self) {
    Point3f ret;
    ret.x[0] = -self.x[0];
    ret.x[1] = -self.x[1];
    ret.x[2] = -self.x[2];
    return ret;
}
inline float operator*(const Point3f& self, const Point3f& other) {
    return self.x[0] * other.x[0] + self.x[1] * other.x[1] + self.x[2] * other.x[2];
}
inline Point3f operator*(float other, const Point3f& self) {
    return self * other;
}
inline Point3f operator*(const Point3f& self, float other) {
    Point3f ret(self);
    ret *= other;
    return ret;
}
inline Point3f operator/(const Point3f& self, float other) {
    Point3f ret(self);
    ret /= other;
    return ret;
}
inline Point3f cross_product(const Point3f& self, const Point3f& other) {
    Point3f ret;
    ret.x[0] = self.x[1] * other.x[2] - self.x[2] * other.x[1];
    ret.x[1] = self.x[2] * other.x[0] - self.x[0] * other.x[2];
    ret.x[2] = self.x[0] * other.x[1] - self.x[1] * other.x[0];
    return ret;
}

inline float cos(const Point3f& self, const Point3f& other) {
    return self * other / self.length() / other.length();
}

class Point2f;

inline Point2f operator+(const Point2f& self, const Point2f& other);
inline Point2f operator+(const Point2f& self, float other);
inline Point2f operator+(float other, const Point2f& self);
inline Point2f operator-(const Point2f& self, const Point2f& other);
inline Point2f operator-(const Point2f& self, float other);
inline Point2f operator-(const Point2f& self);
inline float operator*(const Point2f& self, const Point2f& other);
inline Point2f operator*(const Point2f& self, float other);
inline Point2f operator*(float other, const Point2f& self);
inline Point2f operator/(const Point2f& self, float other);

class Point2f {
public:
    float x[2];
    Point2f(float a, float b) {
        x[0] = a;
        x[1] = b;
    }
    Point2f() {
        x[0] = 0;
        x[1] = 0;
    }
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
    float length_square() const {
        return *this * *this;
    }
    Point2f& normalize() {
        float l = length();
        if (l == 0.0f)
            return *this;
        x[0] /= l;
        x[1] /= l;
        return *this;
    }
    ~Point2f() {
    }
};

inline Point2f operator+(const Point2f& self, const Point2f& other) {
    Point2f ret(self);
    ret += other;
    return ret;
}
inline Point2f operator+(const Point2f& self, float other) {
    Point2f ret(self);
    ret += other;
    return ret;
}
inline Point2f operator+(float other, const Point2f& self) {
    return self + other;
}
inline Point2f operator-(const Point2f& self, const Point2f& other) {
    Point2f ret(self);
    ret -= other;
    return ret;
}
inline Point2f operator-(const Point2f& self, float other) {
    Point2f ret(self);
    ret -= other;
    return ret;
}
inline Point2f operator-(const Point2f& self) {
    Point2f ret;
    ret.x[0] = self.x[0];
    ret.x[1] = self.x[1];
    return ret;
}
inline float operator*(const Point2f& self, const Point2f& other) {
    return self.x[0] * other.x[0] + self.x[1] * other.x[1];
}
inline Point2f operator*(const Point2f& self, float other) {
    Point2f ret(self);
    ret *= other;
    return ret;
}
inline Point2f operator*(float other, const Point2f& self) {
    return self * other;
}
inline Point2f operator/(const Point2f& self, float other) {
    Point2f ret(self);
    ret /= other;
    return ret;
}

inline float cos(const Point2f& self, const Point2f& other) {
    return self * other / self.length() / other.length();
}

#endif
