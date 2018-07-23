#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace std;

typedef double real;

const real infinity = 1e30;
const real pi = acos(-1);

struct Vec3 {
    real x, y, z;

    Vec3(real _x=0, real _y=0, real _z=0) :
        x(_x), y(_y), z(_z) {}

    real len() { return sqrt(x * x + y * y + z * z); }
    Vec3 operator + (Vec3 b) { return Vec3(x + b.x, y + b.y, z + b.z); }
    Vec3 operator - (Vec3 b) { return Vec3(x - b.x, y - b.y, z - b.z); }
    Vec3 operator * (Vec3 b) { return Vec3(x * b.x, y * b.y, z * b.z); }
    Vec3 operator / (Vec3 b) { return Vec3(x / b.x, y / b.y, z / b.z); }
    Vec3 operator + (real l) { return Vec3(x + l, y + l, z + l); }
    Vec3 operator - (real l) { return Vec3(x - l, y - l, z - l); }
    Vec3 operator * (real l) { return Vec3(x * l, y * l, z * l); }
    Vec3 operator / (real l) { return Vec3(x / l, y / l, z / l); }

    Vec3 operator % (Vec3 b) {
        return Vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }

    real max() {
        return std::max(std::max(x, y), z);
    }

    real min() {
        return std::min(std::min(x, y), z);
    }

    real sum() {
        return x + y + z;
    }

    Vec3 normalize() {
        return *this / this->len();
    }
};

typedef Vec3 Color;

struct Ray {
    Vec3 loc, dir;

    Ray(Vec3 loc_, Vec3 dir_) :
        loc(loc_), dir(dir_.normalize()) {}
};

Vec3 min(Vec3 a, Vec3 b) {
    return Vec3(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}

Vec3 max(Vec3 a, Vec3 b) {
    return Vec3(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}

real dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

real sqr(real t) {
    return t * t;
}

real catheti(real h, real x) {
    if (h < x) {
        if (x < h + 1e-6)
            return 0;
        printf("catheti: h = %.6f, x = %.6f\n", h, x);
        return 0;
    }
    return sqrt(h * h - x * x);
}
