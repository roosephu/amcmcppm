#include "vec3.h"
#include "brdf.h"

const real intersect_threshold = 1e-6;
const real triangle_extend = -1e-7;
const int MIN_SIZE = 5;

struct Intersection {
    Vec3 loc, normal;
};

struct Object {
    BRDF *brdf;
    int id;

    Object() {}

    virtual pair<real, real> projection(int d) = 0;
    virtual Vec3 lower() = 0;
    virtual Vec3 upper() = 0;
    virtual real intersect(Ray *ray) = 0;
    virtual Vec3 get_normal(Vec3 intersection) = 0;
    virtual BRDF *get_brdf(PathSpace &path) {
        return brdf->get_brdf(path);
    }
};

struct Sphere : Object {
    Vec3 loc;
    real rad;

    Sphere(Vec3 loc_, real r_, BRDF *brdf_) :
        loc(loc_), rad(r_)
    {
        brdf = brdf_;
    }

    real intersect(Ray *r) {
        real proj = dot(loc - r->loc, r->dir);
        real dist_center = dot(loc - r->loc, loc - r->loc);
        // real dist = sqrt(dist_center - proj * proj);
        real r_proj = rad * rad - dist_center + proj * proj;
        if (r_proj < 0)
            return infinity;
        r_proj = sqrt(r_proj);
        if (proj - r_proj > intersect_threshold)
            return proj - r_proj;
        if (proj + r_proj > intersect_threshold)
            return proj + r_proj;
        return infinity;
    }

    pair<real, real> projection(int d) {
        if (d == 0)
            return {loc.x - rad, loc.x + rad};
        else if (d == 1)
            return {loc.y - rad, loc.y + rad};
        else
            return {loc.z - rad, loc.z + rad};
    }

    Vec3 lower() {
        return loc - rad;
    }

    Vec3 upper() {
        return loc + rad;
    }

    Vec3 get_normal(Vec3 intersection) {
        return (intersection - loc).normalize();
    }
};

struct Triangle : Object {
    Vec3 x, y, z, normal, L, R;
    double area;

    Triangle(Vec3 x_, Vec3 y_, Vec3 z_, BRDF *brdf_) : x(x_), y(y_), z(z_) {
        Vec3 v = (y - x) % (z - x);
        area = v.len();
        normal = v / area;

        L = min(min(x, y), z);
        R = max(max(x, y), z);
        brdf = brdf_;
    }

    real intersect(Ray *ray) {
        if (area < 1e-7)
            return infinity;
        real t = dot(x - ray->loc, normal) / dot(ray->dir, normal);

        if (t < intersect_threshold)
            return infinity;
        Vec3 v = ray->loc + ray->dir * t;
        real a = dot((x - v) % (y - v), normal);
        real b = dot((y - v) % (z - v), normal);
        real c = dot((z - v) % (x - v), normal);

        if ((a >= -triangle_extend) == (b >= -triangle_extend) && (b >= -triangle_extend) == (c >= -triangle_extend))
            return t;
        return infinity;
    }

    Vec3 lower() {
        return L;
    }

    Vec3 upper() {
        return R;
    }

    Vec3 get_normal(Vec3 intersection){
        return normal;
    }

    pair<real, real> projection(int d) {
        if (d == 0)
            return {L.x, R.x};
        else if (d == 1)
            return {L.y, R.y};
        else
            return {L.z, R.z};
    }
};

struct TriangleNormal : Object {
    Vec3 nx, ny, nz;
    real area;

    Vec3 x, y, z, normal, L, R;

    TriangleNormal(pair<Vec3, Vec3> x_, pair<Vec3, Vec3> y_, pair<Vec3, Vec3> z_, BRDF *brdf_) {
        x = x_.first, nx = x_.second;
        y = y_.first, ny = y_.second;
        z = z_.first, nz = z_.second;
        Vec3 v = (y - x) % (z - x);
        area = v.len();
        // printf("%.8f\n", area);
        normal = v / area;
        // printf("%.3f\n", normal.y);

        L = min(min(x, y), z);
        R = max(max(x, y), z);
        brdf = brdf_;
    }

    Vec3 lower() {
        return L;
    }

    Vec3 upper() {
        return R;
    }

    pair<real, real> projection(int d) {
        if (d == 0)
            return {L.x, R.x};
        else if (d == 1)
            return {L.y, R.y};
        else
            return {L.z, R.z};
    }

    real intersect(Ray *ray) {
        if (area < 1e-7)
            return infinity;
        real t = dot(x - ray->loc, normal) / dot(ray->dir, normal);

        if (t < intersect_threshold)
            return infinity;
        Vec3 v = ray->loc + ray->dir * t;
        real a = dot((x - v) % (y - v), normal);
        real b = dot((y - v) % (z - v), normal);
        real c = dot((z - v) % (x - v), normal);

        if ((a >= -triangle_extend) == (b >= -triangle_extend) && (b >= -triangle_extend) == (c >= -triangle_extend))
            return t;
        return infinity;
    }

    Vec3 get_normal(Vec3 v){
        real a = dot((x - v) % (y - v), normal);
        real b = dot((y - v) % (z - v), normal);
        real c = dot((z - v) % (x - v), normal);
        auto ret = (nz * a + nx * b + ny * c) / area;
        assert(abs(a + b + c - area) <= 1e-7);

        if (debug) {
            printf("\narea = %.8f, x = (%f %f %f) y = (%f %f %f), z = (%f %f %f)\n", area, x.x, x.y, x.z, y.x, y.y, y.z, z.x, z.y, z.z);
            printf("nx = (%f %f %f), ny = (%f %f %f), nz = (%f %f %f)\n", nx.x, nx.y, nx.z, ny.x, ny.y, ny.z, nz.x, nz.y, nz.z);
        }
        return ret;
    }
    
};
