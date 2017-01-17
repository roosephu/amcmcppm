#include "vec3.h"
#include "mcmc.h"

struct Light {
    Vec3 flux;
    Light() {}

    virtual Ray gen_ray(PathSpace &path) = 0;
};

struct PointLight : Light {
    Vec3 loc;
    real rad;

    PointLight(Vec3 loc_, Vec3 flux_, real rad_)
        : loc(loc_), rad(rad_)
    {
        flux = flux_;
    }

    Ray gen_ray(PathSpace &path) {
        Vec3 dir = gen_sphere(path);
        return Ray(loc + dir * rad, dir);
    }

};

struct SemisphereLight : Light {
    Vec3 dir, loc;
    real rad;

    SemisphereLight(Vec3 loc_, Vec3 dir_, Vec3 flux_, real rad_)
        : loc(loc_), dir(dir_), rad(rad_)
    {
        flux = flux_;
    }

    Ray gen_ray(PathSpace &path) {
        Vec3 d = gen_sphere(path);
        if (dot(dir, d) <= 0)
            d = d * -1;
        return Ray(loc + d * rad, d);
    }

};

struct PlaneLight : Light {
    Vec3 loc, dir, cx, cy;
    real dist;

    PlaneLight(Vec3 loc_, Vec3 dir_, real dist_, Vec3 flux_)
        : loc(loc_), dir(dir_), dist(dist_)
    {
        flux = flux_;
        if (abs(dir.x) < 0.5)
            cx = Vec3(1, 0, 0);
        else if (abs(dir.y) < 0.5)
            cx = Vec3(0, 1, 0);
        else
            cx = Vec3(0, 0, 1);
        cy = dir % cx;
    }

    Ray gen_ray(PathSpace &path) {
        real u1 = 2 * path.next() - 1, u2 = 2 * path.next() - 1;
        return Ray(loc + cx * u1 + cy * u2, dir);
    }
};
