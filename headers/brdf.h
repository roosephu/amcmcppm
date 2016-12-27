#pragma once
#include "mcmc.h"
#include "vec3.h"
#include <algorithm>
using namespace std;

struct BRDF {
    bool is_delta;

    virtual pair<Vec3, Vec3> sample(PathSpace &path, Vec3 dir, Vec3 normal, bool into) = 0;

    virtual BRDF *get_brdf(PathSpace &path) {
        return this;
    }
};

struct Refraction : BRDF {
    Vec3 weight;
    real ratio, decay;

    Refraction(Vec3 w_, real r_, real d_=0.1) : weight(w_), ratio(r_), decay(d_) {
        is_delta = true;
    }

    pair<Vec3, Vec3> sample(PathSpace &path, Vec3 dir, Vec3 normal, bool into) {
        real in_dot_normal = dot(normal, dir);
        real real_ratio = into ? 1 / ratio : ratio;

        real sin_old = catheti(1, in_dot_normal);
        real sin_new = sin_old * real_ratio;
        Vec3 ref = (dir - normal * (2 * in_dot_normal)).normalize();

        if (sin_new > 1) { // only reflexion
            return {ref, Vec3(1, 1, 1)};
        } else {
            real cos_new = catheti(1, sin_new), cos_old = -in_dot_normal;

            Vec3 out = (dir * real_ratio - normal * (real_ratio * in_dot_normal + cos_new)).normalize();

            real rs = (cos_old - real_ratio * cos_new) / (cos_old + real_ratio * cos_new);
            real rp = (real_ratio * cos_old - cos_new) / (real_ratio * cos_old + cos_new);
            real prob_reflexion = (real) (rs * rs + rp * rp) * 0.5;

            // assert(abs(dot(dir % normal, out)) < 1e-6);

            // real R0 = sqr((real_ratio - 1) / (real_ratio + 1));
            // real sub = 1 + (into ? in_dot_normal : dot(normal, out));
            // real R = R0 + (1 - R0) * pow(sub, 5);
            // if (!into)
            //     printf("%d %.6f %.6f\n", (int) into, R, prob_reflexion);

            if (path.next() <= prob_reflexion) {
                return {ref, weight};
            } else {
                return {out, weight};
            }
        }
    }
};

struct Reflexion : BRDF {
    Vec3 weight;

    Reflexion(Vec3 w_) : weight(w_) {
        is_delta = true;
    }

    pair<Vec3, Vec3> sample(PathSpace &path, Vec3 dir, Vec3 normal, bool into) {
        Vec3 ref = (dir - normal * (2 * dot(normal, dir))).normalize();
        return {ref, weight};
    }
};

struct Diffusion : BRDF {
    Vec3 weight;

    Diffusion(Vec3 w_) : weight(w_) {
        is_delta = false;
    }

    pair<Vec3, Vec3> sample(PathSpace &path, Vec3 dir, Vec3 normal, bool into) {
        Vec3 out = gen_sphere(path);
        if (dot(out, normal) < 0)
            out = out * -1;
        return {out, weight / pi};
    }
};

struct MixedBRDF : BRDF {
    vector<pair<BRDF *, real>> brdf;
    real sum;

    MixedBRDF(vector<pair<BRDF *, real>> input) {
        brdf = input;
        sum = 0;
        for (auto it : brdf) {
            sum += it.second;
        }
    }

    BRDF *get_brdf(PathSpace &path) {
        real x = path.next() * sum;
        for (auto &it : brdf) {
            if (x < it.second) {
                return it.first->get_brdf(path);
            }
            x -= it.second;
        }
        return nullptr;
    }

    pair<Vec3, Vec3> sample(PathSpace &path, Vec3 dir, Vec3 normal, bool into) {
        return get_brdf(path)->sample(path, dir, normal, into);
    }
};

struct Blur : BRDF {
    BRDF *brdf;

    Blur(BRDF *b_) : brdf(b_)
    {
        is_delta = false;
    }

    pair<Vec3, Vec3> sample(PathSpace &path, Vec3 dir, Vec3 normal, bool into) {
        auto ret = brdf->get_brdf(path)->sample(path, dir, normal, into);
        ret.first = (ret.first + Vec3(path.next(), path.next(), path.next()) * 0.1).normalize();
        return ret;
    }
};