#pragma once
#include "vec3.h"
#include "mcmc.h"
#include "brdf.h"
#include "mcmc.h"
#include "object.h"
#include <algorithm>
#include <cassert>
using namespace std;

// const real eps = 1e-9;

// PathSpace rng(false);

const int N = 1000000;

struct bvh_node {
    bvh_node *lc, *rc;
    int dir;
    Vec3 L, R;
    Object *obj[MIN_SIZE];
    int size;

    bvh_node() {}

    real intersect(Ray *ray) {
        Vec3 tl = (L - ray->loc) / ray->dir, tr = (R - ray->loc) / ray->dir;
        real enter = min(tl, tr).max(), leave = max(tl, tr).min();
        if (enter > leave) {
            return infinity;
        }
        return max(enter, (real) 0.);
    }
};

struct BVH {
    bvh_node *root;
    real hit;
    Object *target;
    Object *objects[N];
    int n;
    long count, count_intersect;

    BVH() : n(0) {
    }

    void add(Object *obj) {
        assert(n < N);
        objects[n] = obj;
        n += 1;
    }

    void build() {
        root = build(objects, n);
        printf("[bvh]: finished. %d objects\n", n);
    }

    bvh_node *build(Object **objs, int n) {
        // if (n > 10000) {
        //     printf("[bvh] constructing: %d\n", n);
        // }
        bvh_node *ret = new bvh_node();
        ret->size = n;

        ret->L = Vec3(1, 1, 1) * +infinity;
        ret->R = Vec3(1, 1, 1) * -infinity;
        for (int i = 0; i < n; ++i) {
            ret->L = min(ret->L, objs[i]->lower());
            ret->R = max(ret->R, objs[i]->upper());
        }
        real area = ret->R.sum() - ret->L.sum();

        double min_cost = infinity;
        int direction, position;
        vector<tuple<real, Object *, real>> coord;
        coord.reserve(n);

        for (int d = 0; d < 3; ++d) {
            coord.clear();
            for (int i = 0; i < n; ++i)
                coord.push_back((tuple<real, Object *, real>) {objs[i]->projection(d).second, objs[i], (real) 1});
            sort(coord.begin(), coord.end());

            for (int _ = 0; _ < 2; ++_) {
                Vec3 L = Vec3(1, 1, 1) * +infinity;
                Vec3 R = Vec3(1, 1, 1) * -infinity;
                for (int i = 0; i < n; ++i) {
                    Object *obj = get<1>(coord[i]); 
                    L = min(L, obj->lower());
                    R = max(R, obj->upper());
                    get<2>(coord[i]) += (R.sum() - L.sum()) / area * (i + 1);
                }
                reverse(coord.begin(), coord.end());
            }

            for (int i = 0; i < n - 1; ++i) {
                real cost = get<2>(coord[i]);
                if (cost < min_cost) {
                    min_cost = cost;
                    direction = d;
                    position = i + 1;
                }
            }
        }

        if (n <= MIN_SIZE) {
            for (int i = 0; i < n; ++i)
                ret->obj[i] = objs[i];
            ret->lc = ret->rc = nullptr;
        } else {
            coord.clear();
            for (int i = 0; i < n; ++i)
                coord.push_back((tuple<real, Object *, real>) {objs[i]->projection(direction).second, objs[i], (real) 0});
            sort(coord.begin(), coord.end());
            for (int i = 0; i < n; ++i) {
                objs[i] = get<1>(coord[i]);
            }

            ret->dir = direction;
            ret->lc = build(objs           ,     position);
            ret->rc = build(objs + position, n - position);
        }

        return ret;
    }

    void find(bvh_node *x, Ray *ray, real h) {
        if (h >= hit)
            return;
        ++count;
        if (x->size <= MIN_SIZE) {
            for (int i = 0; i < x->size; ++i) {
                real t = x->obj[i]->intersect(ray);
                if (t < hit) {
                    hit = t;
                    target = x->obj[i];
                }
            }
            return;
        }

        real L = x->lc->intersect(ray), R = x->rc->intersect(ray);
        if (L < R) {
            find(x->lc, ray, L);
            find(x->rc, ray, R);
        } else {
            find(x->rc, ray, R);
            find(x->lc, ray, L);
        }
    }

    void profile() {
        printf("[bvh] count = %ld, intersect = %ld, rate = %.3f\n", count, count_intersect, 1. * count / (count_intersect + 1e-6));
        count = count_intersect = 0;
    }

    pair<real, Object *> slow_intersect(Ray ray) {
        hit = infinity;
        target = nullptr;
        for (int i = 0; i < n; ++i) {
            // if (objects[i]->id == 7779 || objects[i]->id == 1888) {
            //     printf("y\n");
            // }
            real t = objects[i]->intersect(&ray);
            if (t < hit) {
                hit = t;
                target = objects[i];
            }
        }

        return {hit, target};
    }

    pair<real, Object *> intersect(Ray ray) {
        count_intersect += 1;
        hit = infinity;
        target = nullptr;

        find(root, &ray, root->intersect(&ray));

        return {hit, target};
    }
};