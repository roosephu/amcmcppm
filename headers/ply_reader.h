#include <cstdio>
#include <cassert>
#include <cstdlib>
#include "vec3.h"
#include "bvh.h"
#include <map>
using namespace std;

void ply_reader(const char *filename, Env *env) {
    FILE *f = fopen(filename, "r");
    int n, m;

    fscanf(f, "%*s");
    fscanf(f, "%*s %*s %*s %*s");
    fscanf(f, "%*s %*s %*s");
    fscanf(f, "%*s %*s %d", &n);
    fscanf(f, "%*s %*s %*s");
    fscanf(f, "%*s %*s %*s");
    fscanf(f, "%*s %*s %*s");
    fscanf(f, "%*s %*s %d", &m);
    fscanf(f, "%*s %*s %*s %*s %*s");
    fscanf(f, "%*s");

    vector<Vec3> vertices;
    for (int i = 0; i < n; ++i) {
        double x, y, z;
        fscanf(f, "%lf%lf%lf", &x, &y, &z);
        vertices.push_back(Vec3(x, y, z) * 200 + Vec3(75, -10, 80));
    }
    // map<pair<int, int>, vector<Triangle *>> edges;

    // BRDF *brdf = new Diffusion(Vec3(0.25, 0.25, 0.75));
    BRDF *brdf = new Refraction(Vec3(0.996, 0.996, 0.996), 1.5);
    // BRDF *brdf2 = new Diffusion(Vec3(0.25, 0.25, 0.75));

    // int count = 0;
    for (int i = 0; i < m; ++i) {
        int t, a, b, c;
        fscanf(f, "%d%d%d%d", &t, &a, &b, &c);
        assert(t == 3);

        // BRDF *brdf = new Diffusion(Vec3(1. * rand() / RAND_MAX, 1. * rand() / RAND_MAX, 1. * rand() / RAND_MAX));
        Triangle *tri = new Triangle(vertices[a], vertices[b], vertices[c], brdf);
        tri->id = i;
        // if (max(max(vertices[a], vertices[b]), vertices[c]).x <= -0.1) {
        env->bvh.add(tri);
        //     ++count;
        // }

        // edges[{min(a, b), max(a, b)}].push_back(tri);
        // edges[{min(b, c), max(b, c)}].push_back(tri);
        // edges[{min(c, a), max(c, a)}].push_back(tri);

        // edges[a].push_back(tri);
        // edges[b].push_back(tri);
        // edges[c].push_back(tri);
    }
    // for (auto it : edges) {
    //     if (it.second.size() > 2) {
    //         // printf("xxx\n");
    //     }
    // }

    int t;
    assert(fscanf(f, "%d", &t) == -1);
    printf("read ply file: %s with %d vertices and %d triangles\n", filename, n, m);
}
