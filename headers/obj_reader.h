#include "vec3.h"
#include <cstdio>
#include <cstring>
#include "bvh.h"
#include "brdf.h"
using namespace std;

// see https://www.wikiwand.com/en/Wavefront_.obj_file

void obj_reader(const char *filename, Env *env, bool faceNormals) {
    vector<Vec3> vertices, normals;
    char type[5];

    FILE *f = fopen(filename, "r");

    BRDF *brdf = new Refraction(Vec3(1, 1, 1) * 0.99, 1.8);
    // BRDF *brdf = new Diffusion(Vec3(0.75, 0.25, 0.25));

    Vec3 L = Vec3(1, 1, 1) * +infinity;
    Vec3 R = Vec3(1, 1, 1) * -infinity;

    int num_triangles = 0;
    while (~fscanf(f, "%s", type)) {
        if (strcmp(type, "v") == 0) {
            double x, y, z;
            fscanf(f, "%lf%lf%lf", &x, &y, &z);
            y *= 1.5; // transform;
            vertices.push_back(Vec3(x, y, z));
        } else if (strcmp(type, "f") == 0) {
            int a, b, c, na, nb, nc;
            // ignore normal information
            fscanf(f, "%d/%*d/%d %d/%*d/%d %d/%*d/%d", &a, &na, &b, &nb, &c, &nc);
            --a, --b, --c, --na, --nb, --nc;
            if (faceNormals) {
                env->bvh.add(new TriangleNormal({vertices[a], normals[na]}, {vertices[b], normals[nb]}, {vertices[c], normals[nc]}, brdf));
                // L = min(L, vertices[a]);
                // L = min(L, vertices[b]);
                // L = min(L, vertices[c]);
                // R = max(R, vertices[a]);
                // R = max(R, vertices[b]);
                // R = max(R, vertices[c]);
            } else {
                env->bvh.add(new Triangle(vertices[a], vertices[b], vertices[c], brdf));
            }
            num_triangles += 1;
        } else if (strcmp(type, "vt") == 0) { // texture
            double x, y;
            fscanf(f, "%lf%lf", &x, &y);
            // what is that?
        } else if (strcmp(type, "vn") == 0) {
            double x, y, z;
            fscanf(f, "%lf%lf%lf", &x, &y, &z);
            normals.push_back(Vec3(x, y, z));
        } else {
            printf("unknown type: %s\n", type);
        }
    }
    // printf("(%.3f, %.3f, %.3f) (%.3f, %.3f, %.3f)\n", L.x, L.y, L.z, R.x, R.y, R.z);

    fclose(f);
    printf("read file %s: # vertices = %d, # triangles = %d\n", filename, (int) vertices.size(), num_triangles);
}