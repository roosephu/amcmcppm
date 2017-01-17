#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>
#include <random>
#include <list>
// #include <boost/unordered_map.hpp>
#include <cassert>
#define STB_IMAGE_WRITE_IMPLEMENTATION

typedef double real;

int debug = 0;

#include "stb_image_write.h"
// #include <opencv2/opencv.hpp>
#include "headers/vec3.h"
#include "headers/hash.h"
#include "headers/mcmc.h"
#include "headers/bvh.h"
#include "headers/light.h"

using namespace std;
// using namespace cv;

#define HASH_LEVEL 5

// const real threshold = 1e-6;
const int scale = 1, width = 1024, height = 1024;
// const int width = 16, height = 12;
const int num_stages = 10000, num_photons = 1000000;
const real alpha = 0.3; // ppm
const real initial_cube_size = 0.05;
// const int super_sampling = 1;

enum LightSource {
    EYE,
    LIGHT,
};

struct HitPoint {
    Vec3 loc, normal;
    Color weight, flux;
    int count, clock;
    real rad_sqr;

    HitPoint() {}

    HitPoint(Vec3 loc_, Vec3 normal_, Vec3 weight_)
        : loc(loc_), normal(normal_), weight(weight_), count(0), clock(-1)
    {
    }
};

// struct AreaLight : Light {
//     AreaLight(Vec3 loc_, Vec3 dir_, real dist_, Vec3 flux_)
//     {
//         : loc(loc_), dir(dir_), dist(dist_)
//         flux = flux_;
//         if (abs(dir.x) < 0.5)
//             cx = Vec3(1, 0, 0);
//         else if (abs(dir.y) < 0.5)
//             cx = Vec3(0, 1, 0);
//         else
//             cx = Vec3(0, 0, 1);
//         cy = dir % cx;
//     }

//     Ray gen_ray(PathSpace &path) {
//         real u1 = 2 * path.next() - 1, u2 = 2 * path.next() - 1;
//         return Ray(loc + cx * u1 + cy * u2, gen_sphere(path));
//     }
// };

struct Env {
    int clock = 0;
    BVH bvh;
    Light *light;
    // Buffer<HitPoint, int(1.2e6)> hit_points;
    HitPoint hit_points[width * height];
    Buffer<ListNode<HitPoint *>, int(1e7)> hash_buf;
    Hash<HitPoint *> H[HASH_LEVEL];
    real cube_sizes[HASH_LEVEL];
    pair<int, int> pixel;
    bool visited;
    long hash_test[HASH_LEVEL] = {0}, hash_hit[HASH_LEVEL] = {0}, hash_match[HASH_LEVEL] = {0};
    double hash_rate[HASH_LEVEL] = {1, 1, 1, 1, 1};

    unsigned long hash(int level, Vec3 v) {
        return hash(floor(v.x / cube_sizes[level]), floor(v.y / cube_sizes[level]), floor(v.z / cube_sizes[level]));
    }

    unsigned long hash(unsigned long x, unsigned long y, unsigned long z) {
        // return (x << 32) + (y << 16) + z;
        return (x * 1000000007 + y) * 1000000007 + z;
    }

    Env() {
        for (int i = 0; i < HASH_LEVEL; ++i) {
            hash_rate[i] = 1.;
            cube_sizes[i] = initial_cube_size;
        }
    }

    void hash_profile() {
        // if (hash_test == 0)
        //     return;
        printf("[hash] \n");
        long total_test = 0, total_hit = 0;
        for (int i = 0; i < HASH_LEVEL; ++i) {
            if (hash_test[i] != 0) {
                double rate = (double) (hash_hit[i] + 1e-8) / (hash_test[i] + 1e-8);
                hash_rate[i] = hash_rate[i] * rate / hash_test[i];
                printf("[hash] {test = %8ld, match = %8ld, hit = %8ld, rate = %.3f, eff = %.3f}\n", hash_test[i], hash_match[i], hash_hit[i], rate, 1. * hash_match[i] / hash_test[i]);
                total_test += hash_test[i], total_hit += hash_hit[i];
                hash_test[i] = hash_hit[i] = hash_match[i] = 0;
            }
        }
        if (total_test != 0) {
            printf("[hash] total: test = %ld, hit = %ld, rate = %.3f\n", total_test, total_hit, 1. * total_hit / total_test);
        }
    }

    void update_hit_points(Vec3 intersection, Vec3 normal, Vec3 flux) {
        for (int h = 0; h < HASH_LEVEL; ++h) {
            auto hash_val = hash(h, intersection);
            for (auto node = H[h].find(hash_val); node; node = node->next) {
                ++hash_test[h];
                if (node->key != hash_val) {
                    // if (rand() <= 100000) {
                    //     printf("%lu %lu\n", hash_val, node->key);
                    // }
                    continue;
                }
                ++hash_match[h];
                HitPoint *hit_point = node->item;
                if (dot(hit_point->normal, normal) > 1e-3) {
                    Vec3 v = intersection - hit_point->loc;
                    // hit_point->flux = Vec3(0.75, 0.25, 0.25);
                    real distance = dot(v, v);
                    assert(distance <= cube_sizes[h] * 1.8);
                    if (distance <= hit_point->rad_sqr) {
                        // if (hit_point == hit_points + 24) {
                        //     printf("intersect with (0, 24): (%f, %f, %f)\n", intersection.x, intersection.y, intersection.z);
                        // }
                        ++hash_hit[h];
                        visited = true;
                        real m = 1;
                        real g = (hit_point->count + m) / (hit_point->count + m / alpha);
                        hit_point->rad_sqr *= g;
                        hit_point->count += m;
                        hit_point->flux = (hit_point->flux + flux * hit_point->weight / pi) * g;
                        assert(flux.x >= 0 && flux.y >= 0 && flux.z >= 0);
                        assert(hit_point->weight.x >= 0 && hit_point->weight.y >= 0 && hit_point->weight.z >= 0);
                    }
                }
            }
        }
    }

    // vector<real> lengths;
    void trace(Ray ray, LightSource mode, PathSpace &path, Vec3 flux) {
        bool inside = false;
        real total = 0;

        for (int depth = 1; depth <= 10; ++depth) {
            auto _i = bvh.intersect(ray);
            // if (mode == LIGHT) {
            //     bvh.profile();
            //     // printf("intersection: %f %f %f\n", intersection.x, intersection.y, intersection.z);
            // }
            auto t = _i.first;
            auto obj = _i.second;
            if (obj == nullptr)
                break;

            Vec3 intersection = ray.loc + ray.dir * t;
            Vec3 normal = obj->get_normal(intersection).normalize();
            real in_dot_normal = dot(normal, ray.dir);
            bool into = true;

            if (in_dot_normal > 0) {
                in_dot_normal *= -1;
                normal = normal * -1;
                into = false;
            }

            auto brdf = obj->get_brdf(path);
            auto ray_out = brdf->sample(path, ray.dir, normal, into);

            flux = flux * ray_out.second;
            if (inside) {
                assert(t > 0);
                real absorption_rate = exp(-t * 0.61);
                // total += t;
                flux = flux * (Vec3(0.63, 0.45, 0.09) * (1 - absorption_rate) + absorption_rate);
            }
            if (dot(ray_out.first, normal) < 0) // get through
                inside = !inside;
            // if (!brdf->is_delta) {
            //     ray_out.first = (Vec3(-1.050568, 2.000000, 0.741164) - intersection).normalize();
            // }

            if (brdf->is_delta) {
                ray = Ray(intersection, ray_out.first);
            } else if (mode == EYE) {
                int index = pixel.first * height + pixel.second;
                auto hit_point = hit_points + index;
                hit_point->loc = intersection;
                hit_point->normal = normal;
                hit_point->weight = flux;
                // hit_point->flux = Vec3();
                hit_point->clock = clock;
                assert(flux.x >= 0 && flux.y >= 0 && flux.z >= 0);
                break;
            } else {
                update_hit_points(intersection, normal, flux);

                // Russian Roulette
                if (depth >= 100) {
                    real p = ray_out.second.max();
                    if (path.next() > p)
                        break;
                    flux = flux / p;
                }
                ray = Ray(intersection, ray_out.first);
            }
        }
        // lengths.push_back(total);
        // if (total != 0 && mode != EYE)
        //     printf("%f\n", total);
    }

    void initialize_hit_points() {
        // Vec3 L = Vec3(1, 1, 1) * +infinity;
        // Vec3 R = Vec3(1, 1, 1) * -infinity;

        // for (int i = 0; i < width * height; ++i)
        //     auto hit_point = hit_points + i;
        //     L = min(L, hit_point->loc);
        //     R = max(R, hit_point->loc);
        // }

        // real rad = 0.00005; //(xmax + ymax + zmax - xmin - ymin - zmin) / 3 / (width + height) * 4;
        // cube_size = rad * 2;
        // real rad = initial_cube_size / 2;

        // for (auto &hit_point : hit_points) {
        //     hit_point.rad_sqr = rad * rad;
        // }
        // printf("%d hit points initialized.\n", (int)hit_points.size());
        // assert(hit_points.size() != 0);
    }


    void add_hit_point_to_hash(HitPoint *hit_point, int level) {
        real rad = sqrt(hit_point->rad_sqr);
        double cube_size = cube_sizes[level];
        // if (cube_size < rad) {
        //     printf("xxx count = %d, level = %d\n", hit_point->count, level);
        // }

        int xl = floor((hit_point->loc.x - rad) / cube_size), xr = floor((hit_point->loc.x + rad) / cube_size);
        int yl = floor((hit_point->loc.y - rad) / cube_size), yr = floor((hit_point->loc.y + rad) / cube_size);
        int zl = floor((hit_point->loc.z - rad) / cube_size), zr = floor((hit_point->loc.z + rad) / cube_size);
        // if (hit_point == hit_points + 24) {
        //     printf("hash 24: %d %d %d, rad = %f, cube size = %f\n", xl, yl, zl, rad, cube_size);
        // }

        for (int x = xl; x <= xr; ++x) {
            for (int y = yl; y <= yr; ++y) {
                for (int z = zl; z <= zr; ++z) {
                    auto key = hash(x, y, z);
                    ListNode<HitPoint *> *node = hash_buf.get();
                    node->key = key;
                    node->item = hit_point;
                    node->next = nullptr;
                    H[level].insert(key, node);
                }
            }
        }
    }

    void build_hash(double theta) {
        hash_profile();
        // cube_size = 0;
        // for (auto hit_point : hit_points) {
        //     cube_size += pow(hit_point->rad_sqr, 2);
        // }
        // cube_size /= hit_points.size();
        // cube_size = pow(cube_size, 1. / 4);
        hash_buf.reset();
        for (int i = 0; i < HASH_LEVEL; ++i)
            H[i].clear();

        vector<pair<int, HitPoint *>> rads;
        rads.reserve(width * height);
        for (int i = 0; i < width * height; ++i) {
            if (hit_points[i].clock == clock) {
                rads.push_back({hit_points[i].count, hit_points + i});
            }
        }
        sort(rads.begin(), rads.end());
        // reverse(rads.begin(), rads.end());

        // assert(HASH_LEVEL > 1);

        // int start = 0;
        // for (; start < rads.size() && rads[start].first == 0; ++start)
        //     add_hit_point_to_hash(rads[start].second, 0);

        // cube_sizes[0] = initial_cube_size * initial_cube_size;
        // for (int i = 0; i < HASH_LEVEL - 1; ++i)
        //     cube_sizes[i + 1] = cube_sizes[i] * (i + 1) / (i + 1 / alpha);
        // for (int i = 0; i < HASH_LEVEL; ++i)
        //     cube_sizes[i] = sqrt(cube_sizes[i]);

        // int step = rads.size() / HASH_LEVEL + 1;
        // for (int i = 0; i < HASH_LEVEL; ++i)
        //     cube_sizes[i] = sqrt(rads[i * step].second->rad_sqr) * 2;

        int count[HASH_LEVEL] = {0}, num_hit_points = rads.size();

        int starts[HASH_LEVEL] = {0};
        double sum = 0;
        for (int i = 0; i < HASH_LEVEL; ++i)
            sum += hash_rate[i];
        for (int i = 0, last = 0; i < HASH_LEVEL; ++i) {
            starts[i] = last;
            last += hash_rate[i] / sum * num_hit_points;
            hash_rate[i] /= sum;
        }

        int level = -1;
        for (int i = 0; i < num_hit_points; ++i) {
            // int level = min(rads[i].first, HASH_LEVEL - 1);
            // int level = pow(1. * i / num_hit_points, 1 / theta) * HASH_LEVEL;
            real rad = sqrt(rads[i].second->rad_sqr);

            // the following code is wrong: because we don't add points with `clock` != clock
            // if (level < HASH_LEVEL - 1 && i >= starts[level + 1] && cube_sizes[level + 1] >= rad * 2) {
            if (i % (num_hit_points / HASH_LEVEL + 1) == 0) {
                ++level;
                cube_sizes[level] = min(cube_sizes[level], rad * 2);
            }
            add_hit_point_to_hash(rads[i].second, level);
            ++count[level];
            // if (i % 100000 == 0) {
            //     printf("%d: %d\n", i, hash_buf.size());
            // }
        }
        // for (int i = 1; i < HASH_LEVEL; ++i) {
        //     int L = (i - 1) * (rads.size() - start) / (HASH_LEVEL - 1) + start, R = i * (rads.size() - start) / (HASH_LEVEL - 1) + start;
        //     cube_sizes[i] = 2 * sqrt(rads[L].second->rad_sqr);

        //     for (int idx = L; idx < R; ++idx) {
        //         auto hit_point = rads[idx].second;
        //         add_hit_point_to_hash(hit_point, i);
        //     }
        // }

        printf("[hash] ");
        for (int i = 0; i < HASH_LEVEL; ++i) {
            printf("(%.8f: %d/%d), ", cube_sizes[i], count[i], H[i].validEntries());
        }
        printf("\n");
    }

    bool visible(PathSpace &path, real scale) {
        visited = false;
        path.reset();
        trace(light->gen_ray(path), LIGHT, path, light->flux * scale);
        return visited;
    }

    void reset() {
        clock += 1;
    }
};
#include "headers/ply_reader.h"
#include "headers/obj_reader.h"

int gamma_correction(real x) {
    return int(pow(max(min(x, (real) 1.), (real) 0.), 1 / 2.2) * 255 + 0.5);
}

Vec3 img[width][height], aux[width][height];
char image_data[width * height * 3];

void draw_pixel(int i, int j, Vec3 c) {
    image_data[(j * width + i) * 3 + 0] = gamma_correction(c.x);
    image_data[(j * width + i) * 3 + 1] = gamma_correction(c.y);
    image_data[(j * width + i) * 3 + 2] = gamma_correction(c.z);
}
Env env;

void add_rectangle(Vec3 a, Vec3 b, Vec3 c, Vec3 d, BRDF *brdf) {
    env.bvh.add(new Triangle(a, b, c, brdf));
    env.bvh.add(new Triangle(b, c, d, brdf));
}

void write_image(Env *env, double scale) {
    double sum = 0;
    for (int i = 0; i < width * height; ++i) {
        int x = i / height, y = i % height;
        // auto x = hit_point.pixel.first, y = hit_point.pixel.second;

        img[x][y] = env->hit_points[i].flux / env->hit_points[i].rad_sqr;
        sum += img[x][y].sum();
        // img[x][y] = hit_point.normal; scale = 1;
        // printf("%f %f %f\n", hit_point.normal.x, hit_point.normal.y, hit_point.normal.z);
    }
    scale = 3 * 0.18 * width * height / sum;
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            draw_pixel(i, j, img[i][j] * scale);
        }
    }

    stbi_write_png("outputs/output.png", width, height, 3, image_data, width * 3);

}

int main() {

    // ply_reader("dragon_recon/dragon_vrip.ply", &env);
    // ply_reader("tmp.ply", &env);
    obj_reader("water-caustic/models/Mesh000.obj", &env, false);
    obj_reader("water-caustic/models/Mesh001.obj", &env, true);

    // for (int i = 0; i < HASH_LEVEL; ++i) {
    //     printf("%f\n", env.cube_sizes[i]);
    // }

    real ymax = 2;
    add_rectangle( // floor
        Vec3(1, 0, 1),
        Vec3(1, 0, -1),
        Vec3(-1, 0, 1),
        Vec3(-1, 0, -1),
        new Diffusion(Vec3(0.725, 0.71, 0.68))
    );
    add_rectangle( // ceiling
        Vec3(-1, 2., -1),
        Vec3(1, 2., -1),
        Vec3(-1, 2., 1),
        Vec3(1, 2., 1),
        new Diffusion(Vec3(0.725, 0.71, 0.68))
    );
    add_rectangle( // back wall
        Vec3(1, ymax, -1),
        Vec3(1, 0., -1),
        Vec3(-1, ymax, -1),
        Vec3(-1, 0., -1),
        new Diffusion(Vec3(0.725, 0.71, 0.68))
    );
    add_rectangle( // right wall
        Vec3(1, ymax, 1),
        Vec3(1, 0., 1),
        Vec3(1, ymax, -1),
        Vec3(1, 0., -1),
        new Diffusion(Vec3(0.14, 0.45, 0.091))
    );
    add_rectangle( // left wall
        Vec3(-1, ymax, -1),
        Vec3(-1, 0., -1),
        Vec3(-1, ymax, 1),
        Vec3(-1, 0., 1),
        new Diffusion(Vec3(0.63, 0.065, 0.05))
    );
    // return 0;

    // BRDF *mixed = new MixedBRDF({
    //     // {new Reflexion(Vec3(1,1,1)*.996), 0.4},
    //     {new Refraction(Vec3(1,1,1)*.996, 1.5), 0.4},
    //     {new Diffusion(Vec3(0.275, 0.612, 0.949)), 0.3},
    //     // {new Blur(new Refraction(Vec3(1,1,1)*.996, 1.5)), 1},
    //     // {new Refraction(Vec3(1,1,1)*.996, 1.5), 1},
    // });

    // env.bvh.add(new Sphere(Vec3(  1e5 + 1, 40.8, 81.6), 1e5, new Diffusion(Vec3(.75,.25,.25))));
    // env.bvh.add(new Sphere(Vec3(-1e5 + 99, 40.8, 81.6), 1e5, new Diffusion(Vec3(.25,.25,.75))));
    // env.bvh.add(new Sphere(Vec3(       50, 40.8,  1e5), 1e5, new Diffusion(Vec3(.75,.75,.75))));
    // // env.bvh.add(new Sphere(Vec3(       50, 40.8, -1e5+170), 1e5, new Diffusion(Vec3())));
    // env.bvh.add(new Sphere(Vec3(       50,  1e5, 81.6), 1e5, new Diffusion(Vec3(.75,.75,.75))));
    // env.bvh.add(new Sphere(Vec3(       50, -1e5+81.6,81.6), 1e5, mixed));
    // env.bvh.add(new Sphere(Vec3(27,16.5,47), 16.5, new Reflexion(Vec3(1,1,1)*.999)));
    // env.bvh.add(new Sphere(Vec3(35,10,100), 10, new Refraction(Vec3(1,1,1)*.999, 1.5)));
    // env.bvh.add(new Sphere(Vec3(50,8.5,60), 8.5, new Diffusion(Vec3(1,1,1)*.999)));
    // env.light = new PointLight(Vec3(50,3000,85), Vec3(1,1,1)*6e8, 1.5);
    // env.bvh.add(new Sphere(Vec3(0.0, 0.1, 0.0), 0.05, new Refraction(Vec3(0.99, 0.99, 0.99), 1.5)));
    // env.bvh.add(new Triangle(Vec3(100, 0.05, -173), Vec3(0, 0.05, 200), Vec3(-100, 0.05, -173), new Diffusion(Vec3(0.75, 0.25, 0.25))));

    //real R=60;
    // real R = 120;     // radius
    // real T = 30 * pi / 180.;
    // real D = R / cos(T);     //distance
    // // real D=60;     //distance
    // // real R=D*sqrt(2);
    // real Z = 62;
    // Vec3 C = Vec3(0.275, 0.612, 0.949);

    // env.bvh.add(new Sphere(Vec3(50, 28, Z) + Vec3(cos(T), sin(T), 0) * D, R, mixed)); // red
    // env.bvh.add(new Sphere(Vec3(50,28,Z)+Vec3(-cos(T),sin(T),0)*D, R, mixed));
    // env.bvh.add(new Sphere(Vec3(50,28,Z)+Vec3(0,-1,0)*D, R, mixed));
    // env.bvh.add(new Sphere(Vec3(50,28,Z)+Vec3(0,0,-1)*R*2*sqrt(2./3.), R, mixed));
    // //  Sphere(1e5, Vec(50,28,Z)+Vec(0,0,1e5+170),   Vec(1,1,1)*0,Vec(1,1,1)*.996, SPEC), //front
    // //  Sphere(2*R*2*sqrt(2./3.)-R*2*sqrt(2./3.)/3., Vec(50,28,Z)+Vec(0,0,-R*2*sqrt(2./3.)/3.),   Vec(1,1,1)*0,Vec(1,1,1)*.3333, SPEC), //front
    // env.bvh.add(new Sphere(Vec3(50,28,Z)+Vec3(0,0,-R*2*sqrt(2./3.)/3.), 2*2*R*2*sqrt(2./3.)-R*2*sqrt(2./3.)/3, mixed));
    // env.bvh.add(new Sphere(Vec3(50, 28, Z), R * 8, new Diffusion(Vec3(0.75, 0.25, 0.25))));
    env.bvh.build();

    // env.spheres.push_back(new Sphere(Vec3(50,681.6-.27,81.6), 600, DIFFUSION, Vec3(12,12,12), Vec3()));

    // env.spheres.push_back(new Sphere(Vec3(50, 50, -150), 50, REFRACTION, Vec3(), Vec3(0.75, 0.25, 0.25)));
    // env.spheres.push_back(new Sphere(Vec3(50, 150, -150), 50, DIFFUSION, Vec3(8, 8, 8), Vec3()));
    // env.spheres.push_back(new Sphere(Vec3(50, -1e5, -150), 1e5, DIFFUSION, Vec3(), Vec3(0.25, 0.25, 0.75)));

    // Ray cam(Vec3(50, 48, 295.6), Vec3(0, -0.042612, -1));
    // Ray cam(Vec3(0.4, 0.25, 0.5), Vec3(-1.0, -0.25, -1));
    // Ray cam(Vec3(0.5, 0.5, 0.5), Vec3(-1.20, -0.714, -1));
    // cam.loc = cam.loc;
    // env.light = new PointLight(Vec3(50,28,Z), Vec3(1, 1, 1) * 1e7, 0);
    // env.light = new PointLight(Vec3(1.5, 0.5, 0.8), Vec3(1, 1, 1) * 3e3, 0);
    // env.light = new PlaneLight(cam.loc, cam.dir, 2, Vec3(1, 1, 1) * 1e7);
    // env.light = Light(Vec3(1.5e1, 2.5e0, 1e1), Vec3(1, 1, 1) * 1e6, 1e1);
    // Vec3 cx(width * 0.5135 / height), cy = (cx % cam.dir).normalize() * 0.5135;
    // m = Mat(height, width, CV_8UC3, 0.);

    Ray cam(Vec3(0, 1,  6.838), Vec3(0, 0, -5.838));
    Vec3 cx = Vec3(1, 0, 0) / 2.9;
    Vec3 cy = Vec3(0, 1, 0) / 2.9;
    // Ray cam(Vec3(0, 1.9999, 0), Vec3(0, -1, 0));
    // Vec3 cx = Vec3(1, 0, 0);
    // Vec3 cy = Vec3(0, 0, 1);

    // const real eps = 1e-6;
    // Ray cam(Vec3(0, 0, 0), Vec3(eps, eps, -1));
    // Vec3 cx(1, eps, eps);
    // Vec3 cy(eps, 1, eps);
    // env.light = new PointLight(Vec3(-0.005, 1.98, -0.03), Vec3(541127, 381972, 127324) * 0.0003, 0.001);
    // env.light = new PointLight(Vec3(-0.005, 1.98, -0.03), Vec3(541127, 381972, 127324) * 1e-4, 0.00252);
    // env.light = new PlaneLight(Vec3(-0.005, 1.98, -0.03), Vec3(0.005, -1.48, 0.03), 0.00252, Vec3(541127, 381972, 127324) * 6e-4);
    env.light = new SemisphereLight(Vec3(0, 1.99998, 0), Vec3(0, -1, 0), Vec3(541127, 381972, 127324), 0.00252);

    bool initialized = false;
    PathSpace current_path;
    real mutation_size = 1;
    long accepted = 1, mutated = 0, uniform_count = 1;
    long rebuild_hash_level = 10000, current_photons = 0;

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            int index = i * height + j;
            env.hit_points[index] = HitPoint(Vec3(), Vec3(), Vec3());
            env.hit_points[index].rad_sqr = initial_cube_size * initial_cube_size / 4;
        }
    }
    PathSpace random_path(false);

    // const double constants[] = {1, 1.1, 1.7, 1.9, 2.4, 2.5, 2.5, 2.5, 2.5, 3, 3, 3, 3, 3, 3, 3};
    for (int stage = 1; stage <= num_stages; ++stage) {
        printf("[clock] iteration #%d: %.3f\n", stage, 1. * clock() / CLOCKS_PER_SEC);
        // int i1 = 255 * 512 + 206;
        // printf("(255, 206): (%.3f, %.3f, %.3f), rad = %.3f\n", env.hit_points[i1].loc.x, env.hit_points[i1].loc.y, env.hit_points[i1].loc.z, env.hit_points[i1].rad_sqr);

        env.reset();
        // build photon mapping
        for (int i = 0; i < width; ++i) {
            fprintf(stderr, "\r[eye phase]: %2.0f%%", 100. * i / width);
            for (int j = 0; j < height; ++j) {
                env.pixel = {i, j};
                // if (i == 0 && j == 24) {
                //     printf("en?\n");
                // }

                Vec3 d = cx * ((i + rand_gen()) / width - 0.5) + cy * (-(j + rand_gen()) / height + 0.5) + cam.dir;
                env.trace(Ray(cam.loc, d), EYE, random_path, Vec3(1, 1, 1));

            }
            // env.bvh.profile();
        }
        // printf("(256, 206): (%.3f, %.3f, %.3f), rad = %.3f\n", env.hit_points[i2].loc.x, env.hit_points[i2].loc.y, env.hit_points[i2].loc.z, env.hit_points[i2].rad_sqr);
        fprintf(stderr, "\r");
        env.bvh.profile();
        env.initialize_hit_points();
        // env.hash_profile();
        env.build_hash(0);

        // printf("(0, 24): %f %f %f\n", env.hit_points[24].loc.x, env.hit_points[24].loc.y, env.hit_points[24].loc.z);

        int init_count = 0;
        while (!initialized) {
            current_path = PathSpace();
            initialized = env.visible(current_path, 0);
            init_count += 1;
            assert(init_count < 1000000);
        }
        if (init_count > 0) {
            printf("initialize MCMC: %d\n", init_count);
        }

        // Vec3 T(-0.1021928149, 0.1420205835, -0.0016782738);

        // sample from light
        for (int i = 1, percentile = max(1, num_photons / 100); i <= num_photons; ++i) {
            // printf("%d\n", i);
            // debug = i;
            // env.hash_profile();
            // env.bvh.profile();
            ++current_photons;
            if (i % percentile == 0)
                fprintf(stderr, "\riteration #%d: progress: %2.0f%%", stage, 100. * i / num_photons);
            // env.trace(Ray(env.light.loc, T - env.light.loc), 0, LIGHT, current_path, Vec3(1, 1, 1) * 1e7);

            auto uniform_path = PathSpace();
            // env.visible(uniform_path);
            if (env.visible(uniform_path, 1)) {
                uniform_count += 1;
                current_path = uniform_path;
            } else {
                auto candidate_path = current_path.mutate(mutation_size);
                mutated += 1;
                // real scale = (real)uniform_count / (uniform_count + mutated);
                real scale = 1;
                if (env.visible(candidate_path, scale)) {
                    current_path = candidate_path;
                    accepted += 1;
                } else {
                    env.visible(current_path, scale);
                }
                real ratio = 1. * accepted / mutated;
                mutation_size += (ratio - 0.234) / mutated;
                mutation_size = max(mutation_size, (real) 1e-3);
            }

            if (current_photons == rebuild_hash_level) {
                env.build_hash(0);
                rebuild_hash_level = rebuild_hash_level * 2 + 1;
                // env.hash_profile();
            }
            // env.visible(current_path, 1. * uniform_count / (uniform_count + mutated));
        }
        fprintf(stderr, "\r");
        // real scale = 1. * uniform_count / current_photons;
        real scale = 1.;
        write_image(&env, scale / current_photons);
        // for (int i = 0; i < width * height; ++i) {
        //     if (env.hit_points[i].count == 0 && env.hit_points[i].clock == env.clock) {
        //         printf("%d %d\n", i / height, i % height);
        //     }
        // }

        // sort(env.lengths.begin(), env.lengths.end());
        // printf("%f %f\n", env.lengths[env.lengths.size() / 3], env.lengths[2 * env.lengths.size() / 3]);

        env.hash_profile();
        printf("iteration #%d: uniform count = %.f, accepted = %.f, mutated = %.f, mutation_size = %.4f\n",
            stage, 1. * uniform_count / stage, 1. * accepted / stage, 1. * mutated / stage, mutation_size);
    }

    return 0;
}
