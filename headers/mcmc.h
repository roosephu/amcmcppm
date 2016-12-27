#pragma once

#include <vector>
#include <cassert>
using namespace std;

mt19937 engine(100);
uniform_real_distribution<real> uniform;
normal_distribution<real> gaussian;

real rand_gen() {
    return uniform(engine);
}

real rand_gen_normal() {
    return gaussian(engine);
}

struct PathSpace {
    vector<real> values;
    int index;
    bool save;

    PathSpace(bool save_=true) : index(0), save(save_) {}

    real next() {
        real ret;
        if (save) {
            while ((int) values.size() <= index) {
                values.push_back(rand_gen());
            }
            ret = values[index];
            index ++;
        } else {
            ret = rand_gen();
        }
        return ret;
    }

    void reset() {
        index = 0;
    }

    PathSpace mutate(real mutation_size) {
        assert(mutation_size > 0);
        PathSpace ret;
        ret.values = values;
        for (auto &value : ret.values) {
            real delta = pow(rand_gen(), 1 / mutation_size + 1);
            if (rand_gen() < 0.5) {
                value += delta;
                if (value > 1)
                    value -= 1;
            } else {
                value -= delta;
                if (value < 0)
                    value += 1;
            }
            // ret.values.push_back(value);
        }
        return ret;
    }

};

Vec3 gen_sphere(PathSpace &path) {
    real z = path.next() * 2 - 1, r = catheti(1, z);
    real theta = 2 * pi * path.next();
    return Vec3(r * cos(theta), r * sin(theta), z);
}

Vec3 _gen_sphere() {
    real z = rand_gen() * 2 - 1, r = catheti(1, z);
    real theta = 2 * pi * rand_gen();
    return Vec3(r * cos(theta), r * sin(theta), z);
}
