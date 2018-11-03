#ifndef SPHDATASTRUCTURES_H
#define SPHDATASTRUCTURES_H

#include "CGALTYPEDEFS.h"

struct Kernel {

    std::vector<float> W;

    std::vector<Vector> dWdx;
};

struct ParticleNeighbours {

    // TODO This should be sorted hence replace by number of entries per
    // particle
    std::vector<size_t> origId;

    std::vector<size_t> neighId;

    std::vector<Vector> distances;

    std::vector<Vector> normalised_distances;

    std::vector<float> squared_length;

    std::vector<Vector> subtract (const std::vector<Vector> vals) const {

        const size_t res_size = neighId.size();
        std::vector<Vector> res (res_size);

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t nid = neighId[ctr];
            const size_t oid =  origId[ctr];

            res[ctr] = {vals[oid] - vals[nid]};
        }
        return res;
    }

    std::vector<float> subtract (const std::vector<float> vals) const {

        const size_t res_size = neighId.size();
        std::vector<float> res (res_size);

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t nid = neighId[ctr];
            const size_t oid =  origId[ctr];

            res[ctr] = {vals[oid] - vals[nid]};
        }
        return res;
    }

    void sum (const std::vector<float>& vals, std::vector<float>& ret) const {

        const size_t res_size = neighId.size();

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t oid =  origId[ctr];
            const size_t nid = neighId[ctr];

            // std::cout
            //     << "[DEBUG] nid " << nid
            //     << " oid " << oid
            //     << " vals[oid] " << vals[oid]
            //     << " ret[oid] " << ret[oid]
            //     << std::endl;

            ret[oid] += vals[nid];
        }
    }

    void sum (const std::vector<Vector>& vals, std::vector<Vector>& ret) const {

        const size_t res_size = neighId.size();

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t oid =  origId[ctr];
            const size_t nid = neighId[ctr];

            // std::cout
            //     << "[DEBUG] nid " << nid
            //     << " oid " << oid
            //     << " vals[oid] " << vals[oid]
            //     << " ret[oid] " << ret[oid]
            //     << std::endl;

            ret[oid] += vals[nid];
        }
    }
};

#endif
