#ifndef PARTICLENEIGHBOURS_H
#define PARTICLENEIGHBOURS_H

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

    std::vector<float> add (const std::vector<float> vals) const {

        const size_t res_size = neighId.size();
        std::vector<float> res (res_size);

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t nid = neighId[ctr];
            const size_t oid =  origId[ctr];

            res[ctr] = {vals[oid] + vals[nid]};
        }
        return res;
    }

    std::vector<Vector> add (const std::vector<Vector> vals) const {

        const size_t res_size = neighId.size();
        std::vector<Vector> res (res_size);

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t nid = neighId[ctr];
            const size_t oid =  origId[ctr];

            res[ctr] = {vals[oid] + vals[nid]};
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

    std::vector<float> get_b (const std::vector<float>& vals) const {

        const size_t res_size = neighId.size();

        // TODO fix this RVO
        std::vector<float> ret (res_size, 0.0);

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t nid = neighId[ctr];

            ret[ctr] = vals[nid];

        }

        return ret;
    }

    void sum (const std::vector<float>& vals, std::vector<float>& ret) const {

        // TODO ASSERT that vals has same size as neighId
        const size_t res_size = neighId.size();

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t oid =  origId[ctr];

            // TODO neighId[ctr] is the id of the  neighbouring cells but not
            // the position in vals
            // const size_t nid = neighId[ctr];

            // std::cout
            //     << "[DEBUG] nid " << nid
            //     << " oid " << oid
            //     << " vals[oid] " << vals[oid]
            //     << " ret[oid] " << ret[oid]
            //     << std::endl;

            ret[oid] += vals[ctr];
        }
    }

    void sum (const std::vector<Vector>& vals, std::vector<Vector>& ret) const {

        const size_t res_size = neighId.size();

        for (size_t ctr = 0; ctr < res_size; ctr++) {

            const size_t oid =  origId[ctr];
            // const size_t nid = neighId[ctr];

            // std::cout
            //     << "[DEBUG] nid " << nid
            //     << " oid " << oid
            //     << " vals[oid] " << vals[oid]
            //     << " ret[oid] " << ret[oid]
            //     << std::endl;

            ret[oid] += vals[ctr];
        }
    }
};


struct LeanParticleNeighbours {

    // TODO This should be sorted hence replace by number of entries per
    // particle
    // std::vector<size_t> origId;

    // std::vector<size_t> first;
    // std::vector<size_t> last;


    std::vector<size_t> nElems;
    std::vector<size_t> neighId;

    // std::vector<Vector> distances;

    // std::vector<Vector> normalised_distances;

    // std::vector<float> squared_length;
    //
    // size_t get_neighbbour_startId(size_t pid) {
    //
    // };
    //
    // size_t get_neighbbour_endId(size_t pid) {
    //
    // };

    // std::vector<Vector> subtract (const std::vector<Vector> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<Vector> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] - vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<float> add (const std::vector<float> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<float> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] + vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<Vector> add (const std::vector<Vector> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<Vector> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] + vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<float> subtract (const std::vector<float> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<float> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] - vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<float> get_b (const std::vector<float>& vals) const {
    //
    //     const size_t res_size = neighId.size();
    //
    //     // TODO fix this RVO
    //     std::vector<float> ret (res_size, 0.0);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //
    //         ret[ctr] = vals[nid];
    //
    //     }
    //
    //     return ret;
    // }
    //
    // void sum (const std::vector<float>& vals, std::vector<float>& ret) const {
    //
    //     // TODO ASSERT that vals has same size as neighId
    //     const size_t res_size = neighId.size();
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t oid =  origId[ctr];
    //
    //         // TODO neighId[ctr] is the id of the  neighbouring cells but not
    //         // the position in vals
    //         // const size_t nid = neighId[ctr];
    //
    //         // std::cout
    //         //     << "[DEBUG] nid " << nid
    //         //     << " oid " << oid
    //         //     << " vals[oid] " << vals[oid]
    //         //     << " ret[oid] " << ret[oid]
    //         //     << std::endl;
    //
    //         ret[oid] += vals[ctr];
    //     }
    // }
    //
    // void sum (const std::vector<Vector>& vals, std::vector<Vector>& ret) const {
    //
    //     const size_t res_size = neighId.size();
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t oid =  origId[ctr];
    //         // const size_t nid = neighId[ctr];
    //
    //         // std::cout
    //         //     << "[DEBUG] nid " << nid
    //         //     << " oid " << oid
    //         //     << " vals[oid] " << vals[oid]
    //         //     << " ret[oid] " << ret[oid]
    //         //     << std::endl;
    //
    //         ret[oid] += vals[ctr];
    //     }
    // }
};


struct ParticleNeighboursCGALFree {

    // TODO This should be sorted hence replace by number of entries per
    // particle
    std::vector<size_t> origId;

    std::vector<size_t> neighId;

    std::vector<Point3D> distances;

    std::vector<Point3D> normalised_distances;

    std::vector<float> squared_length;

    // std::vector<Point3D> subtract (const std::vector<Point3D> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<Vector> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] - vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<float> add (const std::vector<float> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<float> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] + vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<Vector> add (const std::vector<Vector> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<Vector> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] + vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<float> subtract (const std::vector<float> vals) const {
    //
    //     const size_t res_size = neighId.size();
    //     std::vector<float> res (res_size);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //         const size_t oid =  origId[ctr];
    //
    //         res[ctr] = {vals[oid] - vals[nid]};
    //     }
    //     return res;
    // }
    //
    // std::vector<float> get_b (const std::vector<float>& vals) const {
    //
    //     const size_t res_size = neighId.size();
    //
    //     // TODO fix this RVO
    //     std::vector<float> ret (res_size, 0.0);
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t nid = neighId[ctr];
    //
    //         ret[ctr] = vals[nid];
    //
    //     }
    //
    //     return ret;
    // }
    //
    // void sum (const std::vector<float>& vals, std::vector<float>& ret) const {
    //
    //     // TODO ASSERT that vals has same size as neighId
    //     const size_t res_size = neighId.size();
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t oid =  origId[ctr];
    //
    //         // TODO neighId[ctr] is the id of the  neighbouring cells but not
    //         // the position in vals
    //         // const size_t nid = neighId[ctr];
    //
    //         // std::cout
    //         //     << "[DEBUG] nid " << nid
    //         //     << " oid " << oid
    //         //     << " vals[oid] " << vals[oid]
    //         //     << " ret[oid] " << ret[oid]
    //         //     << std::endl;
    //
    //         ret[oid] += vals[ctr];
    //     }
    // }
    //
    // void sum (const std::vector<Vector>& vals, std::vector<Vector>& ret) const {
    //
    //     const size_t res_size = neighId.size();
    //
    //     for (size_t ctr = 0; ctr < res_size; ctr++) {
    //
    //         const size_t oid =  origId[ctr];
    //         // const size_t nid = neighId[ctr];
    //
    //         // std::cout
    //         //     << "[DEBUG] nid " << nid
    //         //     << " oid " << oid
    //         //     << " vals[oid] " << vals[oid]
    //         //     << " ret[oid] " << ret[oid]
    //         //     << std::endl;
    //
    //         ret[oid] += vals[ctr];
    //     }
    // }
};

#endif
