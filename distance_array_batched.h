#pragma once
#include "distance_array_comparison.h"
#include <algorithm>

template <typename T, typename U>
void DistanceArrayBatched(T ref, U conf, double *distances, uint64_t batchsize)
{
    float ref_buffer[batchsize];
    float conf_buffer[batchsize];
    double result_buffer[batchsize];

    uint64_t nref = ref.N;
    uint64_t nconf = conf.N;
    uint64_t bsize_ref = std::min(3 *nref, batchsize);
    uint64_t bsize_conf = std::min(3* nconf, batchsize);
    ref.preload_external(ref_buffer, bsize_ref); // make sure only fill up to number of vals
    conf.preload_external(conf_buffer, bsize_conf);

    uint64_t i, j, b;
    double dx[3];
    double rsq;

    // avoid dealing with coordinate in type signature
    // this should be a reinterpret cast so cheap
    // coordinate *conf_ = (coordinate *)conf_buffer;
    // coordinate *ref_ = (coordinate *)ref_buffer;

    // for (i = 0; i < niter_ref; i++)
    // {
    //     for (j = 0; j < niter_conf; j++)
    //     {
    //         dx[0] = conf_[j][0] - ref_[i][0];
    //         dx[1] = conf_[j][1] - ref_[i][1];
    //         dx[2] = conf_[j][2] - ref_[i][2];
    //         rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
    //         *(distances + i * niter_conf + j) = sqrt(rsq);
    //     }
    // }
}