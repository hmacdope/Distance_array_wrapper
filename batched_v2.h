#pragma once
#include "distance_array_comparison.h"
#include <algorithm>
#include <numeric>

template <typename T, typename U>
void DistanceArrayBatched(T ref, U conf, double *distances, uint64_t batchsize)
{
    // blocked algorithm with a tile repeat for modular overhang
    // printf("batched distance array\n");
    const uint64_t atom_bufsize = 3 * batchsize;
    float ref_buffer[atom_bufsize];
    float conf_buffer[atom_bufsize];
    double result_buffer[atom_bufsize];

    double rsq;
    float dx[3];

    uint64_t nref = ref.N;
    uint64_t nconf = conf.N;
    uint64_t bsize_ref = std::min(nref, batchsize); // is  our batchsize larger than the number of coords?
    uint64_t bsize_conf = std::min(nconf, batchsize);
    uint64_t local_i = 0;
    uint64_t local_j = 0;
    uint64_t jcount = 0;

    // printf("blocksize ref %i blocksize conf %i\n\n", bsize_conf, bsize_ref);
    printf("bsize_conf %ld \n", bsize_conf);
    printf("bsize_ref %ld \n", bsize_ref);
    printf("nconf  %ld \n", nconf);
    printf("nref %ld \n", nref);

    for (int i = 0; i < nref; ++i)
    {
        // printf("\n i  %ld\n", i);

        if (!(i % bsize_ref))
        {
            // printf("\n load_ref \n");
            ref.preload_external(ref_buffer, bsize_ref, i);
            local_i = 0;
        }

        for (int j = 0; j < nconf; ++j)
        {
            // printf("j  %ld\n", j);
            if (!(j % bsize_conf))
            {
                // printf("load_conf \n");

                conf.preload_external(conf_buffer, bsize_conf, j);
                local_j = 0;
            }

            // printf("local j %ld \n", local_j);

            dx[0] = conf_buffer[3 * local_j] - ref_buffer[3 * local_i];
            dx[1] = conf_buffer[3 * local_j + 1] - ref_buffer[3 * local_i + 1];
            dx[2] = conf_buffer[3 * local_j + 2] - ref_buffer[3 * local_i + 2];
            rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
            *(distances + nconf * i + j) = sqrt(rsq);

            local_j += 1;
        }
        local_i += 1;
        // printf("local i %ld \n", local_i);
    }
}
