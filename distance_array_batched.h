#pragma once
#include "distance_array_comparison.h"
#include <algorithm>

template <typename T, typename U>
void DistanceArrayBatched(T ref, U conf, double *distances, uint64_t batchsize)
{
    // blocked algorithm with a tile repeat for modular overhang
    printf("batched distance array\n");
    const uint64_t atom_bufsize = 3 * batchsize;
    float ref_buffer[atom_bufsize];
    float conf_buffer[atom_bufsize];
    double result_buffer[atom_bufsize];

    uint64_t nref = ref.N;
    uint64_t nconf = conf.N;
    uint64_t bsize_ref = std::min(nref, batchsize); // is  our batchsize larger than the number of coords?
    uint64_t bsize_conf = std::min(nconf, batchsize);

    printf("blocksize ref %i blocksize conf %i\n\n", bsize_conf, bsize_ref);

    uint64_t iter_ref = 0;
    uint64_t iter_conf = 0;
    uint64_t i, j;
    double rsq;
    float dx[3];
    int ref_overhang = nref % bsize_ref;
    int conf_overhang = nconf % bsize_conf;


    for (; iter_ref < nref - ref_overhang; iter_ref += bsize_ref)
    {
        printf("ref preload \n");
        ref.preload_external(ref_buffer, bsize_ref);

        for (; iter_conf < nconf - conf_overhang; iter_conf += bsize_conf)
        {
            printf("conf preload\n");
            conf.preload_external(conf_buffer, bsize_conf);

            // avoid dealing with coordinate in type signature
            // this should be a reinterpret cast so cheap
            coordinate *conf_ = (coordinate *)conf_buffer;
            coordinate *ref_ = (coordinate *)ref_buffer;

            for (i = 0; i < bsize_ref; i++)
            {
                for (j = 0; j < bsize_conf; j++)
                {
                    dx[0] = conf_[j][0] - ref_[i][0];
                    dx[1] = conf_[j][1] - ref_[i][1];
                    dx[2] = conf_[j][2] - ref_[i][2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + i * nconf + j) = sqrt(rsq);
                }
            }
        }

        conf.reset_external_buffer_iteration();
        iter_conf = 0;
    }
}
