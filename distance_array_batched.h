#pragma once
#include "distance_array_comparison.h"
#include <algorithm>




template <typename T, typename U>
void DistanceArrayBatched(T ref, U conf, double *distances, uint64_t batchsize)
{
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

    if (nref % bsize_ref | nconf % bsize_conf) // overhang in either dimension?
    {
        printf("overhangs!!\n");
        ref.preload_external(ref_buffer, bsize_ref);
        conf.preload_external(conf_buffer, bsize_conf);
        iter_ref += nref % bsize_ref;
        iter_conf += nconf % bsize_conf;
    }

    for (; iter_ref < nref; iter_ref += bsize_ref)
    {
        printf("ref preload \n");
        ref.preload_external(ref_buffer, bsize_ref);

        for (; iter_conf < nconf; iter_conf += bsize_conf)
        {
            printf("conf preload\n");
            conf.preload_external(conf_buffer, bsize_conf);
        }

        conf.reset_external_buffer_iteration();
        iter_conf = 0;
    }
}
