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
    float dx[3];
    int ref_overhang = nref % bsize_ref;
    int conf_overhang = nconf % bsize_conf;

    // if (ref_overhang | conf_overhang) // overhang in either dimension?
    // {
    //     iter_ref -= ref_overhang
    //     iter_conf -= 
    // }

    for (; iter_ref < nref - ref_overhang; iter_ref += bsize_ref)
    {
        printf("ref preload \n");
        ref.preload_external(ref_buffer, bsize_ref);

        for (; iter_conf < nconf - conf_overhang; iter_conf += bsize_conf)
        {
            printf("conf preload\n");
            conf.preload_external(conf_buffer, bsize_conf);
            _calc_distance_array(ref_buffer, bsize_ref, conf_buffer, bsize_conf, result_buffer);
            distances += nconf - bsize_conf +conf_overhang;
            print_square_mat(result_buffer, bsize_ref, "batched");
        }


        conf.reset_external_buffer_iteration();
        iter_conf = 0;
    }
}
