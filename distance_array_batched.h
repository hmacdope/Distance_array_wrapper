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

    // slight indirection required here to homogenize interface between
    // _AtomGroupIterator and _ArrayIterator where passing stack allocated
    // array as float*& does not decay to a float* as desired.
    float *ref_buffer_ = ref_buffer;
    float *conf_buffer_ = conf_buffer;

    uint64_t nref = ref.N;
    uint64_t nconf = conf.N;
    uint64_t bsize_ref = std::min(nref, batchsize); // is  our batchsize larger than the number of coords?
    uint64_t bsize_conf = std::min(nconf, batchsize);

    // printf("blocksize ref %i blocksize conf %i\n\n", bsize_conf, bsize_ref);

    uint64_t iter_ref = 0;
    uint64_t iter_conf = 0;
    uint64_t i, j;
    double rsq;
    float dx[3];
    int ref_overhang = nref % bsize_ref;
    int conf_overhang = nconf % bsize_conf;

    // printf("MAIN LOOP\n");
    for (iter_ref = 0; iter_ref < nref - ref_overhang; iter_ref += bsize_ref)
    {
        // printf("ref preload iter ref %i \n", iter_ref);
        ref.preload_external(ref_buffer_, bsize_ref);

        for (iter_conf = 0; iter_conf < nconf - conf_overhang; iter_conf += bsize_conf)
        {
            // printf("conf preload iter conf %i\n", iter_conf);
            conf.preload_external(conf_buffer_, bsize_conf);

            for (i = 0; i < bsize_ref; i++)
            {
                for (j = 0; j < bsize_conf; j++)
                {
                    dx[0] = conf_buffer_[3 * j] - ref_buffer_[3 * i];
                    dx[1] = conf_buffer_[3 * j + 1] - ref_buffer_[3 * i + 1];
                    dx[2] = conf_buffer_[3 * j + 2] - ref_buffer_[3 * i + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_ref * nconf + iter_conf + i * nconf + j) = sqrt(rsq);
                    // printf("mem loc %ld ", iter_ref * nconf + iter_conf + i * nconf + j);
                }
            }
        }

        conf.reset_external_buffer_iteration();
    }
    // printf("\n\n OVERHANG TIME\n");
    uint64_t gcd_conf = std::gcd(nconf, bsize_conf);
    uint64_t gcd_ref = std::gcd(nref, bsize_ref);
    // printf(" nref %ld\n", nref);
    // printf(" nconf %ld\n", nconf);
    // printf("iter_conf %ld\n ", iter_conf);
    // printf("iter_ref %ld\n", iter_ref);
    // printf("gcd conf %ld\n", gcd_conf);
    // printf("gcd ref %ld\n", gcd_ref);
    // printf("ref_overhang %ld\n", ref_overhang);
    // printf("conf_overhang %ld\n", conf_overhang);
    if (ref_overhang)
    {
        // printf("REF OVERHANG\n");
        ref.preload_external(ref_buffer_, ref_overhang);

        for (int j = 0; j < nconf; j += gcd_conf)
        {
            conf.preload_external(conf_buffer_, gcd_conf);

            for (int i = 0; i < ref_overhang; i++)
            {

                for (int k = 0; k < gcd_conf; k++)
                {
                    dx[0] = conf_buffer_[3 * k] - ref_buffer_[3 * i];
                    dx[1] = conf_buffer_[3 * k + 1] - ref_buffer_[3 * i + 1];
                    dx[2] = conf_buffer_[3 * k + 2] - ref_buffer_[3 * i + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_ref * nconf + i * nconf + j + k) = sqrt(rsq);
                    // printf(" mem loc %ld\n", iter_ref * nconf + i * nconf + j + k);
                }
            }
        }
    }
    if (conf_overhang)
    {
        // }
        // printf("CONF OVERHANG\n");
        // contiguous in this dimension
        conf.reset_external_buffer_iteration();
        ref.reset_external_buffer_iteration();
        conf.seek(nconf - conf_overhang);

        conf.preload_external(conf_buffer_, conf_overhang);
        for (int j = 0; j < nref; j += gcd_ref)
        {
            ref.preload_external(ref_buffer_, gcd_ref);
            for (int i = 0; i < conf_overhang; i++)
            {
                for (int k = 0; k < gcd_ref; k++)
                {
                    dx[0] = conf_buffer_[3 * i] - ref_buffer_[3 * k];
                    dx[1] = conf_buffer_[3 * i + 1] - ref_buffer_[3 * k + 1];
                    dx[2] = conf_buffer_[3 * i + 2] - ref_buffer_[3 * k + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_conf + j * nconf + i + k * nconf) = sqrt(rsq);
                    // printf("mem loc %ld\n", iter_conf + j * nconf + i + k * nconf);
                }
            }
        }
    }
}

//             *(distances + (nref - ref_overhang) * nconf + i * nconf + j) = sqrt(rsq);