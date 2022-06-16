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

    for (iter_ref = 0; iter_ref < nref - ref_overhang; iter_ref += bsize_ref)
    {
        // printf("ref preload iter ref %i \n", iter_ref);
        ref.preload_external(ref_buffer, bsize_ref);

        for (iter_conf = 0; iter_conf < nconf - conf_overhang; iter_conf += bsize_conf)
        {
            // printf("conf preload iter conf %i\n", iter_conf);
            conf.preload_external(conf_buffer, bsize_conf);

            for (i = 0; i < bsize_ref; i++)
            {
                for (j = 0; j < bsize_conf; j++)
                {
                    dx[0] = conf_buffer[3 * j] - ref_buffer[3 * i];
                    dx[1] = conf_buffer[3 * j + 1] - ref_buffer[3 * i + 1];
                    dx[2] = conf_buffer[3 * j + 2] - ref_buffer[3 * i + 2];
                    rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                    *(distances + iter_ref * nconf + iter_conf + i * nconf + j) = sqrt(rsq);
                    // printf("mem loc %ld ", iter_ref * nconf + iter_conf + i * nconf + j);
                }
            }
        }

        conf.reset_external_buffer_iteration();
    }
    // printf("\n\n OVERHANG TIME\n");
    // now deal with overhang
    uint64_t gcd_conf = std::gcd(nconf, bsize_conf);
    uint64_t gcd_ref = std::gcd(nref, bsize_ref);
    // printf("iter_conf %ld\n ", iter_conf);
    // printf("iter_ref %ld\n", iter_ref);
    // printf("gcd conf %ld\n", gcd_conf);
    // printf("gcd conf %ld\n", gcd_ref);

    // printf("REF OVERHANG\n");
    // strided in this dimension
    ref.preload_external(ref_buffer, ref_overhang);

    for (int i = 0; i < nconf; i += gcd_conf)
    {
        conf.preload_external(conf_buffer, gcd_conf);
        for (int ii = 0; ii < ref_overhang; ii++)
        {
            for (int jj = 0; jj < gcd_conf; jj++)
            {
                // printf("pair %f %f \n %f %f \n %f %f\n",conf_buffer[3 * jj],  ref_buffer[3 * ii],  conf_buffer[3 * ii +1],  ref_buffer[3 * ii +1],   conf_buffer[3 * ii +2],  ref_buffer[3 * ii +2]   );
                dx[0] = conf_buffer[3 * jj] - ref_buffer[3 * ii];
                dx[1] = conf_buffer[3 * jj + 1] - ref_buffer[3 * ii + 1];
                dx[2] = conf_buffer[3 * jj + 2] - ref_buffer[3 * ii + 2];
                rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                // printf("dist %f\n", sqrt(rsq));
                *(distances + iter_ref * nconf + i + ii * nconf + jj) = sqrt(rsq);
                // printf("mem loc %ld \n", iter_ref * nconf + i + ii * nconf + jj);
            }
        }
    }
    // printf("CONF OVERHANG\n");
    // contiguous in this dimension
    conf.reset_external_buffer_iteration();
    ref.reset_external_buffer_iteration();
    conf.seek(nconf - conf_overhang);

    conf.preload_external(conf_buffer, conf_overhang);

    for (int j = 0; j < nref; j += gcd_ref)
    {
        ref.preload_external(ref_buffer, gcd_ref);
        for (int jj = 0; jj < conf_overhang; jj++)
        {
            for (int ii = 0; ii < gcd_ref; ii++)
            {
                dx[0] = conf_buffer[3 * jj] - ref_buffer[3 * ii];
                dx[1] = conf_buffer[3 * jj + 1] - ref_buffer[3 * ii + 1];
                dx[2] = conf_buffer[3 * jj + 2] - ref_buffer[3 * ii + 2];
                rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
                // printf("dist %f\n", sqrt(rsq));
                *(distances + iter_conf +  j*nref + jj ) = sqrt(rsq);
                // printf("mem loc %ld \n", iter_conf + j*nref + jj );
            }
        }
    }
}

//             *(distances + (nref - ref_overhang) * nconf + i * nconf + j) = sqrt(rsq);