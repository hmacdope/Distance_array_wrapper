#include "distance_array_comparison.h"
#include "distance_array_batched.h"
#include <gtest/gtest.h>

int main()
{

    // setup
    constexpr uint64_t N = 5;
    constexpr uint64_t M = 3;
    constexpr bool debug = true;
    constexpr bool print_result = true;

    constexpr int bufsize = 2; // IN ATOMS

    float buffer[3 * bufsize];

    float coords1[3 * N];
    float coords2[3 * M];

    double result[N * M] = {0};
    double result2[N * M] = {0};

    // classes that mock atomgroup
    auto ag_mock1 = AGWrapper(N);
    auto ag_mock2 = AGWrapper(M);

    auto float_mock1 = FloatWrapper(N);
    auto float_mock2 = FloatWrapper(M);

    for (uint64_t i = 0; i < 3 * N; i++)
    {
        coords1[i] = static_cast<float>(i);
    }

    for (uint64_t i = 0; i < 3 * M; i++)
    {
        coords2[i] = static_cast<float>(i);
    }

    if (debug)
    {
        for (uint64_t i = 0; i < 3 * N; i++)
        {
            printf(" %f \n", ag_mock1.coords[i]);
        }

        // ag_mock1.preload_external(buffer, bufsize);
        // for (uint64_t i = 0; i < bufsize; i++)
        //     printf(" %f \n", buffer[i]);
    }

    // raw MDA style
    _calc_distance_array(coords1, N, coords2, M, result);
    if (print_result)
    {
        print_rect_mat(result, N, M, "raw mda");
    }

    DistanceArrayBatched(float_mock1, float_mock2, result2, bufsize);
    if (print_result)
    {
        print_rect_mat(result2, N, M, "batched");
    }

    for (int i = 0; i < M * N; i++)
    {
        EXPECT_FLOAT_EQ(result[i], result2[i]);
    }
}