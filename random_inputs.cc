#include <random>
#include <gtest/gtest.h>
#include "distance_array_comparison.h"
#include "distance_array_batched.h"

int main()
{
    const bool print_result=false;
    constexpr int N_test = 1000;
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<uint64_t> dist(1, 10); // distribution in range [1, 6]

    for (int i = 0; i < N_test; i++)
    {

        uint64_t N = dist(rng);
        uint64_t M = dist(rng);
        int batchsize = dist(rng);

        printf("N  %ld M %ld batchsize %i  \n", N, M, batchsize);

        float coords1[3 * N];
        float coords2[3 * M];

        double result1[N * M] = {0};
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

        // raw MDA style
        _calc_distance_array(coords1, N, coords2, M, result1);
        if (print_result)
        {
            print_rect_mat(result1, N, M, "raw mda");
        }

        DistanceArrayBatched(ag_mock1, ag_mock2, result2, batchsize);
        if (print_result)
        {
            print_rect_mat(result2, N, M, "batched");
        }

        for(int i=0; i<M*N; i++) {
            EXPECT_FLOAT_EQ(result1[i], result2[i]);
        }
    }
}