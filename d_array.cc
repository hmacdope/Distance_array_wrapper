#include "distance_array_comparison.h"


int main()
{



    // setup
    constexpr uint64_t N = 10000;
    constexpr bool debug = false;
    constexpr bool print_result = false;

    float coords1[3 * N];
    float coords2[3 * N];

    double result1[N * N] = {0};
    double result2[N * N] = {0};
    double result3[N * N] = {0};
    double result4[N * N] = {0};

    // classes that mock atomgroup
    auto ag_mock1 = AGWrapper(N);
    auto ag_mock2 = AGWrapper(N);

    auto float_mock1 = FloatWrapper(N);
    auto float_mock2 = FloatWrapper(N);

    for (uint64_t i = 0; i < 3 * N; i++)
    {
        coords1[i] = static_cast<float>(i);
        coords2[i] = static_cast<float>(i);
    }

    if (debug)
    {
        for (uint64_t i = 0; i < 3 * N; i++)
        {
            printf(" %f \n", coords1[i]);
            printf(" %f \n", coords2[i]);
            printf(" %f \n", ag_mock1.coords[i]);
        }
    }

    // raw MDA style
    _calc_distance_array(coords1, N, coords2, N, result1);

    // wrapped version, two FloatWrapper arguments
    DistanceArray(float_mock1, float_mock2, result2);

    // wrapped version mixed args
    DistanceArray(float_mock1, ag_mock1, result3);

    // wrapped version two AGwrapper arguments
    DistanceArray(ag_mock1, ag_mock2, result4);

    if (print_result)
    {
        print_square_mat(result1, N, "raw mda");
        print_square_mat(result2, N, "wrapped float args");
        print_square_mat(result3, N, "wrapped mixed_args");
        print_square_mat(result4, N, "wrapped AG args");
    }
}