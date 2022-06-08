#include "distance_array_comparison.h"
#include "distance_array_batched.h"

int main()
{

    // setup
    constexpr uint64_t N = 10;
    constexpr bool debug = true;
    constexpr bool print_result = true;

    constexpr int bufsize = 3; // IN ATOMS 

    float buffer[3* bufsize];

    float coords1[3 * N];
    float coords2[3 * N];

    double result[N * N] = {0};
    double result2[N * N] = {0};


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
            printf(" %f \n", ag_mock1.coords[i]);
        }

        // ag_mock1.preload_external(buffer, bufsize);
        // for (uint64_t i = 0; i < bufsize; i++)
        //     printf(" %f \n", buffer[i]);

    }

    // raw MDA style
    _calc_distance_array(coords1, N, coords2, N, result);
    if (print_result)
    {
        print_square_mat(result, N, "raw mda");
    }

    DistanceArrayBatched(ag_mock1, ag_mock2, result2, bufsize);
    if (print_result)
    {
        print_square_mat(result2, N, "batched");
    }
}