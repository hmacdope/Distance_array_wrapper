#include <cmath>
#include <cstdint>
#include <cstdio>

typedef float coordinate[3];


void _calc_distance_array(coordinate *ref, uint64_t numref, coordinate *conf,
                                 uint64_t numconf, double *distances)
{
    uint64_t i, j;
    double dx[3];
    double rsq;

    for (i = 0; i < numref; i++)
    {
        for (j = 0; j < numconf; j++)
        {
            dx[0] = conf[j][0] - ref[i][0];
            dx[1] = conf[j][1] - ref[i][1];
            dx[2] = conf[j][2] - ref[i][2];
            rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
            *(distances + i * numconf + j) = sqrt(rsq);
        }
    }
}

int main()
{
    constexpr int N = 20;
    constexpr bool debug = false;
    constexpr bool print_result = true;



    float coords1[3 * N];
    float coords2[3 * N];
    double result[N*N] = {0};

    for (int i = 0; i < 3 * N; i++)
    {
        coords1[i] = static_cast<float>(i);
        coords2[i] = static_cast<float>(i);
    }

    if (debug)
    {
        for (int i = 0; i < 3 * N; i++)
        {
            printf(" %f \n", coords1[i]);
            printf(" %f \n", coords2[i]);
        }
    }

    // c style cast but oh well
    _calc_distance_array((coordinate*)coords1, N, (coordinate*)coords2, N, result);

    if (print_result)
    {
        for (int i = 0; i < N; i++)
        {
            printf(" %f \n", result[i]);
        }
    }
}