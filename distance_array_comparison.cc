#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <typeinfo>

// from mda
typedef float coordinate[3];

// helper
template <typename T>
void print_square_mat(T *buffer, uint64_t buf_len, std::string tag)
{
    printf(" %s \n", tag.c_str());
    for (uint64_t i = 0; i < buf_len; i++)
    {
        for (uint64_t j = 0; j < buf_len; j++)
        {
            printf(" %f ", buffer[buf_len * i + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// mocks AtomGroup
class AGWrapper
{
public:
    uint64_t N;
    std::vector<uint64_t> ix;
    std::vector<float> coords;
    uint64_t i = 0;

    AGWrapper(uint64_t N) : N(N)
    {
        ix.reserve(N);
        coords.reserve(3 * N);
        for (uint64_t i = 0; i < N; i++)
        {
            ix.push_back(static_cast<float>(i));
        }
        for (uint64_t i = 0; i < 3 * N; i++)
        {
            coords.push_back(static_cast<float>(i));
        }
    }

    float *next()
    {
        if (i < N)
        {
            i += 1;
            float *ptr = coords.data() + 3 * (i - 1);
            return ptr;
        }
        else
        {
            reset_iteration();
            return next();
        }
    }

    void reset_iteration()
    {
        i = 0;
    }
};

class FloatWrapper
{

public:
    uint64_t N;
    std::vector<float> coords;
    uint64_t i = 0;

    FloatWrapper(uint64_t N) : N(N)
    {
        coords.reserve(3 * N);
        for (uint64_t i = 0; i < 3 * N; i++)
        {
            coords.push_back(static_cast<float>(i));
        }
    }

    float *next()
    {
        if (i < N)
        {
            i += 1;
            float *ptr = coords.data() + 3 * (i - 1);
            return ptr;
        }
        else
        {
            reset_iteration();
            return next();
        }
    }

    void reset_iteration()
    {
        i = 0;
    }
};

// directly from MDA
void _calc_distance_array(float *ref, uint64_t numref, float *conf,
                          uint64_t numconf, double *distances)
{
    uint64_t i, j;
    double dx[3];
    double rsq;

    // avoid dealing with coordinate in type signature
    // this should be a reinterpret cast so cheap
    coordinate *conf_ = (coordinate *)conf;
    coordinate *ref_ = (coordinate *)ref;

    for (i = 0; i < numref; i++)
    {
        for (j = 0; j < numconf; j++)
        {
            dx[0] = conf_[j][0] - ref_[i][0];
            dx[1] = conf_[j][1] - ref_[i][1];
            dx[2] = conf_[j][2] - ref_[i][2];
            rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
            *(distances + i * numconf + j) = sqrt(rsq);
        }
    }
}

template <typename T, typename U>
void DistanceArray(T ref, U conf, double *distances)
{

    uint64_t i, j;
    double dx[3];
    double rsq;
    float *ref_;
    float *conf_;

    for (i = 0; i < ref.N; i++)
    {
        ref_ = ref.next();
        for (j = 0; j < conf.N; j++)
        {
            conf_ = conf.next();

            dx[0] = conf_[0] - ref_[0];
            dx[1] = conf_[1] - ref_[1];
            dx[2] = conf_[2] - ref_[2];
            rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
            *(distances + i * conf.N + j) = sqrt(rsq);
        }
    }
}

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