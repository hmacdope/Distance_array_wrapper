#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <tuple>

// from mda
typedef float coordinate[3];

// helper
template <typename T>
void print_square_mat(T *buffer, int buf_len, std::string tag)
{
    printf(" %s \n", tag.c_str());
    for (int i = 0; i < buf_len; i++)
    {
        for (int j = 0; j < buf_len; j++)
        {
            printf(" %f ", buffer[buf_len * i + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// mocks AtomGroup
class AtomGroupMock
{
public:
    int N;
    std::vector<int> ix;
    std::vector<float> coords;
    int i = 0;

    AtomGroupMock(int N) : N(N)
    {
        for (int i = 0; i < N; i++)
        {
            ix.push_back(static_cast<float>(i));
        }
        for (int i = 0; i < 3 * N; i++)
        {
            coords.push_back(static_cast<float>(i));
        }
    }

    float *next()
    {
        if (i  < N)
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
    int N;
    std::vector<float> coords;
    int i = 0;

    FloatWrapper(int N) : N(N)
    {
        for (int i = 0; i < 3 * N; i++)
        {
            coords.push_back(static_cast<float>(i));
        }
    }

    float *next()
    {
        if (i  < N)
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

// wrapper for plain float *
template <typename T>
T *iter_coords(T *inp)
{
    inp += 3;
    return inp + 3;
}

// wrapper for AtomGroupMock, MUST BE PASSED BY REF, otherwise new one constructed and
// ref invalid upon exit
auto iter_coords(AtomGroupMock &inp)
{
    float *pointer = inp.next();
    return pointer;
}

// wrapper for FloatWrapper, MUST BE PASSED BY REF, otherwise new one constructed and
// ref invalid upon exit
auto iter_coords(FloatWrapper &inp)
{
    float *pointer = inp.next();
    return pointer;
}


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
void DistanceArray(T ref, uint64_t numref, U conf, uint64_t numconf, double *distances)
{

    uint64_t i, j;
    double dx[3];
    double rsq;


    for (i = 0; i < numref; i++)
    {
        auto ref__ = iter_coords(ref);

        for (j = 0; j < numconf; j++)
        {
            auto conf__ = iter_coords(conf);

            dx[0] = conf__[0] - ref__[0];
            dx[1] = conf__[1] - ref__[1];
            dx[2] = conf__[2] - ref__[2];
            rsq = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
            *(distances + i * numconf + j) = sqrt(rsq);
        }
    }
}

int main()
{
    // setup
    constexpr int N = 10;
    constexpr bool debug = false;
    constexpr bool print_result = true;

    float coords1[3 * N];
    float coords2[3 * N];

    double result1[N * N] = {0};
    double result2[N * N] = {0};
    double result3[N * N] = {0};
    double result4[N * N] = {0};

    // classes that mock atomgroup
    auto ag_mock1 = AtomGroupMock(N);
    auto ag_mock2 = AtomGroupMock(N);

    auto float_mock1 = FloatWrapper(N);
    auto float_mock2 = FloatWrapper(N);

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
            printf(" %f \n", ag_mock1.coords[i]);
        }
    }

    // raw MDA style
    _calc_distance_array(coords1, N, coords2, N, result1);

    // wrapped version, float arguments

    // DistanceArray(coords1, N, coords2, N, result2);

    // //wrapped version mixed args
    // DistanceArray(coords1, N, ag_mock1, N, result3);

    // wrapped version two Atomgroups
    DistanceArray(ag_mock1, N, ag_mock2, N, result4);

    if (print_result)
    {
        print_square_mat(result1, N, "raw mda");
        print_square_mat(result2, N, "wrapped float args");
        print_square_mat(result3, N, "wrapped mixed_args");
        print_square_mat(result4, N, "wrapped AG args");
    }
}