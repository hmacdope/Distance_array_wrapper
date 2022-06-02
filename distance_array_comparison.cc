#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

typedef float coordinate[3];

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

class AtomGroupMock
{
public:
    int N;
    std::vector<float> ix;
    std::vector<float> coords;

    AtomGroupMock(int N_)
    {
        printf("constructed\n");
        N = N_;
        for (int i = 0; i < N; i++)
        {
            ix.push_back(static_cast<float>(i));
        }
        for (int i = 0; i < 3 * N; i++)
        {
            coords.push_back(static_cast<float>(i));
        }
    }
};

template <typename T>
struct IsAg
{
    static constexpr bool value = false;
};
template <>
struct IsAg<AtomGroupMock>
{
    static constexpr bool value = true;
};

template <typename T>
struct IsArr
{
    static constexpr bool value = false;
};
template <>
struct IsArr<float *>
{
    static constexpr bool value = true;
};

template <typename T>
float *_wraps_ag(T inp)
{
    return inp;
}

template <>
float * _wraps_ag(AtomGroupMock inp)
{
    return inp.coords.data();
}

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

template <typename T, typename U>
void DistanceArray(T ref, uint64_t numref, U conf, uint64_t numconf, double *distances)
{

    auto ref_ = _wraps_ag<T>(ref);
    auto conf_ = _wraps_ag<U>(conf);

    _calc_distance_array((coordinate *)ref_, numref, (coordinate *)conf_, numconf, distances);
}

int main()
{
    constexpr int N = 10;
    constexpr bool debug = false;
    constexpr bool print_result = true;

    float coords1[3 * N];
    float coords2[3 * N];

    double result1[N * N] = {0};
    double result2[N * N] = {0};
    double result3[N * N] = {0};

    auto ag_mock1 = AtomGroupMock(N);
    auto ag_mock2 = AtomGroupMock(N);

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
    _calc_distance_array((coordinate *)coords1, N, (coordinate *)coords2, N, result1);

    // wrapped version, float arguments

    DistanceArray(coords1, N, coords2, N, result2);

    // wrapped version mixed args
    DistanceArray(coords1, N, ag_mock1, N, result3);

    if (print_result)
    {
        print_square_mat(result1, N, "raw mda");
        print_square_mat(result2, N, "wrapped float args");
        print_square_mat(result3, N, "mixed_args");
    }
}