#include <benchmark/benchmark.h>
#include "distance_array_comparison.h"

class Distance : public benchmark::Fixture
{
public:
    // setup
    constexpr static uint64_t N = 1000;

    float coords1[3 * N];
    float coords2[3 * N];

    double result[N * N] = {0};

    // classes that mock atomgroup
    AGWrapper ag_mock1 = AGWrapper(N);
    AGWrapper ag_mock2 = AGWrapper(N);

    FloatWrapper float_mock1 = FloatWrapper(N);
    FloatWrapper float_mock2 = FloatWrapper(N);

    void SetUp(const ::benchmark::State &state)
    {

        for (uint64_t i = 0; i < 3 * N; i++)
        {
            coords1[i] = static_cast<float>(i);
            coords2[i] = static_cast<float>(i);
        }
    }

    void TearDown(const ::benchmark::State &state)
    {
    }

    void BM_calc_distance_array(benchmark::State &state)
    {
        for (auto _ : state)
        {
            _calc_distance_array(coords1, N, coords2, N, result);
        }
        state.SetItemsProcessed(N * N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
            N * state.iterations(),
            benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }

    void BM_DistanceArrayAg(benchmark::State &state)
    {
        for (auto _ : state)
        {
            DistanceArray(ag_mock1, ag_mock2, result);
        }
        state.SetItemsProcessed(N * N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
            N * state.iterations(),
            benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }


    void BM_DistanceArrayFloatWrap(benchmark::State &state)
    {
        for (auto _ : state)
        {
            DistanceArray(float_mock1, float_mock2, result);
        }
        state.SetItemsProcessed(N * N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
            N * state.iterations(),
            benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }

};

BENCHMARK_F(Distance, distance_array)
(benchmark::State &state)
{
    BM_calc_distance_array(state);
}


BENCHMARK_F(Distance, DistanceArrayAg)
(benchmark::State &state)
{
    BM_DistanceArrayAg(state);
}


BENCHMARK_F(Distance, DistanceArrayFloatWrap)
(benchmark::State &state)
{
    BM_DistanceArrayFloatWrap(state);
}

BENCHMARK_MAIN();