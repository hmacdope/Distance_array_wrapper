#include <benchmark/benchmark.h>
#include "distance_array_comparison.h"
#include "distance_array_batched.h"

class Distance : public benchmark::Fixture
{
public:
    // setup
    constexpr static uint64_t N = 10000;

    float coords1[3 * N];
    float coords2[3 * N];

    double result[N * N] = {0};

    // classes that mock atomgroup
    AGWrapper ag_mock1 = AGWrapper(N);
    AGWrapper ag_mock2 = AGWrapper(N);
    AGWrapper ag_mock_non_contig1 = AGWrapper(N, false);
    AGWrapper ag_mock_non_contig2 = AGWrapper(N, false);

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

    void BM_DistanceArrayAgNonContig(benchmark::State &state)
    {
        for (auto _ : state)
        {
            DistanceArray(ag_mock_non_contig1, ag_mock_non_contig2, result);
        }
        state.SetItemsProcessed(N * N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
            N * state.iterations(),
            benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }

    // void BM_DistanceArrayFloatWrap(benchmark::State &state)
    // {
    //     for (auto _ : state)
    //     {
    //         DistanceArray(float_mock1, float_mock2, result);
    //     }
    //     state.SetItemsProcessed(N * N * state.iterations());
    //     state.counters["Per Result"] = benchmark::Counter(
    //         N * state.iterations(),
    //         benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    // }

    void BM_DistanceArrayMixed(benchmark::State &state)
    {
        for (auto _ : state)
        {
            DistanceArray(float_mock1, ag_mock1, result);
        }
        state.SetItemsProcessed(N * N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
            N * state.iterations(),
            benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }

    void BM_DistanceArrayAgPreload(benchmark::State &state)
    {
        for (auto _ : state)
        {
            DistanceArrayPreloadAG(ag_mock1, ag_mock2, result);
        }
        state.SetItemsProcessed(N * N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
            N * state.iterations(),
            benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }

    void BM_AgPreload(benchmark::State &state)
    {
        for (auto _ : state)
        {
            ag_mock1.preload();
        }
        state.SetItemsProcessed(N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
            N * state.iterations(),
            benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }

    void BM_DistanceArrayAgExternalPreload(benchmark::State &state)
    {
        for (auto _ : state)
        {
            DistanceArrayBatched(ag_mock1, ag_mock2, result, 64);
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

BENCHMARK_F(Distance, DistanceArrayAgNonContig)
(benchmark::State &state)
{
    BM_DistanceArrayAgNonContig(state);
}

// BENCHMARK_F(Distance, DistanceArrayFloatWrap)
// (benchmark::State &state)
// {
//     BM_DistanceArrayFloatWrap(state);
// }

BENCHMARK_F(Distance, DistanceArrayMixed)
(benchmark::State &state)
{
    BM_DistanceArrayMixed(state);
}

BENCHMARK_F(Distance, DistanceArrayAgPreload)
(benchmark::State &state)
{
    BM_DistanceArrayAgPreload(state);
}

BENCHMARK_F(Distance, AgPreload)
(benchmark::State &state)
{
    BM_AgPreload(state);
}


BENCHMARK_F(Distance, DistanceArrayAgExternalPreload)
(benchmark::State &state)
{
    BM_DistanceArrayAgExternalPreload(state);
}

BENCHMARK_MAIN();