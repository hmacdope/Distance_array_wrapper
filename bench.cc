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

    void SetUp(const ::benchmark::State &state)
    {

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
        state.SetItemsProcessed(N*N * state.iterations());
        state.counters["Per Result"] = benchmark::Counter(
        N * state.iterations(),
        benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
    }
};

BENCHMARK_F(Distance, calc_bonds)
(benchmark::State &state)
{
    BM_calc_distance_array(state);
}

BENCHMARK_MAIN();