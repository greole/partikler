#include "CGALTYPEDEFS.hpp"
#include <benchmark/benchmark.h>
#include "SearchCubes.cpp"
#include <vector>
#include <cstdlib>      // std::rand, std::srand


// template<class T>
// static void CreateSearchCubes(benchmark::State& state) {
//     size_t n_points =  state.range(0);
//     std::vector<Point> points = create_uniform_particle_cube(n_points);
//     float cube_dx = 1.0/((float) n_points);
//     for (auto _: state) {
//         T search (Logger(99), points, 3.0*cube_dx, false);
//     }
// }


// template<class T>
// static void CreateSearchCubesShuffled(benchmark::State& state) {
//   size_t n_points =  state.range(0);
//   std::vector<Point> points = create_uniform_particle_cube(n_points);
//   std::random_shuffle ( points.begin(), points.end() );

//   float cube_dx = 1.0/((float) n_points);
//   for (auto _: state) {
//     T search (Logger(99), points, 3.0*cube_dx, false);
//   }
// }

// template<class T>
// static void CreateSearchCubesRandomised(benchmark::State& state) {
//     size_t n_points =  state.range(0);
//     std::vector<Point> points = create_uniform_particle_cube(n_points);
//     float cube_dx = 1.0/((float) n_points);
//     disperse_particles(points, cube_dx*10.0);
//     for (auto _: state) {
//         T search (Logger(99), points, 3.0*cube_dx, false);
//     }
// }

static void OwnerCubeSearch(benchmark::State& state) {
  std::vector<Point> points = create_uniform_particle_cube(3);
  // float cube_dx = 1.0/((float) n_points);
  for (auto _: state) {
    SortedNeighbours ret {};
    ret.ownId.reserve(27*points.size()/2);
    ret.neighId.reserve(27*points.size()/2);
    owner_cube_search(points, 0, 27, 2.0, ret);
  }
}
static void OwnerCubeSearchVector(benchmark::State& state) {
  std::vector<Point> points = create_uniform_particle_cube(3);
  // float cube_dx = 1.0/((float) n_points);
  for (auto _: state) {
    SortedNeighbours ret {};
    ret.ownId.reserve(27*points.size()/2);
    ret.neighId.reserve(27*points.size()/2);
    vector_owner_cube_search(points, 0, 27, 2.0, ret);
  }
}

static void SortingParticles(benchmark::State& state) {
    size_t n_points =  state.range(0);
    std::vector<Point> points = create_uniform_particle_cube(n_points);
    float cube_dx = 1.0/((float) n_points);
    const SearchCubeDomain scd = initSearchCubeDomain(points, cube_dx);
    for (auto _: state) {
      countingSortParticles(scd, points);
    }
}

static void SortingSortedParticles(benchmark::State& state) {
  size_t n_points =  state.range(0);
  std::vector<Point> points = create_uniform_particle_cube(n_points);
  float cube_dx = 1.0/((float) n_points);
  const SearchCubeDomain scd = initSearchCubeDomain(points, cube_dx);
  SortedParticles sp = countingSortParticles(scd, points);
  for (auto _: state) {
    countingSortParticles(scd, sp.particles);
  }
}

static void MergedSort(benchmark::State& state) {
  size_t n_points =  state.range(0);
  std::vector<Point> points = create_uniform_particle_cube(n_points);
  float cube_dx = 1.0/((float) n_points);
  const SearchCubeDomain scd = initSearchCubeDomain(points, cube_dx);
  SortedParticles sp = countingSortParticles(scd, points);
  for (auto _: state) {
    mergedCountingSortAndNeighbourSearch(scd, sp.particles);
  }
}

static void SortAndCreateNeighbours(benchmark::State& state) {
  size_t n_points =  state.range(0);
  std::vector<Point> points = create_uniform_particle_cube(n_points);
  float cube_dx = 1.0/((float) n_points);
  const SearchCubeDomain scd = initSearchCubeDomain(points, cube_dx);
  for (auto _: state) {
    SortedParticles sp = countingSortParticles(scd, points);
    createNeighbours(scd, sp);
  }
}

static void CreateNeighbours(benchmark::State& state) {
  size_t n_points =  state.range(0);
  std::vector<Point> points = create_uniform_particle_cube(n_points);
  float cube_dx = 1.0/((float) n_points);
  const SearchCubeDomain scd = initSearchCubeDomain(points, cube_dx);
  SortedParticles sp = countingSortParticles(scd, points);
  for (auto _: state) {
    createNeighbours(scd, sp);
  }
}

BENCHMARK(OwnerCubeSearch)->MinTime(2);
BENCHMARK(OwnerCubeSearchVector)->MinTime(2);

// Particle sorting functions
BENCHMARK(SortingParticles)
  ->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);
BENCHMARK(SortingSortedParticles)
->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);

BENCHMARK(CreateNeighbours)
->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);
BENCHMARK(SortAndCreateNeighbours)
->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);
BENCHMARK(MergedSort)
->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);


//  Shuffle Particles
// BENCHMARK_TEMPLATE(CreateSearchCubesShuffled, SearchCubes)
// ->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);
// BENCHMARK_TEMPLATE(CreateSearchCubesShuffled, SearchCubes2)
// ->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);
// BENCHMARK_TEMPLATE(CreateSearchCubesShuffled, SearchCubes3)
// ->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);
// BENCHMARK_TEMPLATE(CreateSearchCubesShuffled, SearchCubes4)
// ->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);
// BENCHMARK_TEMPLATE(CreateSearchCubesShuffled, SearchCubes5)
// ->RangeMultiplier(2)->Range(2,2<<6)->Unit(benchmark::kMicrosecond);

//  Randomise Particles
// BENCHMARK_TEMPLATE(CreateSearchCubesRandomised, SearchCubes)
//     ->RangeMultiplier(2)->Range(2,2<<4)->Unit(benchmark::kMicrosecond);
// BENCHMARK_TEMPLATE(CreateSearchCubesRandomised, SearchCubes2)
//     ->RangeMultiplier(2)->Range(2,2<<4)->Unit(benchmark::kMicrosecond);
// BENCHMARK_TEMPLATE(CreateSearchCubesRandomised, SearchCubes5)
// ->RangeMultiplier(2)->Range(2,2<<4)->Unit(benchmark::kMicrosecond);


BENCHMARK_MAIN();
