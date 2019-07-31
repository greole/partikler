#ifndef SEARCHCUBES_H
#define SEARCHCUBES_H

#include "SPHDatastructures.hpp"

SortedNeighbours
createSTLNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Point> & pos,
    std::vector<SearchCube> & searchCubes,
    const std::vector<Facet_handle> & facets
    );

SortedNeighbours
createNeighbours(
    const SearchCubeDomain scd,
    const std::vector<Point> & pos,
    std::vector<SearchCube> & searchCubes
    );

SortedParticles countingSortParticles(
    const SearchCubeDomain scd, const std::vector<Point> &unsorted_particles);

#endif
