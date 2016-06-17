#ifndef QAHWA_DEFINITIONS
#define QAHWA_DEFINITIONS
#include "math/LinearAlgebra.hpp"
#include "DCDReader.hpp"
#include "DCDWriter.hpp"
#include "InputFileReader.hpp"

namespace qahwa
{
    using Vector3d   = ax::Vector3d;
    using SnapShot   = std::vector<Vector3d>;
    using Trajectory = std::vector<SnapShot>;

    using DCDReader = coffeemill::DCDReader;
    using DCDWriter = coffeemill::DCDWriter;
    using InputFileReader = coffeemill::InputFileReader;

    constexpr static double kB = 1.986231313e-3; // kcal/mol/K
    constexpr static double tolerance = 1e-8;
    struct Settings
    {
        double      temperature;
        std::size_t bins;
    };
}

#endif /* QAHWA_DEFINITIONS */
