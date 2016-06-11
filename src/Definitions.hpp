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
}

#endif /* QAHWA_DEFINITIONS */
