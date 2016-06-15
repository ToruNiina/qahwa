#include "UserDefinedFunctions.hpp"

namespace qahwa
{

// this is for example. anchor.
double PerturbingPotential::operator()(const SnapShot& snapshot) const
{
    const double dist = length(snapshot.at(1));
    const double d_nat = length(this->anchored_position_);

    return this->coefficient_ * (dist - d_nat) * (dist - d_nat);
}

// this is for example. distance from origin.
double ReactionCoordinate::operator()(const SnapShot& snapshot) const
{
    return length(snapshot.at(1));
}


}
