#include "UnbiasSolver.hpp"

namespace qahwa
{

ProbabilityDensityFunction UnbiasSolver::unbias(const Trajectory& traj,
                                      const PerturbingPotential perturb,
                                      const ReactionCoordinate& rctcrd) const 
{
    const double beta = 1.0 / (this->temperature_ * kB);
    HistogramMaker histmaker(this->bins_);
    const auto range = histmaker.find_range(traj, rctcrd);
    histmaker.set_range(range);

    ProbabilityDensityFunction unbiased = 
        histmaker.make_prob_dens_func(traj, rctcrd);

    const double min_p = 1.0 / traj.size();
    for(auto iter = traj.cbegin(); iter != traj.cend(); ++iter)
    {
        const double x = rctcrd(*iter);
        const double W = perturb(*iter);
        unbiased.at(x) += min_p * std::exp(beta * W);
    }

    return unbiased;
}


}
