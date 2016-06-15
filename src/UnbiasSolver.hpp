#include "Histogram.hpp"
#include "ProbabilityDensityFunction.hpp"
#include "UserDefinedFunctions.hpp"

namespace qahwa
{

class UnbiasSolver
{
  public:
    UnbiasSolver() = default;
    UnbiasSolver(const std::size_t bins, const double temperature)
        : bins_(bins), temperature_(temperature)
    {}
    ~UnbiasSolver() = default;

    ProbabilityDensityFunction unbias(const Trajectory& traj,
                                      const PerturbingPotential perturb,
                                      const ReactionCoordinate& rctcrd) const;

  private:
    std::size_t bins_;
    double temperature_;
};

}
