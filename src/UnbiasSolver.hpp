#include "Histogram.hpp"
#include "ProbabilityDensityFunction.hpp"
#include "UserDefinedFunctions.hpp"

namespace qahwa
{

class UnbiasSolver
{
  public:
    using window_type = std::pair<Trajectory, PerturbingPotential>;

  public:
    UnbiasSolver() = default;
    UnbiasSolver(const std::size_t bins, const double temperature)
        : bins_(bins), beta_(1.0 / (kB * temperature))
    {}
    ~UnbiasSolver() = default;

    ProbabilityDensityFunction unbias(const Trajectory& traj,
                                      const PerturbingPotential perturb,
                                      const ReactionCoordinate& rctcrd) const;

    ProbabilityDensityFunction
    reconstruct(const std::vector<window_type>& windows,
                const ReactionCoordinate& rctcrd);

    Histogram<double>
    make_pmf(const ProbabilityDensityFunction& pdf) const;
    
  private:

    std::size_t bins_;
    double beta_;
};

}
