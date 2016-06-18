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
        : bins_(bins), beta_(1.0 / (kB * temperature))
    {}
    ~UnbiasSolver() = default;

    ProbabilityDensityFunction unbias(const Trajectory& traj,
                                      const PerturbingPotential perturb,
                                      const ReactionCoordinate& rctcrd) const;

    Histogram<double>
    make_pmf(const ProbabilityDensityFunction& pdf) const
    {
        Histogram<double> pmf(pdf.bins(), pdf.range_begin(), pdf.range_end());
        auto potential = pmf.begin();
        for(auto iter = pdf.cbegin(); iter != pdf.cend(); ++iter)
        {
            *potential = -1.0 * std::log(*iter) / this->beta_;
            ++potential;
        }
        return pmf;
    }

  private:

    std::size_t bins_;
    double beta_;
};

}
