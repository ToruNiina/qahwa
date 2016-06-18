#ifndef QAHWA_WHAM_SOLVER
#define QAHWA_WHAM_SOLVER
#include "ProbabilityDensityFunction.hpp"
#include "UserDefinedFunctions.hpp"

namespace qahwa
{

class WhamSolver
{
  public:
    using window_type = std::pair<Trajectory, PerturbingPotential>;

  public:
    WhamSolver(const std::size_t bins, const double temperature)
        : bins_(bins), beta_(1.0 / (temperature * kB))
    {}
    ~WhamSolver() = default;

    /*!
     * @brief calculate free energy parameters
     * solve the WHAM equation self-consistently with iteration.
     * @return array of free energy parameters
     * @sa WhamSolver::reconstruct
     */
    std::vector<double> solve(const std::vector<window_type>& windows) const;

    /*!
     * @brief reconstruct unbiased probability density
     * using free_energy_parameter calculated by solve() method
     * @return unbiased probability density function
     * @sa WhamSolver::solve
     */
    ProbabilityDensityFunction
    reconstruct(const std::vector<double>& free_energy_parameter,
                const std::vector<window_type>& windows,
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

    double  beta() const {return this->beta_;}
    double& beta()       {return this->beta_;}

  private:

    /*!
     * @brief calculate denominator of wham-equation
     * calculate sum of (total_step * exp(-beta * W_k(R_ij)) * Fk) for k.
     * Fk == exp(beta * fk);
     */
    double calc_denominator(const std::vector<window_type>& windows,
                            const std::vector<double>& expfs,
                            const SnapShot& snapshot) const;

    bool is_converged(const std::vector<double>& f_prev,
                      const std::vector<double>& f) const;
    
  private:

    std::size_t bins_;
    double beta_;
};
   
}

#endif /* QAHWA_WHAM_SOLVER */
