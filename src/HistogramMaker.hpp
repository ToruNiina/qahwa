#ifndef QAHWA_HISTOGRAM_MAKER
#define QAHWA_HISTOGRAM_MAKER
#include "Histogram.hpp"
#include "UserDefinedFunctions.hpp"

namespace qahwa
{

class HistogramMaker
{
  public:
    using histogram_type = Histogram<std::size_t>;
    using range_type     = histogram_type::range_type;

  public:
    HistogramMaker(){}
    HistogramMaker(const std::size_t bins): bins_(bins){}
    HistogramMaker(const std::size_t bins, const double begin, const double end)
        : bins_(bins), range_(begin, end)
    {}
    ~HistogramMaker() = default;

    range_type find_range(const Trajectory& traj,
                          const ReactionCoordinate& reaction_coord) const;

    histogram_type make_histogram(const Trajectory& traj,
                       const ReactionCoordinate& reaction_coord) const;

    Histogram<double> make_prob_dens_func(const Trajectory& traj,
                       const ReactionCoordinate& reaction_coord) const;

    void set_range(const range_type& range){range_ = range;}

    std::size_t  bin()    const {return bins_;}
    std::size_t& bin()          {return bins_;}
    double  range_begin() const {return range_.first;}
    double& range_begin()       {return range_.first;}
    double  range_end()   const {return range_.second;}
    double& range_end()         {return range_.second;}

  private:
    std::size_t bins_  = 100;
    range_type  range_ = {0.0, 0.0};
};



}

#endif /* QAHWA_HISTOGRAM_MAKER */
