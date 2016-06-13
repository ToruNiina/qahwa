#ifndef QAHWA_HISTOGRAM_MAKER
#define QAHWA_HISTOGRAM_MAKER
#include "Histogram.hpp"
#include "UserDefinedFunctions.hpp"

namespace qahwa
{

// one-dimentional only...
class HistogramMaker
{
  public:
    using histogram_type = Histogram<std::size_t>;

  public:
    HistogramMaker(){}
    HistogramMaker(const std::size_t bins): bins_(bins){}
    HistogramMaker(const std::size_t bins, const double begin, const double end)
        : bins_(bins), begin_(begin), end_(end)
    {}
    ~HistogramMaker() = default;

    std::pair<double, double> find_range(const Trajectory& traj,
        const ReactionCoordinate<double>& reaction_coord) const;

    histogram_type make_histogram(const Trajectory& traj,
        const ReactionCoordinate<double>& reaction_coord) const;

    void set_range(const std::pair<double, double>& range)
    {
        this->begin_ = range.first;
        this->end_ = range.second;
        return;
    }

    std::size_t  bin() const {return bins_;}
    std::size_t& bin()       {return bins_;}
    double  range_begin() const {return begin_;}
    double& range_begin()       {return begin_;}
    double  range_end() const {return end_;}
    double& range_end()       {return end_;}

  private:
    std::size_t bins_  = 100;
    double      begin_ = 0e0;
    double      end_   = 0e0;
};



}

#endif /* QAHWA_HISTOGRAM_MAKER */
