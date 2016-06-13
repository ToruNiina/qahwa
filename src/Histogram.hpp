#ifndef QAHWA_HISTOGRAM
#define QAHWA_HISTOGRAM
#include "Definitions.hpp"
#include <iostream>

namespace qahwa
{

template<typename T_value>
class Histogram
{
  public:

    using size_type      = std::size_t;
    using index_type     = size_type;
    using value_type     = T_value;
    using container_type = std::vector<T_value>;
    using iterator       = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

  public:
    Histogram(){}
    Histogram(const std::size_t bin)
        : bins_(bin), values_(bin)
    {}
    Histogram(const std::size_t bin, const double begin, const double end)
        : bins_(bin), range_begin_(begin), range_end_(end), values_(bin)
    {}
    ~Histogram(){}

    void        clear() {values_.clear();}
    bool        empty() const {return values_.empty();}
    std::size_t size()  const {return values_.size();}
    value_type  at(index_type i) const {return values_.at(i);}
    value_type& at(index_type i)       {return values_.at(i);}
    value_type  operator[](index_type i) const {return values_[i];}
    value_type& operator[](index_type i)       {return values_[i];}

    iterator begin() {return values_.begin();}
    iterator end()   {return values_.end();}
    const_iterator cbegin() const {return values_.cbegin();}
    const_iterator cend()   const {return values_.cend();}

    value_type  at(double coord) const
    {
        if(coord < this->begin_) return 0e0;
        if(coord > this->end_) return 0e0;

        const double dx  = (this->end_ - this->begin_) / bins_;
        const double pos = coord - this->begin_;
        const std::size_t index = static_cast<std::size_t>(pos / dx);

        return values_.at(index);
    }

    value_type& at(double coord)
    {
        if(coord < this->begin_)
           throw std::out_of_range("histgram out of range");
        if(coord > this->end_)
           throw std::out_of_range("histgram out of range");

        const double dx  = (this->end_ - this->begin_) / bins_;
        const double pos = coord - this->begin_;
        const std::size_t index = static_cast<std::size_t>(pos / dx);

        return values_.at(index);
    }

    double       range_begin() const {return range_begin_;}
    double&      range_begin()       {return range_begin_;}
    double       range_end()   const {return range_end_;}
    double&      range_end()         {return range_end_;}
    std::size_t  bins()  const {return bins_;}
    std::size_t& bins()        {return bins_;}
    const container_type& values() const {return values_;}
          container_type& values()       {return values_;}

  private:

    std::size_t bins_  = 0;
    double      range_begin_ = 0e0;
    double      range_end_   = 0e0;
    container_type values_;
};

template<typename T_val>
std::ostream& operator<<(std::ostream& os, const Histogram<T_val>& hist)
{
    os << "reaction-coord  value" << std::endl;
    const double beg = hist.range_begin();
    const double end = hist.range_end();
    const double dx  = (end - beg) / hist.bins();

    std::size_t index = 0;
    for(auto iter = hist.cbegin(); iter != hist.cend(); ++iter)
    {
        os << beg + dx * index << " " << *iter << std::endl;
    }

    return os;
}

}

#endif /* QAHWA_HISTOGRAM */
