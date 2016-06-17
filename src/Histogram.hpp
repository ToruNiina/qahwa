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

    using size_type       = std::size_t;
    using index_type      = size_type;
    using value_type      = T_value;
    using container_type  = std::vector<T_value>;
    using coordinate_type = double;
    using range_type      = std::pair<coordinate_type, coordinate_type>;
    using iterator        = typename container_type::iterator;
    using const_iterator  = typename container_type::const_iterator;

  public:

    Histogram(){}
    Histogram(const std::size_t bin)
        : bins_(bin), values_(bin, 0)
    {}
    Histogram(const std::size_t bin, const double begin, const double end)
        : bins_(bin), range_(begin, end), values_(bin, 0)
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

    value_type at(coordinate_type coord) const
    {
        if(coord < this->range_.first)  return 0e0;
        if(coord > this->range_.second) return 0e0;

        const double dx  = (this->range_.second - this->range_.first) / bins_;
        const double pos = coord - this->range_.first;
        const std::size_t index = static_cast<std::size_t>(pos / dx);

        return values_.at(index);
    }

    value_type& at(coordinate_type coord)
    {
        if(coord < this->range_.first)
           throw std::out_of_range("histogram out of range");
        if(coord > this->range_.second)
           throw std::out_of_range("histogram out of range");

        const double dx  = (this->range_.second - this->range_.first) / bins_;
        const double pos = coord - this->range_.first;
        const std::size_t index = static_cast<std::size_t>(pos / dx);

        return values_.at(index);
    }

    range_type   range() const {return range_;}
    range_type&  range()       {return range_;}

    coordinate_type   range_begin() const {return range_.first;}
    coordinate_type&  range_begin()       {return range_.first;}
    coordinate_type   range_end()   const {return range_.second;}
    coordinate_type&  range_end()         {return range_.second;}

    std::size_t  bins()  const {return bins_;}
    std::size_t& bins()        {return bins_;}
    const container_type& values() const {return values_;}
          container_type& values()       {return values_;}

  private:

    std::size_t bins_  = 0;
    range_type  range_ = {0.0, 0.0};
    container_type values_;
};

template<typename T_val>
std::ostream& operator<<(std::ostream& os, const Histogram<T_val>& hist)
{
    os << "reaction-coord  value" << std::endl;
    const double beg = hist.range().first;
    const double end = hist.range().second;
    const double dx  = (end - beg) / hist.bins();

    std::size_t index = 0;
    for(auto iter = hist.cbegin(); iter != hist.cend(); ++iter)
    {
        os << beg + dx * index << " " << *iter << std::endl;
        ++index;
    }

    return os;
}

}

#endif /* QAHWA_HISTOGRAM */
