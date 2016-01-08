#ifndef H_HISTO_DEF
#define H_HISTO_DEF

#include <vector>
#include <string>

#include "named_func.hpp"

class HistoDef{
public:
  HistoDef(const std::vector<double> &bins,
           const NamedFunc &var,
           const std::string &x_title = "",
           const NamedFunc &cut = NamedFunc(1.),
           const NamedFunc &weight = FUNC(b.weight()));
  HistoDef(size_t nbins,
           double xmin,
           double xmax,
           const NamedFunc &var,
           const std::string &x_title = "",
           const NamedFunc &cut = NamedFunc(1.),
           const NamedFunc &weight = FUNC(b.weight()));
  HistoDef(const HistoDef &) = default;
  HistoDef& operator=(const HistoDef &) = default;
  HistoDef(HistoDef &&) = default;
  HistoDef& operator=(HistoDef &&) = default;
  ~HistoDef() = default;

  size_t GetNbins() const;
  std::string GetName() const;
  
  std::vector<double> bins_;
  NamedFunc var_, cut_, weight_;
  std::string x_title_;
private:
  static std::vector<double> GetEdges(size_t nbins, double xmin, double xmax);
};

#endif
