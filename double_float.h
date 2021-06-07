#ifndef INC_DOUBLE_FLOAT_H
#define INC_DOUBLE_FLOAT_H

#include <string>
#include <utility>
#include <cassert>
#include <limits>
#include <sstream>

template <typename FloatingPointType>
class double_float {
public:
  using float_type = FloatingPointType;

  // Constructors
  double_float() {}
  double_float(const float_type& a_, const float_type& b_) : a(a_), b(b_) {}
  double_float(const std::pair<float_type, float_type>& p) : a(p.first), b(p.second) {}

  // Methods
  //std::string get_str() const;
  //void set_str(const std::string& str);
  std::string to_string(int precision) const;
  double_float<float_type> negative() const { return double_float<float_type>(std::make_pair(-a, -b)); }

  // Getters/Setters
  //float_type& first() { return a; }
  //float_type& second() { return b; }
  const float_type& first() const { return a; }
  const float_type& second() const { return b; }

  // Ops
  //double_float operator+(const float_t& y);
  //double_float operator*(const float_t& y);

  // Helper functions
  //static double_float<float_type> add(const double_float<FloatingPointType>& a, const double_float<float_type>& b);
  static std::pair<float_type, float_type> fast_exact_sum(const float_type& a, const float_type& b);
  static std::pair<float_type, float_type> exact_sum(const float_type& a, const float_type& b);
  static std::pair<float_type, float_type> exact_product(const float_type& a, const float_type& b);

  static void normalize_pair(std::pair<float_type, float_type>* p, bool fast = true);
  static std::pair<float_type, float_type> split(const float_type& a);

  // -- DEBUGGING
  std::string to_string_raw() const;
  // --
private:
  float_type a, b;
};

// Exact edition of two floating point numbers, given |a| > |b|
template<typename FloatingPointType>
std::pair<FloatingPointType, FloatingPointType>
double_float<FloatingPointType>::fast_exact_sum(const float_type& a, const float_type& b)
{
  assert(std::abs(a) >= std::abs(b));

  std::pair<float_type, float_type> out;
  out.first = a + b;
  out.second = b - (out.first - a);

  return out;
}

template<typename FloatingPointType>
std::pair<FloatingPointType, FloatingPointType>
double_float<FloatingPointType>::exact_sum(const float_type& a, const float_type& b)
{
  std::pair<float_type, float_type> out;

  out.first = a + b;
  float_type v = out.first - a;
  out.second = (a - (out.first - v)) + (b - v);

  return out;
}
template<typename FloatingPointType>
inline void
double_float<FloatingPointType>::normalize_pair(std::pair<float_type, float_type>* p, bool fast)
{
  *p = (fast ? fast_exact_sum(p->first, p->second) : exact_sum(p->first, p->second));
}

// TODO Merge into one function using std::numeric_limits()
template<typename FloatingPointType>
std::pair<FloatingPointType, FloatingPointType>
inline double_float<FloatingPointType>::split(const FloatingPointType& a)
{
  static_assert(std::numeric_limits<FloatingPointType>::is_iec559,
    "double_float<> invoked with non-native floating-point unit");

  constexpr int MantissaBits = sizeof(FloatingPointType) == 4 ? 23
    : sizeof(FloatingPointType) == 8 ? 53 : 0;
  constexpr int SplitBits = MantissaBits / 2 + 1;
  constexpr FloatingPointType SplitThreshold =
    std::numeric_limits<FloatingPointType>::max() / FloatingPointType(1 << (SplitBits + 1));
  constexpr FloatingPointType Splitter = FloatingPointType((1 << SplitBits) + 1);

  FloatingPointType temp, a_ = a, hi, lo;
  std::pair<FloatingPointType, FloatingPointType> out;

  if (a > SplitThreshold || a < -SplitThreshold) {
    a_ *= 1. / FloatingPointType(1 << (SplitBits + 1));
    temp = Splitter * a;
    hi = temp - (temp - a_);
    lo = a_ - hi;
    hi *= FloatingPointType(1 << (SplitBits + 1));
    lo *= FloatingPointType(1 << (SplitBits + 1));
  }
  else {
    temp = Splitter * a;
    hi = temp - (temp - a);
    lo = a - hi;
  }

  return std::make_pair(hi, lo);
}

template<typename FloatingPointType>
std::pair<FloatingPointType, FloatingPointType>
double_float<FloatingPointType>::exact_product(const float_type& a, const float_type& b)
{
  auto a_split = split(a);
  auto b_split = split(b);
  std::pair<float_type, float_type> p;

  p.first = a * b;
  p.second = ((a_split.first * b_split.first - p.first)
    + a_split.first * b_split.second + a_split.second * b_split.first)
    + a_split.second * b_split.second;

  return p;
}

template <typename FloatingPointType>
inline double_float<FloatingPointType>
operator+(const double_float<FloatingPointType>& a, const FloatingPointType& b)
{
  using double_float_t = double_float<FloatingPointType>;

  auto s = double_float_t::exact_sum(a.first(), b);

  s.second += a.second();
  double_float_t::normalize_pair(&s);

  return double_float_t(s);
}

template <typename FloatingPointType>
inline double_float<FloatingPointType>
operator-(const double_float<FloatingPointType>& a, const FloatingPointType& b)
{
  return double_float<FloatingPointType>(a + (-b));
}

template <typename FloatingPointType>
inline double_float<FloatingPointType>
operator+(const double_float<FloatingPointType>& a, const double_float<FloatingPointType>& b)
{
  using double_float_t = double_float<FloatingPointType>;

  std::pair<FloatingPointType, FloatingPointType> s, t;

  s = double_float_t::exact_sum(a.first(), b.first());
  t = double_float_t::exact_sum(a.second(), b.second());

  s.second += t.first;
  double_float_t::normalize_pair(&s);
  s.second += t.second;
  double_float_t::normalize_pair(&s);

  return double_float_t(s);
}

template <typename FloatingPointType>
inline double_float<FloatingPointType>
operator-(const double_float<FloatingPointType>& a, const double_float<FloatingPointType>& b)
{
  return double_float<FloatingPointType>(a + b.negative());
}

template <typename FloatingPointType>
inline double_float<FloatingPointType>
operator*(const double_float<FloatingPointType>& a, const FloatingPointType& b)
{
  using double_float_t = double_float<FloatingPointType>;

  auto p = double_float_t::exact_product(a.first(), b);
  p.second += a.second() * b;

  double_float_t::normalize_pair(&p);

  return double_float_t(p);
}


template <typename FloatingPointType>
inline double_float<FloatingPointType>
operator*(const double_float<FloatingPointType>& a, const double_float<FloatingPointType>& b)
{
  using double_float_t = double_float<FloatingPointType>;

  auto p = double_float_t::exact_product(a.first(), b.first());
  p.second += a.first() * b.second() + a.second() * b.first();

  double_float_t::normalize_pair(&p);

  return double_float_t(p);
}

template<typename FloatingPointType>
inline std::string double_float<FloatingPointType>::to_string(int precision) const
{
  std::ostringstream ss;

  double_float<FloatingPointType> d(*this);

  while (precision-- >= 0 || d.first() > 1.) {
    using namespace std;

    if (d.first() > 1.) {
      int d_exp = (int)log10(d.first());
      int digit = floor(d.first() / pow(10, d_exp));

      ss << digit;
      d = d - digit * pow(10, d_exp);
    }
    else {
      if (ss.str().find('.') == std::string::npos)
        ss << '.';

      d = d * 10.;
      int digit = int(d.first());
      ss << digit;
      d = d - (double)digit;
    }
  }

  return ss.str();
}

template<typename FloatingPointType>
inline std::string double_float<FloatingPointType>::to_string_raw() const
{
  std::stringstream ss;
  ss << "a:\t" << std::hexfloat << a << "\n";
  ss << "b:\t" << std::hexfloat << b << "\n";
  return ss.str();
}
#endif // ! INC_DOUBLE_FLOAT_H
