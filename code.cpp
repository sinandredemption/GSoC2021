#include <iostream>
#include <cmath>
#include <limits>
#include <boost/multiprecision/cpp_bin_float.hpp>

template <typename T>
T f1(T x) {
  return x * x - 1;
}

template <typename T>
T derivative(T(*f)(T), const T& x) {
  const T dx = sqrt(std::numeric_limits<T>::epsilon()) * x;
  return ((*f)(x + dx) - (*f)(x)) / dx;
}

template <typename T>
T newtons_iteration(T(*f)(T), const T& x) {
  return x - f(x) / derivative(f, x);
}

int main() {
  std::cout.precision(100);
  std::cout << std::fixed;

  std::cout << "DEMO #1: Solving x^2 == 1 using Newton's iteration\n"
    << "Using initial approximation x = 1.1...\n" << std::endl;

  boost::multiprecision::cpp_bin_float_100 x("1.1"), x_prime = 0;

  for (int i = 0; x != x_prime; ++i) {
    std::cout << "Iteration #" <<  i << ": x = " << x;

    x_prime = x;
    x = newtons_iteration(f1, x);

    std::cout << " (press ENTER to continue)";
    std::cin.get();
  }

  std::cout << "\nDEMO #2: Gauss-Lengendre algorithm for Pi\n"
    << "see < https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_algorithm >" << std::endl;

  using mp_float = boost::multiprecision::cpp_bin_float_100;

  mp_float a(1), b = 1 / sqrt(mp_float(2)), t("0.25"), pi = 0;
  int correct_digits = 0;
  for (int i = 1; correct_digits < 99; ++i) {
    std::cout << "Iteration #" << i << ": ";

    mp_float a_prime = (a + b) / 2;
    mp_float b_prime = sqrt(a * b);
    mp_float t_prime = t - exp2(i - 1) * (a - a_prime) * (a - a_prime);

    a = a_prime;
    b = b_prime;
    t = t_prime;

    pi = (a + b) * (a + b) / (4 * t);

    correct_digits = (int)-log10(abs(pi - boost::math::constants::pi<mp_float>()));

    std::cout << "pi = " << pi;
    std::cout << " | " << correct_digits << " digits correct | ";
    std::cout << " (press ENTER to continue)";
    std::cin.get();
  }

  return 0;
}
