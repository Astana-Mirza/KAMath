#ifndef EXTRA_FIELDS_H
#define EXTRA_FIELDS_H

#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <complex>
#include <sstream>


namespace KAMath {

using boost::multiprecision::cpp_int;
using Rational = boost::rational<cpp_int>;
using Complex = std::complex<long double>;

// функции проверки на ноль соответствующего поля
bool is_zero(Rational r);
bool is_zero(Complex c);

// строковое представление чисел
std::string to_string(const Rational& r);
std::string to_string(const Complex& c);

// возведение в степень
Rational pow(const Rational& v, const int64_t p);
Complex pow(const Complex& v, const int64_t p);

};

#endif	// #ifndef EXTRA_FIELDS_H

