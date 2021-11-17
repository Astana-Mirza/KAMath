#ifndef EXTRA_FIELDS_H
#define EXTRA_FIELDS_H

#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <complex>

using boost::multiprecision::cpp_int;
using Rational = boost::rational<cpp_int>;
using Complex = std::complex<long double>;


// функции проверки на ноль соответствующего поля
bool is_zero(Rational r) { return r == 0; }
bool is_zero(Complex c)  { return c.real()==0 && c.imag()==0; }


// строковое представление рационального числа
std::string to_string(const Rational& r) {
	if (r.denominator() != 1)
		return boost::lexical_cast<std::string>(r);
	else return boost::lexical_cast<std::string>(r.numerator());
}


// строковое представление комплексного числа
std::string to_string(const Complex& c) {
	if (c.real() == 0 && c.imag() == 0) return "0";
	if (c.real() == 0)
		return boost::lexical_cast<std::string>(c.imag()) + "i";
	if (c.imag() == 0)
		return boost::lexical_cast<std::string>(c.real());
	if (c.imag() > 0)
		return "(" +  boost::lexical_cast<std::string>(c.real()) + "+" +  boost::lexical_cast<std::string>(c.imag()) + "i)";
	return "(" +  boost::lexical_cast<std::string>(c.real()) + boost::lexical_cast<std::string>(c.imag()) + "i)";
}

#endif	// #ifndef EXTRA_FIELDS_H

