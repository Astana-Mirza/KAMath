#include "extra_fields.h"


namespace KAMath {


/**************************************************************************
** Функции проверки на ноль
**************************************************************************/
bool is_zero(Rational r) { return r == 0; }
bool is_zero(Complex c)  { return c.real()==0 && c.imag()==0; }


/**************************************************************************
** Строковое представление рационального числа
**************************************************************************/
std::string to_string(const Rational& r) {
	if (r.denominator() != 1)
		return boost::lexical_cast<std::string>(r);
	else return boost::lexical_cast<std::string>(r.numerator());
}


/**************************************************************************
** Строковое представление комплексного числа
**************************************************************************/
std::string to_string(const Complex& c) {
	std::ostringstream ss1, ss2;
	long double r = c.real(), i = c.imag();
	// округление до 0
	if (c.real() < 0.0001 && c.real() > -0.0001) r = 0;
	if (c.imag() < 0.0001 && c.imag() > -0.0001) i = 0;
	if (r == 0 && i == 0) return "0";
	if (r == 0) {
		ss1 << i;
		std::string im = ss1.str();
		if (im == "1")  return  "i";
		if (im == "-1") return "-i";
		return im + "i";
	}
	if (i == 0) {
		ss1 << r;
		return ss1.str();
	}
	if (i > 0) {
		ss1 << r;
		ss2 << i;
		std::string im = ss2.str();
		return "(" +  ss1.str() + "+" + ((im=="1") ? "" : im)  + "i)";
	}
	ss1 << r;
	ss2 << i;
	std::string im = ss2.str();
	return "(" +  ss1.str() + ((im=="-1") ? "-" : im) + "i)";
}


/**************************************************************************
** Возведение в степень рационального числа
**************************************************************************/
Rational pow(const Rational& v, const int64_t p) {
	return Rational(boost::multiprecision::pow(v.numerator(), (size_t)p),
		boost::multiprecision::pow(v.denominator(), (size_t)p));
}


/**************************************************************************
** Возведение в степень комплексного числа
**************************************************************************/
Complex pow(const Complex& v, const int64_t p) {
	return std::pow(v, p);
}


};

