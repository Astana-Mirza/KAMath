#include "field_reader.h"


namespace KAMParser {


/**************************************************************************
** Считать элемент поля Галуа характеристики 2
**************************************************************************/
KAMath::Field_FGF2 read_from_stream(std::istream* ip, KAMath::Field_FGF2 one) {
	char ch;
	size_t sz;
	ip->get(ch);
	if (isdigit(ch))
	{
		ip->putback(ch);
		*ip >> sz;
		if (sz > 1)
			throw std::runtime_error{"Wrong Galois field value"};
		return KAMath::Field_FGF2(one.prim_pol(), sz);
	}
	else if (ch != PRIM_CHAR)
		throw std::runtime_error{"Wrong Galois field value"};
	return KAMath::Field_FGF2(one.prim_pol(),
		KAMath::Field_FGF2::to_signature(1, one.prim_pol()));
}


/**************************************************************************
** Считать элемент поля рациональных чисел
**************************************************************************/
KAMath::Rational read_from_stream(std::istream* ip, KAMath::Rational one) {
	char ch;
	ip->get(ch);
	cpp_int num = isdigit(ch) ? ch-'0' : 0;
	cpp_int denom = 1;
	do { ip->get(ch); } while (isspace(ch) && ch != '\n');
	while (isdigit(ch)) {
		num = num * 10 + (ch-'0');
		ip->get(ch);
	}
	while (isspace(ch) && ch != '\n') { ip->get(ch); }
	if (ch == '/') {
		do { ip->get(ch); } while (isspace(ch) && ch != '\n');
		if (!isdigit(ch) || ch == '0')
			throw std::runtime_error{"Wrong rational value"};
		denom = 0;
		while (isdigit(ch)) {
			denom = denom * 10 + (ch-'0');
			ip->get(ch);
		}
	}
	ip->putback(ch);
	return KAMath::Rational(num, denom);
}


/**************************************************************************
** Считать элемент поля комплексных чисел
**************************************************************************/
KAMath::Complex read_from_stream(std::istream* ip, KAMath::Complex one) {
	char ch;
	long double ld;
	ip->get(ch);
	if (ch == PRIM_CHAR)
		throw std::runtime_error{"Wrong complex value"};
	if (ch == IMAGINARY_UNIT)
		return KAMath::Complex(0, 1);
	ip->putback(ch);
	*ip >> ld;
	ip->get(ch);
	if (ch == IMAGINARY_UNIT)
		return KAMath::Complex(0, ld);
	ip->putback(ch);
	return KAMath::Complex(ld, 0);
}


/**************************************************************************
** Считать элемент поля вычетов по модулю простого числа
**************************************************************************/
KAMath::Field_ZP read_from_stream(std::istream* ip, KAMath::Field_ZP one) {
	int64_t val;
	*ip >> val;
	return KAMath::Field_ZP(one.charact(), val);
}


};

