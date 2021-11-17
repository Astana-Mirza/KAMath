#ifndef MODULO_PRIME_H
#define MODULO_PRIME_H

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace KAMath {

class Field_ZP {
public:
	Field_ZP(const size_t P=2, const int64_t v=0, const bool prime=false):
			ch{P}, val{(v > 0) ? v%P : P-((-v)%P)}
		{ if (!prime) check_prime(); }
	Field_ZP(const Field_ZP& v) : ch{v.ch}, val{v.val} {}
	void operator = (const Field_ZP& v) { val = v.val % ch; }
	void operator = (const size_t v) { val = v % ch; }

	Field_ZP inv() const {	// обратный элемент по умножению
		if (val == 0)
			throw std::runtime_error{
				"Attempt to find inverse of 0"};
		return Field_ZP(ch, std::pow(val, ch-2));
	}
	Field_ZP operator - () const
		{ return Field_ZP(ch, ch-val, true); }
	Field_ZP operator + (const Field_ZP& v) const
		{ return Field_ZP(ch, val+v.val, true); }
	Field_ZP operator - (const Field_ZP& v) const
		{ return Field_ZP(ch, val+(-v).val, true); }
	Field_ZP operator * (const Field_ZP& v) const
		{ return Field_ZP(ch, val*v.val, true); }
	Field_ZP operator / (const Field_ZP& v) const
		{ return Field_ZP(ch, val*v.inv().val, true); }
	
	Field_ZP operator += (const Field_ZP& v) { *this = *this+v; return *this; }
	Field_ZP operator -= (const Field_ZP& v) { *this = *this-v; return *this; }
	Field_ZP operator *= (const Field_ZP& v) { *this = *this*v; return *this; }
	Field_ZP operator /= (const Field_ZP& v) { *this = *this/v; return *this; }
	bool operator == (const Field_ZP& v) const { return ch==v.ch && val==v.val; }

	size_t charact() const { return ch; }
	size_t get_val() const { return val; }
private:
	size_t ch;	// характеристика поля	
	size_t val;	// значение элемента
	void check_prime() const;
};


Field_ZP pow(const Field_ZP& v, const int64_t p);
bool is_zero(const Field_ZP& val);
std::string to_string(const Field_ZP& val);
};

#endif	// #ifndef MODULO_PRIME_H
