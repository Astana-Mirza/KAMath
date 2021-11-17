#ifndef MONOM_H
#define MONOM_H

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <map>
#include <set>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/lexical_cast.hpp>

namespace KAMath {

using boost::multiprecision::cpp_int;

struct Monom {
	mutable std::map<std::string, cpp_int> vars;

	Monom(const std::map<std::string, cpp_int>& v={}) : vars{v} {}
	Monom(const std::vector<std::string>& v, const std::vector<cpp_int>& deg);
	Monom(const Monom& m) : vars{m.vars} {}
	void  operator = (const Monom& m) { vars = m.vars; }

	cpp_int full_deg() const;		// полная степень монома
	void remove_var(const std::string& var);
	Monom operator *  (const Monom& m) const;
	Monom operator /  (const Monom& m) const;
	bool  operator == (const Monom& m) const;
	bool  operator != (const Monom& m) const { return !(*this == m); }
	Monom operator *= (const Monom& m) { *this = *this*m; return *this; }

	bool is_divisible(const Monom& m) const;	// делится ли на моном m
	size_t var_n() const { return vars.size(); }	// количество переменных

	// НОК мономов a и b
	static Monom LCM(const Monom& a, const Monom& b);
	// сравнения в разных упорядочениях
	static bool cmp_lex(const Monom& a, const Monom& b);
	static bool cmp_invlex(const Monom& a, const Monom& b);
	static bool cmp_grlex(const Monom& a, const Monom& b);
	static bool cmp_grevlex(const Monom& a, const Monom& b);
	static bool cmp_rinvlex(const Monom& a, const Monom& b);
};


std::string to_string(const Monom& m);

};

#endif	// #ifndef MONOM_H

