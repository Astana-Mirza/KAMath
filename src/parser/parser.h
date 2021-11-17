#ifndef PARSER_H
#define PARSER_H

#include "../polynom/polynom.h"
#include "field_reader.h"

namespace KAMParser {

enum class Kind : char {
	polynom, end,
 	plus='+', minus='-', mul='*', eq='=', lp='(', rp=')'
};


template <typename F>
struct Token {
	Kind kind;
	KAMath::Polynom<F> value;
};


template <typename F>
class Token_stream {
public:
	Token_stream(std::istream& s, F n, KAMath::Ordering order):
		ip{&s}, owns{false}, one{n}, ord{order} {}
	Token_stream(std::istream* p, F n, KAMath::Ordering order):
		ip{p}, owns{true}, one{n}, ord{order} {}

	~Token_stream() { close(); }

	Token<F> get();
	const Token<F>& current() { return ct; }
	cpp_int read_deg();
	std::string read_index();

	void set_input(std::istream& s) { close(); ip = &s; owns = false; }
	void set_input(std::istream* p) { close(); ip = p; owns = true; }
private:
	void close() { if (owns) delete ip; }

	std::istream* ip;
	bool owns;
	F one;
	KAMath::Ordering ord;
	Token<F> ct {Kind::end, {}};
};

template <typename F>
KAMath::Polynom<F> expr(const Token_stream<F>& ts, bool get);


/**************************************************************************
** Получение один входной элемент (скобку, операцию или полином)
**************************************************************************/
template <typename F>
Token<F> Token_stream<F>::get() {
	char ch;
	do { 						// пропуск пробелов
		if (!ip->get(ch))
			return ct={Kind::end, {}};
	} while (ch != '\n' && isspace(ch));

	switch (ch) {
		case '\n':
		case '\0': return ct = {Kind::end, {}};
		case '*':
		case '+':
		case '-':
		case '(':
		case ')':
		case '=':
			return ct={static_cast<Kind>(ch), {}};
	}
	F coeff = one;
	if (isdigit(ch) || isalpha(ch) || ch == '.') {
		// элемент поля F
		if (isdigit(ch) || ch == '.' || ch == IMAGINARY_UNIT ||
						ch == PRIM_CHAR) {
			ip->putback(ch);
			coeff = read_from_stream(ip, one);
			int64_t deg = read_deg();
			coeff = KAMath::pow(coeff, deg);
			ip->get(ch);
		}
		KAMath::Polynom<F> pol(ord);
		if (!KAMath::is_zero(coeff))
			pol.add_term({coeff, {}});
		// моном
		while (isalpha(ch) && ch != IMAGINARY_UNIT && ch != PRIM_CHAR) {
			std::string s = std::string()+ch, ind = read_index();
			if (ind != "")
				s += "_{" + ind + "}";
			cpp_int deg = read_deg();
			pol *= KAMath::Monom{{{s, deg}}};
			ip->get(ch);
		}
		ip->putback(ch);
		return ct={Kind::polynom, pol};
	}
	throw std::runtime_error{std::string("Invalid character: ") + ch};
}


/**************************************************************************
** Считывание степени
**************************************************************************/
template <typename F>
cpp_int Token_stream<F>::read_deg() {
	char ch;
	cpp_int deg = 1;
	char bracket = '\0';
	do { ip->get(ch); } while (isspace(ch) && ch != '\n');
	if (ch == '^') {
		do { ip->get(ch); } while (isspace(ch) && ch != '\n');
		if (ch == '{' || ch == '(')
			bracket = (ch == '{') ? '}' : ')';
		else if (!isdigit(ch))
			throw std::runtime_error{
				std::string("Unknown symbol: ") + ch};
		deg = isdigit(ch) ? ch-'0' : 0;
		do { ip->get(ch); } while (isspace(ch) && ch != '\n');
		while (isdigit(ch)) {
			deg = deg * 10 + (ch-'0');
			ip->get(ch);
		}
		ip->putback(ch);
		if (deg < 0)
			throw std::runtime_error{"Wrong degree"};
		if (bracket != '\0') {
			do { ip->get(ch); } while (isspace(ch) && ch != '\n');
			if (ch != bracket)
				throw std::runtime_error{
					std::string("\'") + bracket + "\' expected"};
		}
	} else ip->putback(ch);
	return deg;
}


/**************************************************************************
** Считывание индекса переменной
**************************************************************************/
template <typename F>
std::string Token_stream<F>::read_index() {
	char ch, bracket = '\0';
	std::string res = "";
	ip->get(ch);
	if (ch != '_') {
		ip->putback(ch);
		return res;
	}
	ip->get(ch);
	if (ch == '{' || ch == '(') {
		bracket = (ch == '{') ? '}' : ')';
		ip->get(ch);
		while (!ip->eof() && ch != bracket) {
			res += ch;
			ip->get(ch);
		}
	}
	else res += ch;
	return res;
}


/**************************************************************************
** Примитивный элемент считывания (не делится на сомножители или слагаемые)
**************************************************************************/
template <typename F>
KAMath::Polynom<F> prim(Token_stream<F>& ts, bool get) {
	if (get) ts.get();
	switch (ts.current().kind) {
	case Kind::polynom:
	{
		KAMath::Polynom<F> v = ts.current().value;
		ts.get();
		if (ts.current().kind == Kind::lp ||
		ts.current().kind == Kind::polynom)
			v *= prim(ts, false);
		return v;
	}
	case Kind::minus:
		return -prim(ts, true);
	case Kind::lp:
	{
		KAMath::Polynom<F> e = expr(ts, true);
		if (ts.current().kind != Kind::rp)
			throw std::runtime_error{"\')\' expected"};
		cpp_int d = ts.read_deg();
		e = e.pow(d);
		ts.get();
		if (ts.current().kind == Kind::lp ||
		ts.current().kind == Kind::polynom)
			e *= prim(ts, false);
		return e;
	}
	default:
		throw std::runtime_error{"Unexpected sequence of characters"};
	}
}


/**************************************************************************
** Элемент считывания (произведение примитивных элементов)
**************************************************************************/
template <typename F>
KAMath::Polynom<F> term(Token_stream<F>& ts, bool get) {
	KAMath::Polynom<F> left = prim(ts, get);
	for (;;) {
		switch (ts.current().kind) {
		case Kind::mul:
			left *= prim(ts, true);
			break;
		default:
			return left;
		}
	}
}


/**************************************************************************
** Выражение считывания (сумма или разность элементов считывания)
**************************************************************************/
template <typename F>
KAMath::Polynom<F> expr(Token_stream<F>& ts, bool get) {
	KAMath::Polynom<F> left = term(ts, get);
	for (;;) {
		switch (ts.current().kind) {
		case Kind::plus:
			left += term(ts, true);
			break;
		case Kind::minus:
		case Kind::eq:
			left -= term(ts, true);
			break;
		default:
			return left;
		}
	}
}

};

#endif	// #ifndef PARSER_H

