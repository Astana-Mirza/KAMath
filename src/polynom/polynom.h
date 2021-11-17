#ifndef POLYNOM_H
#define POLYNOM_H

#include "monom.h"

namespace KAMath {

enum class Ordering { lex, invlex, grlex, grevlex, rinvlex};


// одночлен с коэффициентом, полином хранит вектор с ними
template <typename F>
struct Term {
	F coeff;
	Monom mon;

	Term operator -  () const { return Term{-coeff, mon}; }
	bool operator == (const Term& t) const {
		return coeff==t.coeff && mon==t.mon;
	}
};


/**************************************************************************
** Строковое представление одночлна с коэффициентом
**************************************************************************/
template <typename F>
std::string to_string(const Term<F>& t) {
	std::string cf = to_string(t.coeff), mn = to_string(t.mon);
	if (cf == "1" && mn != "")
		return mn;
	if (cf == "-1" && mn != "")
		return "-"+mn;
	return cf + mn;
}


template <typename F>		// многочлен над полем F
class Polynom {
public:
	std::vector<Term<F>> terms;

	Polynom(Ordering order=Ordering::lex) : terms{}, ord{order} {}
	Polynom(const Term<F>& t, Ordering order=Ordering::lex):
					terms{t}, ord{order} {}

	// тождественен ли многочлен нулю
	bool empty() const { return terms.empty(); }
	Ordering get_ordering() const { return ord; }
	void add_term(const Term<F>& t);
	// отсортировать члены по убыванию в упорядочении ord
	void sort_terms();
	void change_ordering(Ordering n) { ord = n; sort_terms(); }
	// вернуть полином с другим упорядочением
	Polynom changed_ordering(Ordering n) const
		{ Polynom pol{*this}; pol.change_ordering(n); return pol; }
	std::set<std::string> variables() const;

	bool operator == (const Polynom& p) const
		{ return terms==p.terms && ord==p.ord; }
	bool operator != (const Polynom& p) { return !(*this == p); }
	Polynom operator + (const Polynom& p) const;
	Polynom operator - (const Polynom& p) const;
	Polynom operator - () const;
	Polynom operator * (const Polynom& p) const;
	Polynom operator * (const F& c) const;
	Polynom operator * (const Monom& m) const;
	Polynom operator / (const F& c) const;

	Polynom operator += (const Polynom& p) { *this = *this + p; return *this; }
	Polynom operator -= (const Polynom& p) { *this = *this - p; return *this; }
	Polynom operator *= (const Polynom& p) { *this = *this * p; return *this; }
	Polynom operator *= (const F& c) { *this = *this * c; return *this; }
	Polynom operator *= (const Monom& m) { *this = *this * m; return *this; }
	Polynom operator /= (const F& c) { *this = *this / c; return *this; }
	Polynom pow(cpp_int p) const;	// возвести в степень p
	// старший член
	Term<F> LT() const{ return empty() ? Term<F>{{}, {}} : terms[0]; }
	Monom LM() const { return LT().mon; }	// старший моном
	F LC() const { return LT().coeff; }	// старший коэффициент
	Polynom mod(const std::vector<Polynom>& divs) const;
	// подстановка полинома pol в переменную var
	Polynom subst(const Polynom<F>& pol, const std::string& var) const;
private:
	Ordering ord;	// используемое упорядочение

	static bool cmp_lex(const Term<F>& a, const Term<F>& b)
		{ return Monom::cmp_lex(a.mon, b.mon); }
	static bool cmp_invlex(const Term<F>& a, const Term<F>& b)
		{ return Monom::cmp_invlex(a.mon, b.mon); }
	static bool cmp_grlex(const Term<F>& a, const Term<F>& b)
		{ return Monom::cmp_grlex(a.mon, b.mon); }
	static bool cmp_grevlex(const Term<F>& a, const Term<F>& b)
		{ return Monom::cmp_grevlex(a.mon, b.mon); }
	static bool cmp_rinvlex(const Term<F>& a, const Term<F>& b)
		{ return Monom::cmp_rinvlex(a.mon, b.mon); }
};


/**************************************************************************
** Добавление одночлена к многочлену
**************************************************************************/
template <typename F>
void Polynom<F>::add_term(const Term<F>& t) {
	bool found = false;	// ищем такой моном
	for (auto tm = terms.begin(); tm != terms.end(); ++tm) {
		if (tm->mon == t.mon) {
			// если нашли, сложить коэффициенты
			tm->coeff += t.coeff;
			if (is_zero(tm->coeff)) {
				terms.erase(tm);
				return;
			}
			found = true;
			break;
		}
	}
	if (!found) {	// если не нашли, добавим его в полином
		terms.push_back(t);
	}
}


/**************************************************************************
** Сортировка членов полинома в соответствии с выбранным упорядочением
**************************************************************************/
template <typename F>
void Polynom<F>::sort_terms() {
	switch (ord) {
		case Ordering::lex:
			std::sort(terms.begin(), terms.end(), cmp_lex);
			break;
		case Ordering::invlex:
			std::sort(terms.begin(), terms.end(), cmp_invlex);
			break;
		case Ordering::grlex:
			std::sort(terms.begin(), terms.end(), cmp_grlex);
			break;
		case Ordering::grevlex:
			std::sort(terms.begin(), terms.end(), cmp_grevlex);
			break;
		case Ordering::rinvlex:
			std::sort(terms.begin(), terms.end(), cmp_rinvlex);
			break;
	}
}


/**************************************************************************
** Получить множество всех переменных полинома
**************************************************************************/
template <typename F>
std::set<std::string> Polynom<F>::variables() const {
	std::set<std::string> res;
	for (auto& tm : terms) {
		for (auto& v : tm.mon.vars) {
			res.insert(v.first);
		}
	}
	return res;
}


/**************************************************************************
** Сложение полиномов
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::operator + (const Polynom<F>& p) const {
	Polynom<F> res = *this;
	for (auto& tm : p.terms) {
		res.add_term(tm);
	}
	res.sort_terms();
	return res;
}


/**************************************************************************
** Вычитание полинома
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::operator - (const Polynom<F>& p) const {
	Polynom<F> res = *this;
	for (auto& tm : p.terms) {
		res.add_term(-tm);
	}
	res.sort_terms();
	return res;
}


/**************************************************************************
** Унарный минус
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::operator - () const {
	Polynom<F> res(ord);
	for (auto& tm : terms) res.add_term(-tm);
	res.sort_terms();
	return res;
}


/**************************************************************************
** Произведение полиномов
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::operator * (const Polynom<F>& p) const {
	Polynom<F> res(ord);
	for (auto& t1 : terms) {
		for (auto& t2 : p.terms) {
			res.add_term({t1.coeff*t2.coeff, t1.mon*t2.mon});
		}
	}
	res.sort_terms();
	return res;
}


/**************************************************************************
** Произведение полинома на константу
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::operator * (const F& c) const {
	Polynom<F> res(ord);
	for (auto& t : terms) {
		res.add_term({t.coeff*c, t.mon});
	}
	return res;
}


/**************************************************************************
** Произведение полинома на моном
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::operator * (const Monom& m) const {
	Polynom<F> res(ord);
	for (auto& t : terms) {
		res.add_term({t.coeff, t.mon*m});
	}
	return res;
}


/**************************************************************************
** Деление полинома на константу
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::operator / (const F& c) const {
	if (is_zero(c)) throw std::runtime_error{"Division by zero"};
	Polynom<F> res(ord);
	for (auto& t : terms) {
		res.add_term({t.coeff/c, t.mon});
	}
	return res;
}


/**************************************************************************
** Возведение полинома в степень
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::pow(cpp_int p) const {
	if (empty())
		return *this;
	if (p == 0)	// представление единицы
		return Polynom<F>{{terms[0].coeff/terms[0].coeff, {}}, ord};
	Polynom<F> res = *this;
	for (cpp_int i = 1; i < p; i++) {
		res *= *this;
	}
	return res;
}


/**************************************************************************
** Деление на набор полиномов с остатком, возвращает остаток
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::mod(const std::vector<Polynom<F>>& divs) const {
	Polynom<F> p = *this, res(ord);
	while (!p.empty()) {
		size_t i = 0;
		bool have_division = false;
		while (i < divs.size() && !have_division) {
			if (p.LM().is_divisible(divs[i].LM())) {
				// единица в поле F находится
				// делением коэффициента на себя
				Polynom<F> tmp({p.LC() / divs[i].LC(),
					p.LM() / divs[i].LM()}, ord);
				p -= tmp*divs[i];
				have_division = true;
			}
			else i++;
		}
		if (!have_division) {
			Polynom<F> tmp(p.LT(), ord);
			res += tmp;
			p -= tmp;
		}
	}
	return res;
}


/**************************************************************************
** Замена переменной var на полином pol (pol может быть числом - частный случай)
**************************************************************************/
template <typename F>
Polynom<F> Polynom<F>::subst(const Polynom<F>& pol,
				const std::string& var) const {
	if (empty()) return *this;
	Polynom<F> res(ord);
	for (auto tm : terms) {
		// сохранить степень, в которой была переменная
		cpp_int pwr = tm.mon.vars[var];
		// удалить заменяемую переменную
		tm.mon.remove_var(var);
		res += pol.pow(pwr) * Polynom<F>(tm, ord);
	}
	return res;
}


/**************************************************************************
** Строковое представление полинома
**************************************************************************/
template <typename F>
std::string to_string(const Polynom<F>& pol) {
	if (pol.empty()) return "0";
	std::string str("");
	for (const Term<F>& t : pol.terms) {
		std::string add = to_string(t);
		if (str != "") {
			if (add[0] == '-') {
				add[0] = ' ';
				str += " -";
			}
			else
				str += " + ";
		}
		str += add;
	}
	return str;
}

};

#endif	// #ifndef POLYNOM_H
