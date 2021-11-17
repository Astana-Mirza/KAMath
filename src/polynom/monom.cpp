#include "monom.h"


namespace KAMath {

/**************************************************************************
** Конструктор, создающий моном из векторов переменных v и их степеней deg
**************************************************************************/
Monom::Monom(const std::vector<std::string>& v,
	     const std::vector<cpp_int>& deg) {
	if (deg.size() != v.size())
		throw std::runtime_error{"Wrong collection of variables"};
	for (size_t i = 0; i < deg.size(); i++) {
		vars[v[i]] = deg[i];
	}
}


/**************************************************************************
** Удалить переменную из монома. Если её нет, ничего не делать
**************************************************************************/
void Monom::remove_var(const std::string& var) {
	for (auto it = vars.begin(); it != vars.end(); ++it) {
		if (it->first == var) {
			vars.erase(it);
			return;
		}
	}
}


/**************************************************************************
** Полная степень (сумма степеней переменных)
**************************************************************************/
cpp_int Monom::full_deg() const {
	cpp_int deg = 0;
	for (auto& v : vars) {
		deg += v.second;
	}
	return deg;
}


/**************************************************************************
** Оператор умножения. Степени переменных складываются
**************************************************************************/
Monom Monom::operator * (const Monom& m) const {
	std::map<std::string, cpp_int> new_vars = vars;
	for (auto& v : m.vars) {
		new_vars[v.first] += v.second;
	}
	return Monom(new_vars);
}


/**************************************************************************
** Оператор деления. Степени переменных вычитаются
**************************************************************************/
Monom Monom::operator / (const Monom& m) const {
	std::map<std::string, cpp_int> new_vars = vars;
	for (auto& v : m.vars) {
		new_vars[v.first] -= v.second;
	}
	return Monom(new_vars);
}


/**************************************************************************
** Оператор сравнения полиномов
**************************************************************************/
bool  Monom::operator == (const Monom& m) const {
	for (auto& v : m.vars) {
		if (vars[v.first] != v.second)
			return false;
	}
	for (auto& v : vars) {
		if (m.vars[v.first] != v.second)
			return false;
	}
	return true;
}


/**************************************************************************
** Проверка на делимость
**************************************************************************/
bool Monom::is_divisible(const Monom& m) const {
	for (auto& v : m.vars) {
		if (vars[v.first] < v.second)
			return false;
	}
	return true;
}


/**************************************************************************
** Наименьшее общее кратное. Для каждой переменной берётся максимальная степень
**************************************************************************/
Monom Monom::LCM(const Monom& a, const Monom& b) {
	std::map<std::string, cpp_int> new_vars = a.vars;
	for (auto& v : b.vars) {
		new_vars[v.first]  = (new_vars[v.first] < v.second) ?
					v.second : new_vars[v.first];
	}
	return Monom(new_vars);
}


/**************************************************************************
** Сравнение мономов в различных упорядочениях. Возвращают значение a > b
**************************************************************************/
bool Monom::cmp_lex(const Monom& a, const Monom& b) {
	cpp_int cmp = 0;
	auto ita = a.vars.begin();
	auto itb = b.vars.begin();
	while (ita != a.vars.end() && itb != b.vars.end()) {
		// больше тот моном, у которого первая переменная меньше
		cmp = (ita->first).compare(itb->first);
		if (cmp != 0)
			return cmp < 0;
		cmp = ita->second - itb->second;
		if (cmp != 0)
			return cmp > 0;
		++ita; ++itb;
	}
	return false;
}


bool Monom::cmp_invlex(const Monom& a, const Monom& b) {
	cpp_int cmp = 0;
	auto ita = a.vars.rbegin();
	auto itb = b.vars.rbegin();
	while (ita != a.vars.rend() && itb != b.vars.rend()) {
		// больше тот моном, у которого последняя переменная больше
		cmp = (ita->first).compare(itb->first);
		if (cmp != 0)
			return cmp > 0;
		cmp = ita->second - itb->second;
		if (cmp != 0)
			return cmp > 0;
		++ita; ++itb;
	}
	return false;
}


bool Monom::cmp_grlex(const Monom& a, const Monom& b) {
	cpp_int ddeg = a.full_deg() - b.full_deg();
	if (ddeg > 0) return true;
	if (ddeg < 0) return false;
	return cmp_lex(a, b);
}


bool Monom::cmp_grevlex(const Monom& a, const Monom& b) {
	cpp_int ddeg = a.full_deg() - b.full_deg();
	if (ddeg > 0) return true;
	if (ddeg < 0) return false;
	return cmp_lex(a, b) && cmp_invlex(b, a);
}


bool Monom::cmp_rinvlex(const Monom& a, const Monom& b) {
	cpp_int cmp = 0;
	auto ita = a.vars.rbegin();
	auto itb = b.vars.rbegin();
	while (ita != a.vars.rend() && itb != b.vars.rend()) {
		// больше тот моном, у которого последняя переменная меньше
		cmp = (ita->first).compare(itb->first);
		if (cmp != 0)
			return cmp < 0;
		cmp = itb->second - ita->second;
		if (cmp != 0)
			return cmp > 0;
		++ita; ++itb;
	}
	return false;
}


/**************************************************************************
** Строковое представление монома
**************************************************************************/
std::string to_string(const Monom& m) {
	std::string repr("");
	for (auto& v : m.vars) {
		if (v.second > 0)  {
			repr += v.first;
			if (v.second > 1)
				repr += "^{" +
				boost::lexical_cast<std::string>(v.second)+"}";
		}
	}
	return repr;
}

};

