#include "nonlinear_2n.h"

namespace KAMath {

/**************************************************************************
** Перевести коэффициенты полинома в нормальный базис
**************************************************************************/
Polynom<Field_FGF2> polynom_to_normal_basis(const std::vector<uint64_t>& matr,
					const Polynom<Field_FGF2>& pol) {
	Polynom<Field_FGF2> res(pol.get_ordering());
	for (auto tm : pol.terms) {
		tm.coeff.to_normal_basis(matr);
		res.add_term(Term<Field_FGF2>{tm.coeff, tm.mon});
	}
	res.sort_terms();
	return res;
}


/**************************************************************************
** Заменить переменные на их представление в нормальном базисе
**************************************************************************/
Polynom<Field_FGF2> extend_vars(const Polynom<Field_FGF2>& pol,
				const std::vector<uint64_t>& matr) {
	if (pol.empty()) return pol;
	int64_t prim_deg = pol.LC().deg();
	uint64_t prim = pol.LC().prim_pol();
	std::set<std::string> original_vars = pol.variables();
	Polynom<Field_FGF2> res = pol;
	for (auto& var : original_vars) {
		Field_FGF2 basis_el(prim, 1,
				GF2_UNUSED_DEGREE, true);
		basis_el.from_normal_basis(matr);
		Monom mon({{var+"_{0}", 1}});
		Polynom<Field_FGF2> var_basis({basis_el, mon},
					pol.get_ordering());
		for (int64_t i = 1; i < prim_deg; i++) {
			basis_el = Field_FGF2(prim, 1<<i,
					GF2_UNUSED_DEGREE, true);
			basis_el.from_normal_basis(matr);
			mon = Monom({{var+"_{"+to_string(i)+"}", 1}});
			var_basis += Polynom<Field_FGF2>({basis_el, mon},
						pol.get_ordering());
		}
		res = res.subst(var_basis, var);
	}
	return res;
}


/**************************************************************************
** Выставить степени переменных в 1
**************************************************************************/
Polynom<Field_ZP> simplify_degrees(const Polynom<Field_ZP>& pol) {
	Polynom<Field_ZP> res(pol.get_ordering());
	for (auto tm : pol.terms) {
		for (auto& v : tm.mon.vars) {	
			if (v.second > 1) v.second = 1;
		}
		res.add_term(tm);
	}
	res.sort_terms();
	return res;
}


/**************************************************************************
** Создать систему, приведя слагаемые при элементах базиса
**************************************************************************/
std::vector<Polynom<Field_ZP>> make_system(const Polynom<Field_FGF2>& pol) {
	if (pol.empty()) return std::vector<Polynom<Field_ZP>>();
	int64_t prim_deg = pol.LC().deg();
	std::vector<Polynom<Field_ZP>> system;
	std::set<std::string> vars = pol.variables();
	system.reserve(prim_deg + vars.size());
	Field_ZP one(2, 1, true);
	for (int64_t i = 0; i < prim_deg; i++) {
		// создание полинома системы
		system.push_back(Polynom<Field_ZP>(pol.get_ordering()));
		// заполнение его термами
		for (auto& t : pol.terms) {
			if (t.coeff.coords() & (1 << i)) {
				system[i].add_term(Term<Field_ZP>{one, t.mon});
			}
		}
		// упрощение степеней
		system[i] = simplify_degrees(system[i]);
	}
	// добавить уравнения вида x^2 + x = 0
	for (auto& v : vars) {
		Monom mon({{v, 2}});
		Polynom<Field_ZP> p({one, mon}, pol.get_ordering());
		mon.vars[v] = 1;
		p.add_term({one, mon});
		system.push_back(p);
	}
	// отсортировать мономы по упорядочению
	for (auto& pol : system) {
		pol.sort_terms();
	}
	return system;
}

};

