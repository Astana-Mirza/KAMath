#ifndef FAUGERE_H
#define FAUGERE_H

#include "buchberger.h"
#include <list>

namespace KAMath {

// сигнатура - вектор полиномов, в котором единственное
// ненулевое значение - моном t под индексом index
struct Signature {
	size_t index;
	Monom t;

	Signature(size_t i=0, Monom mon={{}}) : index{i}, t{mon} {}
	Signature operator * (const Monom& m) const
		{ return {index, t*m}; }
	bool operator == (const Signature& s) const
		{ return index==s.index && t==s.t; }
	bool operator != (const Signature& s) const
		{ return !(*this==s); }
	// true, если a > b
	static bool cmp_sign(const Signature& a, const Signature& b,
							Ordering ord) {
		if (a.index == b.index) {
			switch(ord) {
			case Ordering::lex:
				return Monom::cmp_lex(a.t, b.t);
			case Ordering::invlex:
				return Monom::cmp_invlex(a.t, b.t);
			case Ordering::grlex:
				return Monom::cmp_grlex(a.t, b.t);
			case Ordering::grevlex:
				return Monom::cmp_grevlex(a.t, b.t);
			case Ordering::rinvlex:
				return Monom::cmp_rinvlex(a.t, b.t);
			}
		}
		return a.index < b.index;
	}
};


// помеченный полином - пара значений: полином и его сигнатура
template <typename F>
struct Label_pol {
	Signature sig;
	Polynom<F> pol;

	Label_pol(const Signature& s={}, const Polynom<F>& p={}):
		sig{s}, pol{p} {}
	size_t index() const { return sig.index; } // индекс сигнатуры
	Monom LM() const { return pol.LM(); }	// старший моном
	F LC() const { return pol.LC(); }	// старший коэффициент
	Label_pol mod(const std::vector<Polynom<F>>& G) const
		{ return {sig, pol.mod(G)}; }
	Label_pol operator * (const F& el) const
		{ return {sig, pol*el}; }
	Label_pol operator / (const F& el) const
		{ return {sig, pol/el}; }	// деление на элемент поля F
	Label_pol operator * (const Monom& m) const
		{ return {sig*m, pol*m}; }	// умножение на моном
	bool operator == (const Label_pol& lp) const
		{ return sig==lp.sig && pol==lp.pol; }
};


struct Rule {
	size_t k;	// правило упрощения полинома под номером k
	Monom t;	// вспомогательный моном
};


// критическая пара полиномов
template <typename F>
struct Crit_pair {
	bool is_valid;	// является ли пара непустой
	Monom lcm;	// наименьшее общее кратное
	Monom u1, u2;	// u1 * LT(r1) == u2 * LT(r2) == lcm
	size_t r1, r2;	// индексы полиномов в массиве

	cpp_int deg() const { return lcm.full_deg(); }
	// конструктор с проверкой критерия
	Crit_pair(size_t p1, size_t p2, size_t k,
			const std::vector<Polynom<F>>& G,
			const std::vector<Label_pol<F>>& pols) {
		Label_pol<F> pol1 = pols[p1], pol2 = pols[p2];
		Monom t = Monom::LCM(pol1.LM(), pol2.LM());
		Monom _u1 = t/pol1.LM(), _u2 = t/pol2.LM();
		if (Signature::cmp_sign(pol2.sig*_u2,
				pol1.sig*_u1, pol1.pol.get_ordering())) {
			// если первый больше, поменять местами
			std::swap(pol1, pol2);
			std::swap(_u1, _u2);
			std::swap(p1, p2);
		}
		if (pol1.sig.index > k) {
			is_valid = false;
			return;
		}
		auto mon = Polynom<F>({pol1.LC()/pol1.LC(), _u1*pol1.sig.t},
				      pol1.pol.get_ordering()).mod(G).LM();
		if (_u1*pol1.sig.t != mon) {
			is_valid = false;
			return;
		}
		mon = Polynom<F>({pol2.LC()/pol2.LC(), _u2*pol2.sig.t},
				 pol2.pol.get_ordering()).mod(G).LM();
		if (pol2.sig.index == k && _u2*pol2.sig.t != mon) {
			is_valid = false;
			return;
		}
		is_valid = true;	lcm = t;
		u1 = _u1;		u2 = _u2;
		r1 = p1;		r2 = p2;
	}
	bool operator >  (const Crit_pair& p) const { return deg() > p.deg(); }
	bool operator <  (const Crit_pair& p) const { return deg() < p.deg(); }
	bool operator == (const Crit_pair& p) const {
		return is_valid==p.is_valid && lcm==p.lcm &&
			u1==p.u1 && u2==p.u2 && r1==p.r1 && r2==p.r2; }
};


/**************************************************************************
** Удалить сигнатуры у вектора Label_pol и вернуть вектор полиномов
**************************************************************************/
template <typename F>
std::vector<Polynom<F>> poly(const std::vector<Label_pol<F>>& G) {
	std::vector<Polynom<F>> res;
	res.reserve(G.size());
	for (auto& g : G)
		res.push_back(g.pol);
	return res;
}


/**************************************************************************
** Нормальная форма полинома по отношению к набору G (остаток от деления,
** старший коэффициент равен 1)
**************************************************************************/
template <typename F>
Label_pol<F> normal_form(const Label_pol<F>& f,
			const std::vector<Label_pol<F>>& G) {
	auto p = f.mod(poly(G));
	if (p.pol.empty())
		return p;
	return p / p.LC();
}


/**************************************************************************
** Добавить правило для k-го полинома, т.е. в список вектора rules, который
** определяется индексом сигнатуры, будет добавлено правило {k, t} (t - моном
** в сигнатуре данного полинома)
**************************************************************************/
template <typename F>
void add_rule(std::vector<std::list<Rule>>& rules,
		const Label_pol<F>& lp, size_t k) {
	if (lp.sig.index == rules.size())
		rules.push_back({});
	// добавить новое правило для k-го полинома
	rules[lp.sig.index].push_back({k, lp.sig.t});
}


/**************************************************************************
** Проверка критерия перезаписывания для правил упрощения данного полинома
**************************************************************************/
template <typename F>
bool rewritten(const std::vector<std::list<Rule>>& rules, const Monom& u,
					const Label_pol<F>& lp, size_t k) {
	size_t saved_k = k;
	for (auto& rule : rules[lp.sig.index]) {
		// если уже есть правило упрощения
		// то оно перезаписано
		if ((u*lp.sig.t).is_divisible(rule.t)) {
			saved_k = rule.k;
			break;
		}
	}
	return saved_k != k;
}


/**************************************************************************
** Добавление S-полиномов и правил для них
**************************************************************************/
template <typename F>
std::vector<Label_pol<F>> s_pol(const std::vector<Crit_pair<F>>& P,
				std::vector<std::list<Rule>>& rules,
				std::vector<Label_pol<F>>& pols) {
	std::vector<Label_pol<F>> res;
	size_t n = pols.size()-1;
	for (size_t i = 0; i < P.size(); i++) {
		// проверка критерия перезаписывания
		if (P[i].is_valid &&
		!rewritten(rules, P[i].u1, pols[P[i].r1], P[i].r1) &&
		!rewritten(rules, P[i].u2, pols[P[i].r2], P[i].r2)) {
			n++;
			Label_pol<F> r{pols[P[i].r1].sig * P[i].u1,
					pols[P[i].r1].pol * P[i].u1 -
					pols[P[i].r2].pol * P[i].u2};
			add_rule(rules, r, n);	// добавить новые правила
			res.push_back(r);
			pols.push_back(r);
		}
	}
	std::sort(res.begin(), res.end(), [&](Label_pol<F> p1, Label_pol<F> p2) {
		return Signature::cmp_sign(p1.sig, p2.sig,
					p1.pol.get_ordering());
	});
	return res;
}


/**************************************************************************
** Проверка, редуцируется ли элемент (если да, вернёт редуктор,
** иначе пустой полином)
**************************************************************************/
template <typename F>
Label_pol<F> is_reducible(const std::vector<std::list<Rule>>& rules,
			const Label_pol<F>& r0,
			const std::vector<Label_pol<F>>& G,
			const std::vector<Label_pol<F>>& pols) {
	for (size_t i = 0; i < G.size(); i++) {
		if (!r0.LM().is_divisible(G[i].LM()))	// проверка делимости
			continue;
		Monom u = r0.LM() / G[i].LM();
		// проверка на нормализованность
		if (Polynom<F>({r0.LC()/r0.LC(), u*G[i].sig.t},
		r0.pol.get_ordering()).mod(poly(pols)).LM() != u*G[i].sig.t)
			continue;
		// можно ли использовать предыдущий результат для экономии
		if (rewritten(rules, u, G[i], i))
			continue;
		if (G[i].sig*u == r0.sig)	// устранить повторы
			continue;
		return G[i];
	}
	return Label_pol<F>{};
}


/**************************************************************************
** Редукция одного помеченного полинома
**************************************************************************/
template <typename F>
std::pair<Label_pol<F>, std::vector<Label_pol<F>>> top_reduction(
					std::vector<std::list<Rule>>& rules,
					const Label_pol<F>& r0,
					const std::vector<Label_pol<F>>& G,
					std::vector<Label_pol<F>>& pols) {
	// Многочлен редуцируется к нулю.
	// Не доказано, что в этом случае алгоритм завершится
	if (r0.pol.empty()) {
		throw std::runtime_error{
			"Warning: polynom is reduced to zero, may be infinite"};
	}
	size_t n = pols.size()-1;
	Label_pol<F> r = is_reducible(rules, r0, G, pols);
	if (r.pol.empty()) {
		// вернуть с пустым множеством (r0 не редуцируется к 0)
		return std::pair<Label_pol<F>, std::vector<Label_pol<F>>>(
							r0/r0.LC(), {});
	}
	Monom u = r0.LM()/r.LM();
	if (Signature::cmp_sign(r0.sig, r.sig*u, r0.pol.get_ordering())) {
		Label_pol<F> tmp(r0.sig, r0.pol-r.pol*u);
		 // вернуть непустое множество для его будущей редукции
		return std::pair<Label_pol<F>, std::vector<Label_pol<F>>>(
								{}, {tmp});
	}
	n++;
	Label_pol<F> tmp(r.sig*u, r.pol*u-r0.pol);
	pols.push_back(tmp);
	add_rule(rules, tmp, n);
	return std::pair<Label_pol<F>, std::vector<Label_pol<F>>>({}, {tmp, r0});
}


/**************************************************************************
** Редукция набора помеченных полиномов
**************************************************************************/
template <typename F>
std::vector<Label_pol<F>> reduction(std::vector<std::list<Rule>>& rules,
					const std::vector<Label_pol<F>>& G,
					std::vector<Label_pol<F>>& pols,
					std::vector<Label_pol<F>>& todo) {
	std::vector<Label_pol<F>> done;
	while (!todo.empty()) {
		Label_pol<F> h = normal_form(todo.back(), pols);
		todo.pop_back();
		std::vector<Label_pol<F>> tmp = G;
		tmp.insert(tmp.end(), done.begin(), done.end());
		auto pair = top_reduction(rules, h, tmp, pols);
		if (!pair.first.pol.empty())
			done.push_back(pair.first);
		else
			continue;
		todo.insert(todo.end(), pair.second.begin(), pair.second.end());
		std::sort(todo.begin(), todo.end(), [&](Label_pol<F> p1,
							Label_pol<F> p2) {
			return Signature::cmp_sign(p1.sig, p2.sig,
						  p1.pol.get_ordering());
		});
	}
	return done;
}


/**************************************************************************
** Функция добавляет в базис G новые элементы, используя массивы правил
** и критических пар
**************************************************************************/
template <typename F>
std::vector<Label_pol<F>> f5(std::vector<std::list<Rule>>& rules, size_t k,
			const Polynom<F>& f, std::vector<Label_pol<F>>& G) {
	Label_pol<F> r_i(Signature{k, {}}, f);
	std::vector<Label_pol<F>> Gi = G;
	Gi.push_back(r_i);
	std::vector<Crit_pair<F>> P;
	P.reserve(G.size());
	for (size_t i = 0; i < G.size(); i++) {
		Crit_pair<F> pair(G.size(), i, k, poly(G), Gi);
		if (pair.is_valid)
			P.push_back(pair);
	}
	std::sort(P.begin(), P.end());
	while (!P.empty()) {
		cpp_int d = P[0].deg();
		std::vector<Crit_pair<F>> Pd;
		auto it = P.begin();
		while(it != P.end()) {
			if (it->deg() == d) {
				Pd.push_back(*it);
				P.erase(it);
			} else ++it;
		}
		std::vector<Label_pol<F>> Fs = s_pol(Pd, rules, Gi);
		std::vector<Label_pol<F>> Rd = reduction(rules, Gi, G, Fs);
		for (auto& r : Rd) {
			Gi.push_back(r);
			for (size_t i = 0; i < Gi.size()-1; i++) {
				Crit_pair<F> pair(Gi.size()-1, i,
						  k, poly(Gi), Gi);
				if (pair.is_valid)
					P.push_back(pair);
			}
		}
		std::sort(P.begin(), P.end());
	}
	return Gi;
}


/**************************************************************************
** Функция для нахождения базиса Грёбнера алгоритмом F5 (F5 Incremental)
**************************************************************************/
template <typename F>
std::vector<Polynom<F>> basis_f5(const std::vector<Polynom<F>>& pols) {
	size_t m = pols.size();
	std::vector<std::list<Rule>> rules(m);
	Label_pol<F> r{{m-1, {}}, pols[m-1]};
	// базис для подмножеств из m ... 1 полиномов
	std::vector<Label_pol<F>> G{r};
	for (int64_t i = m-2; i >= 0; i--) {
		G = f5(rules, i, pols[i], G);
	}
	auto system = poly(G);
    reduce_basis(system);
	return system;
}


/**************************************************************************
** Функция для нахождения базиса Грёбнера алгоритмом F5C (F5C incremental)
**************************************************************************/
template <typename F>
std::vector<Polynom<F>> basis_f5c(const std::vector<Polynom<F>>& pols) {
	size_t m = pols.size();
	std::vector<std::list<Rule>> rules(m);
	Label_pol<F> r{{m-1, {}}, pols[m-1]};
	// базис для подмножеств из m ... 1 полиномов
	std::vector<Label_pol<F>> G{r};
	std::vector<Polynom<F>> B{pols[m-1]};
	for (int64_t i = m-2; i >= 0; i--) {
		std::vector<Label_pol<F>> G_new = f5(rules, i, pols[i], G);
		B = poly(G_new);	// редуцировать промежуточый базис
		reduce_basis(B);
		G = std::vector<Label_pol<F>>();
		G.reserve(B.size());
		rules = std::vector<std::list<Rule>>(B.size()+1);
		for (int64_t j = B.size()-1; j >= 0; j--) {
			Label_pol<F> pol{{B.size()-j, {}}, B[j]};
			G.push_back(pol);
		}
	}
	reduce_basis(B);
	return B;
}

};

#endif	// #ifndef FAUGERE_H

