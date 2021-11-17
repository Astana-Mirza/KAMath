#ifndef BUCHBERGER_H
#define BUCHBERGER_H

#include "../polynom/polynom.h"

namespace KAMath {

/**************************************************************************
** Создание S-полиномов
**************************************************************************/
template <typename F>
Polynom<F> s_pol(const Polynom<F>& f, const Polynom<F>& g) {
	if (f.empty() || g.empty())
		throw std::runtime_error{"Cannot find S-polynomial"};
	Polynom<F> flt(f.LT(), f.get_ordering());
	Polynom<F> glt(g.LT(), g.get_ordering());
	Polynom<F> x({f.LC() / f.LC(),
		Monom::LCM(f.LM(), g.LM())}, f.get_ordering());
	flt = Polynom<F>({x.LC()/flt.LC(), x.LM()/flt.LM()},
			 x.get_ordering());
	glt = Polynom<F>({x.LC()/glt.LC(), x.LM()/glt.LM()},
			 x.get_ordering());
	return (flt * f) - (glt * g);
}


/**************************************************************************
** Функция для нахождения базиса Грёбнера алгоритмом Бухбергера
**************************************************************************/
template <typename F>
std::vector<Polynom<F>> basis_buchberger(const std::vector<Polynom<F>>& pols) {
	std::vector<Polynom<F>> G = pols;
	std::vector<Polynom<F>> G_b;
	do {
		G_b = G;
		for (size_t i = 0; i < G_b.size()-1; i++) {
			for (size_t j = i+1; j < G_b.size(); j++) {
				Polynom<F> S = s_pol<F>(G_b[i], G_b[j]).mod(G_b);
				if (!S.empty() &&
				std::find(G.begin(), G.end(), S) == G.end()) {
					G.push_back(S);
				}
			}
		}
	} while (G != G_b);
	return G;
}


/**************************************************************************
** Критерий для улучшенного алгоритма Бухбергера
**************************************************************************/
template <typename F>
bool buchberger_criteria(const std::vector<std::pair<size_t, size_t>>& B,
			const size_t i, const size_t j, const Monom& lcm,
			const std::vector<Polynom<F>>& G) {
	for (size_t l = 0; l < G.size(); l++) {
		bool no_i = false, no_j = false;
		if (l > i)
			no_i = std::find(B.begin(), B.end(),
				std::pair<size_t, size_t>(i, l))==B.end();
		else if (l < i)
			no_i = std::find(B.begin(), B.end(),
				std::pair<size_t, size_t>(l, i))==B.end();
		if (!no_i)
			continue;

		if (l > j)
			no_j = std::find(B.begin(), B.end(),
				std::pair<size_t, size_t>(j, l))==B.end();
		else if (l < j)
			no_j = std::find(B.begin(), B.end(),
				std::pair<size_t, size_t>(l, j))==B.end();
		if (!no_j)
			continue;

		if (lcm.is_divisible(G[l].LM()))
			return true;
	}
	return false;
}


/**************************************************************************
** Функция для нахождения базиса Грёбнера улучшенным алгоритмом Бухбергера
**************************************************************************/
template <typename F>
std::vector<Polynom<F>> basis_buchberger_better(
			const std::vector<Polynom<F>>& pols) {
	size_t t = pols.size();
	std::vector<Polynom<F>> G = pols;
	std::vector<std::pair<size_t, size_t>> B;
	B.reserve((t*t-t)/2);
	for (size_t i = 0; i < t-1; i++) {
		for (size_t j = i+1; j < t; j++) {
			B.push_back(std::pair<size_t, size_t>(i, j));
		}
	}
	while (!B.empty()) {
		auto p = B.back();
		B.pop_back();
		Polynom<F> fi(G[p.first].LT(),  G[p.first].get_ordering());
		Polynom<F> fj(G[p.second].LT(), G[p.second].get_ordering());
		auto lcm =  Monom::LCM(fi.LM(), fj.LM());
		if (Polynom<F>({(fi*fj).LC(), lcm},
				fi.get_ordering()) != fi * fj &&
				!buchberger_criteria(B, p.first, p.second, lcm, G)) {
			Polynom<F> S = s_pol(G[p.first], G[p.second]).mod(G);
			if (!S.empty()) {
				t++;
				G.push_back(S);
				for (size_t i = 0; i < t-1; i++)
					B.push_back(std::pair<size_t, size_t>(i, t-1));
			}
		}
	}
	return G;
}


/**************************************************************************
** Редуцирование базиса
**************************************************************************/
template <typename F>
void reduce_basis(std::vector<Polynom<F>>& pols) {
	size_t p = 0;
	do {
		// делим каждый полином на все остальные
		std::vector<Polynom<F>> G;
		G.reserve(pols.size()-1);
		for (size_t i = 0; i < p; i++)
			G.push_back(pols[i]);
		for (size_t i = p + 1; i < pols.size(); i++)
			G.push_back(pols[i]);
		Polynom<F> pol = pols[p];
		Polynom<F> m = pol.mod(G);
		if (m.empty()) {
			// если делится без остатка, удаляем
			pols.erase(pols.begin()+p);
		} else {
			// иначе делим остаток на старший коэффициент
			pols[p] = m / m.LC();
			p++;
		}
	} while (p < pols.size());
}

};

#endif	// #ifndef BUCHBERGER_H

