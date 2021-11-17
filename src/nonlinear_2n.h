#ifndef NONLINEAR_2N
#define NONLINEAR_2N

#include "polynom/polynom.h"
#include "fields/galois_fast.h"
#include "fields/modulo_prime.h"

namespace KAMath {

// требуется матрица перехода от нормального базиса к степенному
Polynom<Field_FGF2> polynom_to_normal_basis(const std::vector<uint64_t>& matr,
					    const Polynom<Field_FGF2>& pol);
// требуется матрица перехода от степенного базиса к нормальному
Polynom<Field_FGF2> extend_vars(const Polynom<Field_FGF2>& pol,
				const std::vector<uint64_t>& matr);
Polynom<Field_ZP> simplify_degrees(const Polynom<Field_ZP>& pol);
// возвращает систему полиномов над полем {0; 1}
std::vector<Polynom<Field_ZP>> make_system(const Polynom<Field_FGF2>& pol);

};

#endif	// #ifndef NONLINEAR_2N
