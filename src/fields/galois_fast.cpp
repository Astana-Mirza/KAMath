#include "galois_fast.h"


namespace KAMath {

/**************************************************************************
** Количество незначащих нулей (нужно для нахождения степени полинома)
**************************************************************************/
size_t nlz(uint64_t x) {
	uint64_t y;
	size_t n, c;
	n = 64;
	c = 32;
	do {		// подсчёт двоичным поиском
		y = x >> c;
		if (y != 0) {
			n = n - c;
			x = y;
		}
		c = c >> 1;
	} while (c != 0);
	n -= x;
	return n;
}


/**************************************************************************
** Количество значащих подряд идущих нулей в конце числа (в двоичной системе)
**************************************************************************/
size_t ntz(uint64_t x) {
	return 64 - nlz(~x & (x-1));
}


/**************************************************************************
** Количество единиц в двоичной записи числа
**************************************************************************/
size_t pop(uint64_t x) {	// подсчёт двоичным поиском
	x = x - ((x >> 1) & 0x5555555555555555);
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F;
	x = x + (x >> 8);
	x = x + (x >> 16);
	x = x + (x >> 32);
	return x & 0x7F;
}


/**************************************************************************
** Деление по модулю полиномов, хранящих коэффициенты в виде битов в uint64_t
**************************************************************************/
uint64_t Field_FGF2::polynom_mod(uint64_t a, uint64_t b) {
	int64_t div_deg = polynom_deg(b);
	if (div_deg < 0)
		throw std::runtime_error{"Division by zero polynomial"};
	while (div_deg <= polynom_deg(a)) {
		a ^= b << (polynom_deg(a)-div_deg);
	}
	return a;
}


/**************************************************************************
** Получение координат в базисе, используя степень примитивного элемента поля
**************************************************************************/
uint64_t Field_FGF2::to_signature(int64_t a, uint64_t prim) {
	int64_t pr_deg = polynom_deg(prim);
	if (pr_deg < 0)
		throw std::runtime_error{"Division by zero polynomial"};
	uint64_t sig = 0;
	a %= (1 << pr_deg) - 1;
	// отрицательное a означает нулевой элемент, неотрицательное - степень
	if (a < 0)
		return 0;
	if (pr_deg > a)
		return 1 << a;
	if (a > 63) {	// если сигнатура для a не уместится в 64 бита
		int64_t b = a % 64, sig_deg;
		int64_t powers_left;
		// столько раз нужно умножить на себя результат для b,
		// взять остаток от деления и получить результат для а
		powers_left = a - b;
		sig = polynom_mod(1 << b, prim);
		while (powers_left) {
			sig_deg = polynom_deg(sig);
			if (powers_left >= 63-sig_deg) {
				sig <<= 63-sig_deg;
				powers_left -= 63-sig_deg;
			}
			else {
				sig <<= powers_left;
				powers_left = 0;
			}
			sig = polynom_mod(sig, prim);
		}
	}
	else {		// сигнатура а вмещается в 64 бита
		sig = polynom_mod(1 << a, prim);
	}
	return sig;
}


/**************************************************************************
** Получение степени примитивного элемента поля, используя координаты в базисе
**************************************************************************/
int64_t Field_FGF2::from_signature(uint64_t sig, uint64_t prim) {
	if (sig == 0) return -1;
	uint64_t add_power = 0;
	int64_t prim_lm = ntz(prim);
	while (pop(sig) > 1) {
		// понижать степень sig, чтобы младший член был как у prim
		// степень младшего члена примитивного полинома всегда
		// не меньше степени младшего члена сигнатуры, и pwr >= 0
		int64_t pwr = ntz(sig) - prim_lm;
		sig >>= pwr;
		add_power += pwr;
		sig ^= prim;
	}
	return ntz(sig) + add_power;
}


/**************************************************************************
** Создать матрицу перехода в нормальный базис
**************************************************************************/
std::vector<uint64_t> Field_FGF2::make_normal_matr(const int64_t beta_deg,
						   const uint64_t prim) {
	int64_t dim = polynom_deg(prim);
	if (dim < 2)
		throw std::runtime_error{"Wrong primitive polinomial"};
	std::vector<uint64_t> matr(dim);
	for (int64_t i = 0; i < dim; i++) {
		matr[i] = to_signature(beta_deg<<i, prim);
	}
	return matr;
}


/**************************************************************************
** Найти элемент степень элемента b нормального базиса
**************************************************************************/
int64_t Field_FGF2::find_beta_deg(const uint64_t prim) {
	size_t i;
	bool found = false;
	size_t lim = (1 << polynom_deg(prim))-1;
	for (i = 1; i < lim && !found; i++) {
		std::vector<uint64_t> matr = make_normal_matr(i, prim);
		if (is_lin_indep(matr)) found = true;
	}
	return i-1;
}


/**************************************************************************
** Проверка матрицы на обратимость
**************************************************************************/
bool Field_FGF2::is_lin_indep(const std::vector<uint64_t>& matr) {
	std::vector<uint64_t> m = matr;	// приведение к нижнетреугольной матрице
	size_t rk = m.size();
	for (size_t i = 0; i < rk; i++) {
		// если значение на диагонали не 1, то найти 1 в столбце
		if ((m[i] & (1 << (rk-i-1))) == 0) {
			bool found = false;
			size_t j = i+1;
			for (; j < rk && !found; j++) {
				if (m[j] & (1 << (rk-i-1))) found = true;
			}
			// если нашли, поменять строки, чтобы 1 была на диагонали
			if (found)
				std::swap(m[i], m[j-1]);
			// иначе имеем нулевой столбец, определитель равен нулю
			else
				return false;
		}
		// удаляем единицы столбца с помощью 1 на главной диагонали
		for (size_t j = i + 1; j < rk; j++) {
			if (m[j] & (1 << (rk-i-1))) {
				m[j] ^= m[i];
			}
		}
	}
	return true;
}


/**************************************************************************
** Обратная матрица, хранящаяся в виде битовых полей в uint64_t
**************************************************************************/
std::vector<uint64_t> Field_FGF2::inverse_matr(const std::vector<uint64_t>& matr) {
	size_t rk = matr.size();
	std::vector<uint64_t> inv(rk), m(matr);
	// создание единичной матрицы такого же порядка
	for (size_t i = 0; i < rk; i++) {
		inv[i] = 1 << (rk-i-1);
	}
	for (size_t i = 0; i < rk; i++) {
		// если значение на диагонали не 1, то найти 1 в столбце
		if ((m[i] & (1 << (rk-i-1))) == 0) {
			bool found = false;
			size_t j = i+1;
			for (; j < rk && !found; j++) {
				if (m[j] & (1 << (rk-i-1))) found = true;
			}
			// если нашли, поменять строки, чтобы 1 была на диагонали
			if (found) {
				std::swap(m[i], m[j-1]);
				std::swap(inv[i], inv[j-1]);
			}
			// иначе имеем нулевой столбец, определитель равен нулю
			else
				throw std::runtime_error{
					"Unable to find inverse matrix"};
		}
		// удаляем единицы столбца с помощью 1 на главной диагонали
		for (size_t j = 0; j < rk; j++) {
			if ((m[j] & (1 << (rk-i-1))) && i != j) {
				m[j] ^= m[i];
				inv[j] ^= inv[i];
			}
		}
	}
	return inv;
}


/**************************************************************************
** Оператор сложения
**************************************************************************/
Field_FGF2 Field_FGF2::operator + (const Field_FGF2& a) const {
	if (a.coord == 0) return *this;
	if (prim != a.prim || is_normal_basis != a.is_normal_basis)
		throw std::runtime_error{
		"Operands are in different fields or have different basis"};
	return Field_FGF2(prim, coord ^ a.coord,
			  GF2_UNUSED_DEGREE, is_normal_basis);
}


/**************************************************************************
** обратный элемент по умножению
**************************************************************************/
Field_FGF2 Field_FGF2::inv() const {
	if (!is_normal_basis && val_deg < 0)
		throw std::runtime_error{"Attempt to find inverse of 0"};
	if (is_normal_basis && val_deg == 0)
		return *this;
	uint64_t n = (1 << deg()) - 2;
	return pow(*this, n);
}


/**************************************************************************
** Оператор умножения
**************************************************************************/
Field_FGF2 Field_FGF2::operator * (const Field_FGF2& a) const {
	if (a.coord == 0)
		return Field_FGF2(prim, 0, GF2_UNUSED_DEGREE, is_normal_basis);
	if (prim != a.prim || is_normal_basis != a.is_normal_basis)
		throw std::runtime_error{
		"Operands are in different fields or have different basis"};
	dbit p1(128, coord), pr(128, prim), prod(128);
	uint64_t a_coord = a.coord;
	// умножение координат
	for (int64_t i = 0; a_coord != 0; i++) {
		prod = prod ^ (p1 << polynom_deg(a_coord));
		a_coord ^= (1<<polynom_deg(a_coord));
	}
	int64_t prod_deg = -1, prim_deg = deg();
	for (int64_t i = 127; i >= 0; i--) {
		if (prod[i]) {
			prod_deg = i;
			break;
		}
	}
	// деление по модулю на примитивный многочлен
	while (prod_deg >= prim_deg) {
		prod = prod ^ (pr << (prod_deg - prim_deg));
		for (int64_t i = prod_deg; i >= 0; i--) {
			if (prod.test(i)) {
				prod_deg = i;
				break;
			}
		}
	}
	return Field_FGF2(prim, prod.to_ulong(),
			  GF2_UNUSED_DEGREE, is_normal_basis);
}


/**************************************************************************
** Перевод из степенного базиса в нормальный
**************************************************************************/
void Field_FGF2::to_normal_basis(const std::vector<uint64_t>& matr) {
	if (matr.size() != static_cast<size_t>(deg()))
		throw std::runtime_error{"Wrong transformation matrix"};
	if (is_normal_basis)
		return;
	val_deg = GF2_UNUSED_DEGREE;
	is_normal_basis = true;
	int64_t rk = matr.size();
	uint64_t new_coord = 0;
	// умножение матрицы перехода на вектор-столбец
	for (int64_t i = 0; i < rk; i++) {
		// умножение столбца на столбец, т.к. матрица транспонирована
		for (int64_t j = 0; j < rk; j++) {
			uint64_t tmp = (matr[i] & (1 << j)) && (coord & (1 << (rk-i-1)));
			new_coord ^= tmp << (rk-j-1);
		}
	}
	coord = new_coord;
}


/**************************************************************************
** Перевод из нормального базиса в степенной
**************************************************************************/
void Field_FGF2::from_normal_basis(const std::vector<uint64_t>& matr) {
	if (matr.size() != static_cast<size_t>(deg()))
		throw std::runtime_error{"Wrong transformation matrix"};
	if (!is_normal_basis)
		return;
	is_normal_basis = false;
	int64_t rk = matr.size();
	uint64_t new_coord = 0;
	// умножение матрицы перехода на вектор-столбец
	for (int64_t i = 0; i < rk; i++) {
		// умножение столбца на столбец, т.к. матрица транспонирована
		// и координаты в обратном порядке
		for (int64_t j = 0; j < rk; j++) {
			uint64_t tmp = (matr[i] & (1 << j)) && (coord & (1 << i));
			new_coord ^= tmp << j;
		}
	}
	coord = new_coord;
	val_deg = from_signature(coord, prim);
}


/**************************************************************************
** Возведение в степень
**************************************************************************/
Field_FGF2 pow(const Field_FGF2& v, const int64_t p) {
	if (p < 0)
		return pow(v.inv(), -p);
	int q = p % ((1 << v.deg()) - 1);
	if (q == 0)
		return v / v;
	Field_FGF2 a = v;
	for (int i = 1; i < q; i++) a *= v;
	return a;
}


/**************************************************************************
** Проверка на ноль
**************************************************************************/
bool is_zero(const Field_FGF2& val) {
	return val.coords() == 0;
}


/**************************************************************************
** Строковое представление элемента поля Галуа
**************************************************************************/
std::string to_string(const Field_FGF2& m) {
	if (m.is_normal()) return std::string("B(") +
			std::to_string(m.coords()) + ")";
	if (m.value_power() == -1) return "0";
	if (m.value_power() == 0) return "1";
	if (m.value_power() == 1) return "a";
	return "a^{" + to_string(m.value_power()) + "}";
}

};

