#ifndef GALOIS_FAST_H
#define GALOIS_FAST_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <boost/dynamic_bitset.hpp>


namespace KAMath {

using dbit = boost::dynamic_bitset<uint64_t>;
#define GF2_UNUSED_DEGREE -2

// найти количество ведущих незначащих нулей
size_t nlz(uint64_t x);
// найти количество последних подряд идущих значащих нулей
size_t ntz(uint64_t x);
// найти количество единичных битов в двоичной записи числа
size_t pop(uint64_t x);


// класс для быстрого вычисления в полях Галуа
// (поддерживает поля от GF(2^1) до GF(2^63))
class Field_FGF2 {
public:
	// преобразует степень примитивного члена в закодированные коэффициенты
	static uint64_t to_signature(int64_t a, uint64_t prim);
	// преобразует закодированные коэффициенты в степень примитивного члена
	static int64_t from_signature(uint64_t sig, uint64_t prim);
	// функция создаёт матрицу перехода от степенного базиса к нормальному
	// и возвращает её в транспонированном виде
	// beta_deg - степень примитивного элемента, такая что beta = a^beta_deg.
	static std::vector<uint64_t> make_normal_matr(const int64_t beta_deg,
							const uint64_t prim);
	// проверка на линейную зависимость (не равен ли определитель нулю)
	static bool is_lin_indep(const std::vector<uint64_t>& matr);
	// нахождение обратной матрицы. Для транспонированной матрицы
	// обртная матрица тоже будет транспонированной
	static std::vector<uint64_t> inverse_matr(
				const std::vector<uint64_t>& matr);
	// поиск элемента beta нормального базиса
	// (возвращает степень примитивного элемента)
	static int64_t find_beta_deg(const uint64_t prim);
	// конструкторы
	Field_FGF2(const uint64_t p=7, const uint64_t val=0,
		   const int64_t v_deg=GF2_UNUSED_DEGREE, bool is_nb=false):
		prim{p}, coord{val},
		val_deg{(v_deg == GF2_UNUSED_DEGREE && !is_nb) ?
		from_signature(val, p) : v_deg % ((1 << polynom_deg(p))-1)},
		is_normal_basis{is_nb} {}
	Field_FGF2(const Field_FGF2& n) : prim{n.prim}, coord{n.coord},
		val_deg{n.val_deg}, is_normal_basis{n.is_normal_basis} {}
	void operator = (const Field_FGF2& n) {
		prim = n.prim;
		coord = n.coord;
		val_deg = n.val_deg;
		is_normal_basis = n.is_normal_basis;
	}

	// возвращает степень примитивного полинома
	int64_t deg() const { return 63-nlz(prim); }
	// возвращает значение элемента поля как степень примитивного элемента
	int64_t value_power() const { return val_deg; }
	uint64_t coords() const { return coord; }
	uint64_t prim_pol() const { return prim; }
	bool is_normal() const { return is_normal_basis; }
	// обратный элемент по умножению
	Field_FGF2 inv() const;
	Field_FGF2 operator + (const Field_FGF2& a) const;
	Field_FGF2 operator * (const Field_FGF2& a) const;
	Field_FGF2 operator / (const Field_FGF2& a) const { return *this*a.inv(); }
	// унарный - (из-за характеристики 2 то же, что и +)
	Field_FGF2 operator - () const { return *this; }
	Field_FGF2 operator - (const Field_FGF2& a) const { return *this+a; }
	bool operator == (const Field_FGF2& a) const {
		return (prim == a.prim) && (coord == a.coord) &&
			(val_deg == a.val_deg) &&
			(is_normal_basis == a.is_normal_basis);
	}

	Field_FGF2 operator += (const Field_FGF2& a) {*this = *this+a; return *this; }
	Field_FGF2 operator *= (const Field_FGF2& a) {*this = *this*a; return *this; }
	Field_FGF2 operator /= (const Field_FGF2& a) {*this = *this/a; return *this; }
	Field_FGF2 operator -= (const Field_FGF2& a) {*this = *this-a; return *this; }

// перевести в нормальный базис с помощью транспонированной матрицы перехода
void to_normal_basis(const std::vector<uint64_t>& matr);
// перевести в степенной базис с помощью транспонированной матрицы перехода
void from_normal_basis(const std::vector<uint64_t>& matr);

private:
	/* Коэффициенты многочленов в поле GF(2^n) характеристики 2
	 * принимают значения 0 или 1,
	 * поэтому коэффициенты хранятся в виде бит числа.*/
	uint64_t prim;	// коэффициенты примитивного многочлена данного поля
	uint64_t coord;	// координаты значения элемента поля в
	// степенном базисе 1, x, x^2, ..., или в нормальном
	int64_t val_deg;
	bool is_normal_basis;

	static int64_t polynom_deg(const uint64_t pol) { return 63-nlz(pol); }
	// остаток от деления полинома a на b
	static uint64_t polynom_mod(uint64_t a, uint64_t b);
};

// возведение в степень p
Field_FGF2 pow(const Field_FGF2& v, const int64_t p);
bool is_zero(const Field_FGF2& val);

using std::to_string;
std::string to_string(const Field_FGF2& m);

};

#endif		// #ifndef GALOIS_FAST_H

