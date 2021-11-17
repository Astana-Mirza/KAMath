#ifndef FIELD_READER_H
#define FIELD_READER_H

#include <cctype>
#include <exception>
#include <boost/multiprecision/cpp_int.hpp>
#include "../fields/galois_fast.h"
#include "../fields/modulo_prime.h"
#include "../fields/extra_fields.h"


// функции для считывания элементов разных полей
namespace KAMParser {

#define PRIM_CHAR 'a'
#define IMAGINARY_UNIT 'i'
using boost::multiprecision::cpp_int;


KAMath::Field_FGF2 read_from_stream(std::istream* ip, KAMath::Field_FGF2 one);
KAMath::Rational read_from_stream(std::istream* ip, KAMath::Rational one);
KAMath::Complex read_from_stream(std::istream* ip, KAMath::Complex one);
KAMath::Field_ZP read_from_stream(std::istream* ip, KAMath::Field_ZP one);

};

#endif	// #ifndef FIELD_READER_H
