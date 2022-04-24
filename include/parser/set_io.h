#ifndef SET_IO_H
#define SET_IO_H

#include <sets/galois_2n.h>
#include <sets/residue.h>
#include <polynomial/polynom.h>

#include <istream>
#include <stdexcept>
#include <cctype>

#define PRIM_CHAR 'a'

template < typename CharT, typename Traits >
std::basic_istream< CharT, Traits >&
operator>> ( std::basic_istream< CharT, Traits >& is, Galois2N& val )
{
     char c;
     do
     {
          is.get( c );
     } while( isspace( c ) );

     if ( c == '0' )
     {
          val = Galois2N( val.irreducible_pol() );
     }
     else if ( c == '1' )
     {
          val = Galois2N( val.irreducible_pol(), 0 );
     }
     else if ( c == PRIM_CHAR )
     {
          val = Galois2N( val.irreducible_pol(), 1 );
     }
     else
     {
          throw std::runtime_error{ "wrong Galois field value" };
     }
     return is;
}


template < typename CharT, typename Traits >
std::basic_istream< CharT, Traits >&
operator>> ( std::basic_istream< CharT, Traits >& is, Residue& val )
{
     int64_t value;
     is >> value;
     val = Residue( val.get_modulo(), value );
     return is;
}

std::string to_string( const Galois2N& val )
{
     if ( !val )
     {
          return "0";
     }
     auto prim_pow = val.prim_power();
     if ( prim_pow == 0 )
     {
          return "1";
     }
     if ( prim_pow == 1 )
     {
          return std::string{ PRIM_CHAR };
     }
     return std::string{ PRIM_CHAR } + "^{" + std::to_string( prim_pow ) + "}";
}


std::string to_string( const Residue& val )
{
     return std::to_string( val.get_value() );
}


std::string to_string( const Monom& mon )
{
     const auto& vars = mon.get_vars();
     std::string result;
     for ( const auto& var : vars )
     {
          result += var.first;
          if ( var.second > 1 )
          {
               result += "^{" + std::to_string( var.second ) + "}";
          }
     }
     return result;
}


template < typename Polynom >
std::string to_string( const Polynom& pol )
{
     if ( !pol )
     {
          return "0";
     }
     using std::to_string;
     const auto& terms = pol.get_terms();
     std::string result;
     for ( const auto& term : terms )
     {
          if ( !result.empty() )
          {
               result += " + ";
          }
          std::string coeff = to_string( term.second );
          std::string mon = to_string( term.first );
          if ( coeff == "1" && !mon.empty() )
          {
               result += mon;
          }
          else
          {
               result += coeff;
               if ( !mon.empty() )
               {
                    result += "*" + mon;
               }
          }
     }
     return result;
}

#endif // #ifndef SET_IO_H
