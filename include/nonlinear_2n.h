#ifndef NONLINEAR_2N
#define NONLINEAR_2N

#include <polynomial/polynom.h>
#include <sets/galois_2n.h>
#include <sets/residue.h>

#include <set>

template < typename Ordering >
std::set< std::string > get_all_vars( const Polynom< Galois2N, Ordering >& pol )
{
     std::set< std::string > vars;               // all variables in polynom

     const auto& terms = pol.get_terms();
     for ( const auto& term : terms )
     {
          const auto& monom_vars = term.first.get_vars();
          for ( const auto& var : monom_vars )
          {
               vars.insert( var.first );
          }
     }
     return vars;
}


// for all variables x in GF(2^n) sibstitute x = x_{0} + a*x_{1} + ... + a^{n-1}*x_{n-1}
template < typename Ordering >
Polynom< Galois2N, Ordering > extend_vars( const Polynom< Galois2N, Ordering >& pol )
{
     const auto irr_pol = pol.leading_coeff().irreducible_pol();
     size_t prim_pow = irr_pol.size() - 1;
     auto result = pol;
     std::set< std::string > vars = get_all_vars( pol );
     for ( const auto& var : vars )
     {
          Polynom< Galois2N, Ordering > extended;
          for ( size_t i = 0; i < prim_pow; i++ )
          {
               Monom mon{ { { var + "_{" + std::to_string( i ) + "}", 1 } } };
               extended += { { { mon, Galois2N( irr_pol, i ) } } };
          }
          result = result.subst( extended, var );
     }
     return result;
}


// for all variables x in GF(2) x == x^k, so set all degrees to 1
template < typename Ordering >
Polynom< Residue, Ordering > simplify_degrees( const Polynom< Residue, Ordering >& pol )
{
     Polynom< Residue, Ordering > result;
     const auto& terms = pol.get_terms();
     for ( const auto& term : terms )
     {
          Monom simplified = term.first;
          const auto& vars = simplified.get_vars();
          for ( auto& var : vars )
          {
               simplified.set_deg( var.first, 1 );
          }
          result += { { { simplified, term.second } } };
     }
     return result;
}

template < typename Ordering >
void add_extra_pols( const Polynom< Galois2N, Ordering >& pol, std::vector< Polynom< Residue, Ordering > >& system )
{
     std::set< std::string > vars = get_all_vars( pol );
     Residue one( 2, 1 );
     for ( const auto& var : vars )
     {
          system.push_back( Polynom< Residue, Ordering >{ {
               { Monom{ { { var, 2 } } } , one },
               { Monom{ { { var, 1 } } } , one }
          } } );
     }
}


// make polynomial from GF(2^n)[x] to system of polynomials from GF(2)[x_{0}, ..., x_{n-1}]
template < typename Ordering >
std::vector< Polynom< Residue, Ordering > > make_system( const Polynom< Galois2N, Ordering >& pol )
{
     auto extended = extend_vars( pol );
     std::vector< Polynom< Residue, Ordering > > result;
     size_t prim_pow = pol.leading_coeff().irreducible_pol().size() - 1;

     for ( size_t i = 0; i < prim_pow; i++ )
     {
          Polynom< Residue, Ordering > deg_composed;
          const auto& terms = extended.get_terms();
          for ( const auto& term : terms )
          {
               if ( term.second.coeffs().test( i ) )
               {
                    deg_composed += { { { term.first, Residue( 2, 1 ) } } };
               }
          }
          deg_composed = simplify_degrees( deg_composed );
          if ( deg_composed )
          {
               result.push_back( deg_composed );
          }
     }
     add_extra_pols( extended, result );
     return result;
}

#endif // #ifndef NONLINEAR_2N
