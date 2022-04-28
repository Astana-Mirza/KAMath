#ifndef PARSER_H
#define PARSER_H

#include <polynomial/polynom.h>
#include <parser/set_io.h>

#include <string>
#include <istream>

enum class Kind : char
{
     polynom,
     end,
     plus = '+',
     minus = '-',
     mul = '*',
     eq = '=',
     lp = '(',
     rp = ')'
};


template < typename T >
struct Token
{
	Kind kind;
	T value;
};


template < typename T, typename Ordering >
class TokenStream
{
public:
     TokenStream( std::istream& stream, T one ) noexcept;
     TokenStream( std::istream* stream, T one ) noexcept;
     ~TokenStream();

     Token< Polynom< T, Ordering > > get();
     const Token< Polynom< T, Ordering > >& current();
     size_t read_deg();

     void set_input( std::istream& stream );
     void set_input( std::istream* stream );

private:
     void close();
     std::string read_index();

     Token< Polynom< T, Ordering > > ct_ { Kind::end, {} };
     T one_;
     std::istream* is_;
     bool owns_;
};


template < typename T, typename Ordering >
Polynom< T, Ordering > prim(const TokenStream< T, Ordering >& ts, bool get);

template < typename T, typename Ordering >
Polynom< T, Ordering > term(const TokenStream< T, Ordering >& ts, bool get);

template < typename T, typename Ordering >
Polynom< T, Ordering > expr(const TokenStream< T, Ordering >& ts, bool get);

//-----------------------------------------IMPLEMENTATION------------------------------------------

template < typename T, typename Ordering >
TokenStream< T, Ordering >::TokenStream( std::istream& stream, T one ) noexcept:
     one_{ one }, is_{ &stream }, owns_{ false } {}


template < typename T, typename Ordering >
TokenStream< T, Ordering >::TokenStream( std::istream* stream, T one ) noexcept:
     one_{ one }, is_{ stream }, owns_{ false } {}


template < typename T, typename Ordering >
TokenStream< T, Ordering >::~TokenStream()
{
     close();
}


template < typename T, typename Ordering >
Token< Polynom< T, Ordering > > TokenStream< T, Ordering >::get()
{
     char ch;
     do {
          if ( !is_->get( ch ) )
          {
               ct_ = { Kind::end, {} };
               return ct_;
          }
     } while ( ch != '\n' && isspace( ch ) );

     switch ( ch )
     {
          case '\n':
          case '\0':
          {
               ct_ = { Kind::end, {} };
               return ct_;
          }
          case '*':
          case '+':
          case '-':
          case '(':
          case ')':
          case '=':
          {
               ct_ = { static_cast< Kind >( ch ), {} };
               return ct_;
          }
          default: break;
     }
     T coeff = one_;
     if ( isdigit( ch ) || isalpha( ch ) )
     {
		if ( isdigit( ch ) || ch == PRIM_CHAR )
          {
			is_->putback( ch );
               *is_ >> coeff;
			uint64_t deg = read_deg();
			coeff = pow( coeff, deg );
               is_->get( ch );
          }
          Polynom< T, Ordering > pol( coeff );
          while ( isalpha( ch ) && ch != PRIM_CHAR ) {
               std::string var = std::string{ ch }, ind = read_index();
               if ( !ind.empty() )
               {
                    var += "_{" + ind + "}";
               }
               uint64_t deg = read_deg();
               pol *= Polynom< T, Ordering >{ { { Monom { { { var, deg } } }, one_ } } };
               is_->get( ch );
          }
          is_->putback( ch );
          ct_ = { Kind::polynom, pol };
          return ct_;
     }
     throw std::runtime_error{ std::string( "unexpected symbol: " ) + ch };
}


template < typename T, typename Ordering >
const Token< Polynom< T, Ordering > >& TokenStream< T, Ordering >::current()
{
     return ct_;
}


template < typename T, typename Ordering >
uint64_t TokenStream< T, Ordering >::read_deg()
{
     char ch;
     uint64_t deg = 1;
     char bracket = '\0';
     do
     {
          is_->get( ch );
     } while ( isspace( ch ) && ch != '\n' );

     if ( ch == '^' )
     {
          do {
               is_->get( ch );
          } while ( isspace( ch ) && ch != '\n' );
          if ( ch == '{' || ch == '(' )
          {
               bracket = ( ch == '{' ) ? '}' : ')';
          }
          else if ( !isdigit( ch ) )
          {
               throw std::runtime_error{ std::string( "unexpected symbol: " ) + ch };
          }
          else
          {
               is_->putback( ch );
          }

          *is_ >> deg;
          if ( bracket != '\0' )
          {
               do
               {
                    is_->get( ch );
               } while ( isspace( ch ) && ch != '\n' );
               if ( ch != bracket )
               {
                    throw std::runtime_error{ std::string( "\'" ) + bracket + "\' expected" };
               }
          }
     }
     else
     {
          is_->putback( ch );
     }
     return deg;
}


template < typename T, typename Ordering >
void TokenStream< T, Ordering >::set_input( std::istream& stream )
{
     close();
     is_ = &stream;
     owns_ = false;
}


template < typename T, typename Ordering >
void TokenStream< T, Ordering >::set_input( std::istream* stream )
{
     close();
     is_ = stream;
     owns_ = true;
}


template < typename T, typename Ordering >
void TokenStream< T, Ordering >::close()
{
     if ( owns_ )
     {
          delete is_;
     }
}


template < typename T, typename Ordering >
std::string TokenStream< T, Ordering >::read_index()
{
     char ch, bracket = '\0';
     std::string res;
     is_->get( ch );
     if (ch != '_')
     {
          is_->putback( ch );
          return res;
     }
     is_->get( ch );
     if ( ch == '{' || ch == '(' )
     {
          bracket = ( ch == '{' ) ? '}' : ')';
          is_->get( ch );
          while ( !is_->eof() && ch != bracket )
          {
               res += ch;
               is_->get( ch );
          }
     }
     else
     {
          res += ch;
     }
     return res;
}


template < typename T, typename Ordering >
Polynom< T, Ordering > prim( TokenStream< T, Ordering >& ts, bool get ) {
     if ( get )
     {
          ts.get();
     }
     switch ( ts.current().kind )
     {
          case Kind::polynom:
          {
               Polynom< T, Ordering > v = ts.current().value;
               ts.get();
               if ( ts.current().kind == Kind::lp || ts.current().kind == Kind::polynom )
               {
                    v *= prim( ts, false );
               }
               return v;
          }
          case Kind::minus:
          {
               return -prim( ts, true );
          }
          case Kind::lp:
          {
               Polynom< T, Ordering > e = expr( ts, true );
               if ( ts.current().kind != Kind::rp )
               {
                    throw std::runtime_error{ "\')\' expected" };
               }
               uint64_t d = ts.read_deg();
               e = pow( e, d );
               ts.get();
               if ( ts.current().kind == Kind::lp || ts.current().kind == Kind::polynom )
               {
                    e *= prim( ts, false );
               }
               return e;
          }
          case Kind::end:
          {
               throw std::runtime_error{ "Unexpected end of input" };
          }
          default:
          {
               throw std::runtime_error{ "Unexpected sequence of characters" };
          }
     }
     return ts.current().value;
}


template < typename T, typename Ordering >
Polynom< T, Ordering > term( TokenStream< T, Ordering >& ts, bool get )
{
     Polynom< T, Ordering > left = prim( ts, get );
     while ( true )
     {
          switch ( ts.current().kind )
          {
               case Kind::mul:
               {
                    left *= prim( ts, true );
                    break;
               }
               default:
               {
                    return left;
               }
          }
     }
     return left;
}


template < typename T, typename Ordering >
Polynom< T, Ordering > expr( TokenStream< T, Ordering >& ts, bool get )
{
     Polynom< T, Ordering > left = term( ts, get );
     while ( true )
     {
          switch ( ts.current().kind )
          {
               case Kind::plus:
               {
                    left += term( ts, true );
                    break;
               }
               case Kind::minus:
               case Kind::eq:
               {
                    left -= expr( ts, true );
                    break;
               }
               default:
               {
                    return left;
               }
          }
     }
     return left;
}

#endif // #ifndef PARSER_H
