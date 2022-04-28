#ifndef INTERFACE_H
#define INTERFACE_H

#include <nonlinear_2n.h>
#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class Interface; }
QT_END_NAMESPACE

enum class GrobAlgo {
     Buchberger,
     BuchbergerImproved
};


enum class Mode {
     Gf2nEquation,
     Gf2nBasis,
     ResidueBasis
};


enum class MonomOrd {
     Lex,
     Invlex,
     Grlex,
     Grevlex
};


class Interface : public QMainWindow
{
     Q_OBJECT

public:
     Interface( QWidget *parent = nullptr );
     ~Interface();

private:
     Ui::Interface *ui_;
     Mode mode_ = Mode::Gf2nEquation;
     GrobAlgo algo_ = GrobAlgo::Buchberger;
     MonomOrd ordering_ = MonomOrd::Lex;

     void setActions();
     void readSystem();
     void showErrorMessage( QString errorMessage );
     void writeOutput( QStringList& output );

     template < typename Ordering >
     void readOrdered();

     template < typename Ordering >
     QStringList processGf2nEquation( const Galois2N::polynom_type& prim );

     template < typename Ordering >
     QStringList processGf2nBasis( const Galois2N::polynom_type& prim );

     template < typename Ordering >
     QStringList processResidueBasis( uint64_t mod );
};

#endif // INTERFACE_H
