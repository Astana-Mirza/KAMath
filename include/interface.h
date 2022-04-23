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


enum class SetType {
     Gf2n,
     Residue
};


enum class MonomOrd {
     Lex,
     Invlex,
     Grlex,
     Grevlex,
     Rinvlex
};


class Interface : public QMainWindow
{
     Q_OBJECT

public:
     Interface( QWidget *parent = nullptr );
     ~Interface();

private:
     Ui::Interface *ui_;
     SetType set_type_ = SetType::Gf2n;
     GrobAlgo algo_ = GrobAlgo::Buchberger;
     MonomOrd ordering_ = MonomOrd::Lex;

     void setActions();
     void readSystem();
     void showErrorMessage( QString errorMessage );
     void writeOutput( QStringList& output );

     template < typename Ordering >
     void readOrdered();

     template < typename Ordering >
     QStringList processGf2n( const Galois2N::polynom_type& prim );

     template < typename Ordering >
     QStringList processResidue( uint64_t mod );
};

#endif // INTERFACE_H
