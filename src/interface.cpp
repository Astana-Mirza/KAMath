#include <interface.h>
#include <ui_interface.h>
#include <parser/parser.h>

#include <polynomial/grobner/buchberger.h>

#include <QString>
#include <QStringList>
#include <QMessageBox>
#include <QDebug>

#include <sstream>
#include <vector>
#include <ctime>
#include <stdexcept>

Interface::Interface( QWidget *parent ):
     QMainWindow( parent ), ui_( new Ui::Interface )
{
     ui_->setupUi( this );
     setActions();
     ui_->modInput->setPlaceholderText( tr( "Irreducible polynom" ) );
     ui_->modInput->setReadOnly( false );
}


Interface::~Interface()
{
     delete ui_;
}


void Interface::showErrorMessage( QString errorMessage )
{
     QMessageBox message;
     message.setWindowTitle( tr( "Error" ) );
     message.setIcon( QMessageBox::Warning );
     message.setText( tr( "Error occured when executng:\n" ) + errorMessage );
     message.exec();
}


void Interface::setActions()
{
     connect( ui_->calculateBtn, &QPushButton::clicked, this, &Interface::readSystem );
     connect( ui_->actionFieldGF, static_cast< void ( QAction::* )( bool ) >( &QAction::triggered ), this, [&]()
     {
          set_type_ = SetType::Gf2n;
          ui_->modInput->setText( tr( "" ) );
          ui_->modInput->setPlaceholderText( tr( "Irreducible polynom" ) );
          ui_->modInput->setReadOnly( false );
     } );
     connect( ui_->actionFieldZn, static_cast< void ( QAction::* )( bool ) >( &QAction::triggered ), this, [&]()
     {
          set_type_ = SetType::Residue;
          ui_->modInput->setText( tr("") );
          ui_->modInput->setPlaceholderText( tr( "Natural number" ) );
          ui_->modInput->setReadOnly( false );
     } );
     connect( ui_->actionBuch, static_cast< void ( QAction::* )( bool )>( &QAction::triggered ), this, [&]()
     {
          algo_ = GrobAlgo::Buchberger;
     } );
     connect( ui_->actionBuchI, static_cast< void ( QAction::* )( bool )>( &QAction::triggered ), this, [&]()
     {
          algo_ = GrobAlgo::BuchbergerImproved;
     } );
     connect( ui_->actionLEX, static_cast< void ( QAction::* )( bool ) >( &QAction::triggered ), this, [&]()
     {
          ordering_ = MonomOrd::Lex;
     } );
     connect( ui_->actionINVLEX, static_cast< void ( QAction::* )( bool ) >( &QAction::triggered ), this, [&]()
     {
          ordering_ = MonomOrd::Invlex;
     } );
     connect( ui_->actionRINVLEX, static_cast< void ( QAction::* )( bool ) >( &QAction::triggered ), this, [&]()
     {
          ordering_ = MonomOrd::Rinvlex;
     } );
     connect( ui_->actionGRLEX, static_cast< void ( QAction::* )( bool ) >( &QAction::triggered ), this, [&]()
     {
          ordering_ = MonomOrd::Grlex;
     } );
     connect( ui_->actionGREVLEX, static_cast< void ( QAction::* )( bool ) >( &QAction::triggered ), this, [&]()
     {
          ordering_ = MonomOrd::Grevlex;
     } );

     QActionGroup* group_field = new QActionGroup( this );
     ui_->actionFieldGF->setActionGroup( group_field );
     ui_->actionFieldZn->setActionGroup( group_field );
     group_field->setExclusive( true );

     QActionGroup* group_alg = new QActionGroup( this );
     ui_->actionBuch->setActionGroup( group_alg );
     ui_->actionBuchI->setActionGroup( group_alg );
     group_alg->setExclusive( true );
    
     QActionGroup* group_ord = new QActionGroup( this );
     ui_->actionLEX->setActionGroup( group_ord );
     ui_->actionINVLEX->setActionGroup( group_ord );
     ui_->actionRINVLEX->setActionGroup( group_ord );
     ui_->actionGRLEX->setActionGroup( group_ord );
     ui_->actionGREVLEX->setActionGroup( group_ord );
     group_ord->setExclusive( true );
    
     connect( ui_->quitAction, &QAction::triggered, this, &QApplication::quit );
}


void Interface::readSystem()
{
     switch ( ordering_ )
     {
          case MonomOrd::Lex:
          {
               readOrdered< LexGreater >();
               break;
          }
          case MonomOrd::Invlex:
          {
               readOrdered< InvlexGreater >();
               break;
          }
          case MonomOrd::Rinvlex:
          {
               readOrdered< RinvlexGreater >();
               break;
          }
          case MonomOrd::Grlex:
          {
               readOrdered< GrlexGreater >();
               break;
          }
          case MonomOrd::Grevlex:
          {
               readOrdered< GrevlexGreater >();
               break;
          }
     }
}


void Interface::writeOutput( QStringList& output )
{
     QString html = "<html><head></head><body>";
     QStringList left;
     QStringList rest;
     for ( auto& equation: output )
     {
          html += "<p>";
          equation.replace( "^{", "<span style=\"vertical-align: super; color: red;\">" );
          equation.replace( "_{", "<span style=\"vertical-align: sub; color: green;\">" );
          equation.replace(  "}", "</span>" );
          html += equation;
          html += "</p>";
     }
     html += "</body></html>";
     ui_->outputLabel->setTextFormat( Qt::RichText );
     ui_->outputLabel->setText( html );
}


template < typename Ordering >
void Interface::readOrdered()
{
     //try
     //{
          QString modstr = ui_->modInput->text();
          std::stringstream polst;
          QStringList result;

          switch ( set_type_ )
          {
               case SetType::Gf2n:
               {
                    modstr.replace( "a", "x" );
                    polst << modstr.toStdString() + "\n";
                    TokenStream< Residue, Ordering > ts( polst, Residue{ 2, 1 } );
                    auto pol = expr( ts, true );
                    if ( !pol || pol.get_terms().cbegin()->first.var_deg( "x" ) < 2 )
                    {
                         throw std::runtime_error{ "wrong primitive polynomial" };
                    }

                    qDebug() << QString::fromStdString( to_string( pol ) );

                    const auto& terms = pol.get_terms();
                    Galois2N::polynom_type prim( terms.cbegin()->first.full_deg() + 1 );
                    for ( const auto& term : terms )
                    {
                         prim.set( term.first.var_deg( "x" ) );
                    }
                    result = processGf2n< Ordering >( prim );
                    break;
               }
               case SetType::Residue:
               {
                    result = processResidue< Ordering >( modstr.toULongLong() );
                    break;
               }
          }
          writeOutput( result );
     //}
    //catch ( std::exception& ex )
     //{
      //    showErrorMessage( QString::fromStdString( ex.what() ) );
     //}
}


template < typename Ordering >
QStringList Interface::processGf2n( const Galois2N::polynom_type& prim )
{
     QString content = ui_->systemInput->toPlainText();
     std::stringstream ist;
     ist << content.toStdString() + "\n";
     Galois2N one{ prim, 0 };
     TokenStream< Galois2N, Ordering > ts( ist, one );
     std::vector< Polynom< Galois2N, Ordering > > pols;
     // std::vector< Polynom< Residue , Ordering > > system;
     clock_t begin = 0, end = 0;

     QStringList result;
     while ( ts.get().kind != Kind::end )
     {
          auto pol = expr( ts, false );
          if ( pol )
          {
               pols.push_back( pol );
          }
     }
     if ( pols.empty() )
     {
          return result;
     }

     auto system = make_system( pols.front() );
     switch ( algo_ ) {
          case GrobAlgo::Buchberger:
          {
              begin = clock();
              auto basis = Buchberger< Polynom< Residue, Ordering > >::find_basis_brute_force( system );
              reduce_basis( basis );
              end = clock();
              for ( auto& pol : basis )
              {
                  result.append( QString::fromStdString( to_string( pol ) ) + " = 0" );
              }
              break;
          }
          case GrobAlgo::BuchbergerImproved:
          {
               begin = clock();
               auto basis = Buchberger< Polynom< Residue, Ordering > >::find_basis( system );
               reduce_basis( basis );
               end = clock();
               for ( auto& pol : basis )
               {
                   result.append( QString::fromStdString( to_string( pol ) ) + " = 0" );
               }
               break;
          }
     }
     result.append( "Time elapsed:\n" + QString::number( ( end - begin ) / ( double )CLOCKS_PER_SEC ) );
     return result;
}


template < typename Ordering >
QStringList Interface::processResidue( uint64_t mod )
{
     QString content = ui_->systemInput->toPlainText();
     std::stringstream ist;
     ist << content.toStdString() + "\n";
     Residue one{ mod, 1 };

     TokenStream< Residue, Ordering > ts( ist, one );
     std::vector< Polynom< Residue, Ordering > > pols;
     clock_t begin = 0, end = 0;

     QStringList result;
     while ( ts.get().kind != Kind::end )
     {
          auto pol = expr(ts, false);
          if ( pol )
          {
               pols.push_back( pol );
          }
     }

     if ( pols.empty() )
     {
          return result;
     }

     switch ( algo_ ) {
          case GrobAlgo::Buchberger:
          {
              begin = clock();
              auto basis = Buchberger< Polynom< Residue, Ordering > >::find_basis_brute_force( pols );
              reduce_basis( basis );
              end = clock();
              for ( auto& pol : basis )
              {
                  result.append( QString::fromStdString( to_string( pol ) ) + " = 0" );
              }
              break;
          }
          case GrobAlgo::BuchbergerImproved:
          {
               begin = clock();
               auto basis = Buchberger< Polynom< Residue, Ordering > >::find_basis( pols );
               reduce_basis( basis );
               end = clock();
               for ( auto& pol : basis )
               {
                   result.append( QString::fromStdString( to_string( pol ) ) + " = 0" );
               }
               break;
          }
     }
     result.append( "Time elapsed:\n" + QString::number( ( end - begin ) / ( double )CLOCKS_PER_SEC ) );
     return result;
}
