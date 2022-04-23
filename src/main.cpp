#include <interface.h>

#include <QApplication>

int main( int argc, char *argv[] )
{
     QApplication a( argc, argv );
     Interface w;
     w.setWindowTitle( QObject::tr( "Groebner calc" ) );
     w.show();
     return a.exec();
}
