#ifndef INTERFACE_H
#define INTERFACE_H

#include <QMainWindow>
#include <QString>
#include <QStringList>
#include <sstream>
#include <vector>
#include "fields/extra_fields.h"
#include "grobner_alg/faugere.h"
#include "nonlinear_2n.h"
#include "parser/parser.h"
#include <ctime>
#include <stdexcept>
#include <QMessageBox>

using namespace KAMath;
using namespace KAMParser;

QT_BEGIN_NAMESPACE
namespace Ui { class Interface; }
QT_END_NAMESPACE

enum class GrobAlgo {Buchberger, BuchbergerBetter, F5, F5C};
enum class FieldTypes {RationalT, ComplexT, GF_2T, ZpT};

class Interface : public QMainWindow
{
    Q_OBJECT

public:
    Interface(QWidget *parent = nullptr);
    ~Interface();

    FieldTypes field_type = FieldTypes::RationalT;
    GrobAlgo algo = GrobAlgo::Buchberger;
    Ordering ord_type = Ordering::lex;
private:
    Ui::Interface *ui;

    void setActions();
    void readSystem();
    void showErrorMessage(QString errorMessage);
    void writeOutput(QStringList& output);
};

#endif // INTERFACE_H
