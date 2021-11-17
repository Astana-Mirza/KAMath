#include "interface.h"
#include "ui_interface.h"

Interface::Interface(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::Interface) {
    ui->setupUi(this);
    setActions();
}

Interface::~Interface() {
    delete ui;
}

void Interface::showErrorMessage(QString errorMessage)
{
    QMessageBox message;
    message.setWindowTitle("Ошибка!");
    message.setIcon(QMessageBox::Warning);
    message.setText("В процессе исполнения возникла ошибка:\n" + errorMessage);

    message.exec();
}

void Interface::setActions() {
    connect(ui->calculateBtn, &QPushButton::clicked, this, &Interface::readSystem);
    connect(ui->actionFieldR, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){
        field_type = FieldTypes::RationalT;
        ui->modInput->setText("");
        ui->modInput->setPlaceholderText("");
        ui->modInput->setReadOnly(true);
    });
    connect(ui->actionFieldC, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){
        field_type = FieldTypes::ComplexT;
        ui->modInput->setText("");
        ui->modInput->setPlaceholderText("");
        ui->modInput->setReadOnly(true);
    });
    connect(ui->actionFieldGF, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){
        field_type = FieldTypes::GF_2T;
        ui->modInput->setText("");
        ui->modInput->setPlaceholderText("Неприводимый полином");
        ui->modInput->setReadOnly(false);
    });
    connect(ui->actionFieldZn, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){
        field_type = FieldTypes::ZpT;
        ui->modInput->setText("");
        ui->modInput->setPlaceholderText("Простое число");
        ui->modInput->setReadOnly(false);
    });

    connect(ui->actionBuch, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ algo = GrobAlgo::Buchberger; });
    connect(ui->actionBuchB, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ algo = GrobAlgo::BuchbergerBetter; });
    connect(ui->actionF5, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ algo = GrobAlgo::F5; });
    connect(ui->actionF5C, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ algo = GrobAlgo::F5C; });

    connect(ui->actionLEX, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ ord_type = Ordering::lex; });
    connect(ui->actionINVLEX, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ ord_type = Ordering::invlex; });
    //connect(ui->actionRINVLEX, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ ord_type = Ordering::rinvlex; });
    connect(ui->actionGRLEX, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ ord_type = Ordering::grlex; });
    connect(ui->actionGREVLEX, static_cast<void (QAction::*)(bool)>(&QAction::triggered), this, [&](){ ord_type = Ordering::grevlex; });


    QActionGroup* group_field = new QActionGroup(this);
    ui->actionFieldR->setActionGroup(group_field);
    ui->actionFieldC->setActionGroup(group_field);
    ui->actionFieldGF->setActionGroup(group_field);
    ui->actionFieldZn->setActionGroup(group_field);
    group_field->setExclusive(true);


    QActionGroup* group_alg = new QActionGroup(this);
    ui->actionBuch->setActionGroup(group_alg);
    ui->actionBuchB->setActionGroup(group_alg);
    ui->actionF5->setActionGroup(group_alg);
    ui->actionF5C->setActionGroup(group_alg);
    group_alg->setExclusive(true);
    
    QActionGroup* group_ord = new QActionGroup(this);
    ui->actionLEX->setActionGroup(group_ord);
    ui->actionINVLEX->setActionGroup(group_ord);
    //ui->actionRINVLEX->setActionGroup(group_ord);
    ui->actionGRLEX->setActionGroup(group_ord);
    ui->actionGREVLEX->setActionGroup(group_ord);
    group_ord->setExclusive(true);
    
    connect(ui->quitAction, &QAction::triggered, this, &QApplication::quit);
}


template <typename F>
QStringList rInput(std::istream& ist, F one, Ordering order, GrobAlgo algo = GrobAlgo::F5) {
    Token_stream<F> ts(ist, one, order);
    std::vector<Polynom<F>> vec;
    clock_t begin = 0, end = 0;

    QStringList result;
    while (ts.get().kind != Kind::end) {
      //!ts.get().value.empty()
      auto pol = expr(ts, false);
      pol.sort_terms();
      //std::cout << to_string(pol) << std::endl;
      bool unique = true;
      for (auto& p : vec) {
          if (p == pol) { unique = false; break; }
      }
      if (unique) {
          vec.push_back(pol/pol.LC());
      }
   }

    if (vec.empty()) return result;
    switch (algo) {
         case GrobAlgo::Buchberger:
         {
             begin = clock();
             auto s2 = basis_buchberger<F>(vec);
             reduce_basis(s2);
             end = clock();
             for (auto& pol : s2) {
                 result.append(QString::fromStdString(to_string(pol)) + " = 0");
                 //std::cout << to_string(pol) << std::endl;
             }
             break;
         }
         case GrobAlgo::BuchbergerBetter:
         {
             begin = clock();
             auto s2 = basis_buchberger_better<F>(vec);
             reduce_basis(s2);
             end = clock();
             for (auto& pol : s2) {
                 result.append(QString::fromStdString(to_string(pol)) + " = 0");
                 //std::cout << to_string(pol) << std::endl;
             }
             break;
         }
         case GrobAlgo::F5:
         {
             begin = clock();
             auto s2 = basis_f5<F>(vec);
             reduce_basis(s2);
             end = clock();
             for (auto& pol : s2) {
                 result.append(QString::fromStdString(to_string(pol)) + " = 0");
                 //std::cout << to_string(pol) << std::endl;
             }
             break;
         }
         case GrobAlgo::F5C:
         {
             begin = clock();
             auto s2 = basis_f5c<F>(vec);
             end = clock();      
             for (auto& pol1 : s2) {
                 result.append(QString::fromStdString(to_string(pol1)) + " = 0");
                 //std::cout << to_string(pol1) << std::endl;
             }
             break;
         }
    }
    result.append("Времени затрачено:\n"+QString::number((end-begin)/(double)CLOCKS_PER_SEC));
   return result;
}


QStringList rInput(std::istream& ist, Field_FGF2 one, Ordering order, GrobAlgo algo) {
    Token_stream<Field_FGF2> ts(ist, one, order);
    std::vector<Polynom<Field_FGF2>> vec;
    QStringList result;
    clock_t begin = 0, end = 0;

   while (ts.get().kind != Kind::end) {
      //!ts.get().value.empty()
      auto pol = expr(ts, false);
      pol.sort_terms();
      //std::cout << to_string(pol) << std::endl;
      bool unique = true;
      for (auto& p : vec) {
          if (p == pol) { unique = false; break; }
      }
      if (unique) {
          vec.push_back(pol/pol.LC());
      }
   }

   if (vec.empty()) return result;
   auto beta_deg = Field_FGF2::find_beta_deg(one.prim_pol());
   auto matr = Field_FGF2::make_normal_matr(beta_deg, one.prim_pol());
   //std::cout << Field_FGF2::find_beta_deg(one.prim_pol()) << std::endl;
   result.append("B = a^{" + QString::number(beta_deg) + "}");
   Field_FGF2 a(one.prim_pol(), Field_FGF2::to_signature(beta_deg, one.prim_pol()), beta_deg);
   QString norm_basis = QString::fromStdString(to_string(a));
   for (int64_t d = 1; d < one.deg(); d++) {
       a = Field_FGF2(one.prim_pol(), Field_FGF2::to_signature(beta_deg*(1 << d), one.prim_pol()), beta_deg*(1<< d));
       norm_basis += ", " + QString::fromStdString(to_string(a));
   }
   result.append("Нормальный базис: " + norm_basis);
   auto pol = extend_vars(vec[0], matr);
   pol = polynom_to_normal_basis(Field_FGF2::inverse_matr(matr), pol);
   std::vector<Polynom<Field_ZP>> system = make_system(pol);

   switch (algo) {
        case GrobAlgo::Buchberger:
        {
            begin = clock();
            auto s2 = basis_buchberger<Field_ZP>(system);
            reduce_basis(s2);
            end = clock();
            for (auto& pol : s2) {
                result.append(QString::fromStdString(to_string(pol)) + " = 0");
                //std::cout << to_string(pol) << std::endl;
            }
            break;
        }
        case GrobAlgo::BuchbergerBetter:
        {
            begin = clock();
            auto s2 = basis_buchberger_better<Field_ZP>(system);
            reduce_basis(s2);
            end = clock();
            for (auto& pol : s2) {
                result.append(QString::fromStdString(to_string(pol)) + " = 0");
                //std::cout << to_string(pol) << std::endl;
            }
            break;
        }
        case GrobAlgo::F5:
        {
            begin = clock();
            auto s2 = basis_f5<Field_ZP>(system);
            reduce_basis(s2);
            end = clock();
            for (auto& pol : s2) {
                result.append(QString::fromStdString(to_string(pol)) + " = 0");
                //std::cout << to_string(pol) << std::endl;
            }
            break;
        }
        case GrobAlgo::F5C:
        {
            begin = clock();
            auto s2 = basis_f5c<Field_ZP>(system);
            end = clock();
            for (auto& pol : s2) {
                result.append(QString::fromStdString(to_string(pol)) + " = 0");
                //std::cout << to_string(pol) << std::endl;
            }
            break;
        }
   }
   result.append("Времени затрачено:\n"+QString::number((end-begin)/(double)CLOCKS_PER_SEC));
   return result;

//   for (auto& poly : vec) {
//       std::cout << to_string(poly) << std::endl;
//   }
}

#include <iostream>
void Interface::readSystem() {
    try {
        QString poly = ui->modInput->text();
        QString content = ui->systemInput->toPlainText();
        std::stringstream ist;
        std::stringstream polst;
        ist << content.toStdString()+"\n";
        QStringList result;

        switch (field_type)
        {
            case FieldTypes::RationalT:
                result = rInput<Rational>(ist, Rational(1), ord_type, algo);
                break;
            case FieldTypes::ComplexT:
                result = rInput<Complex>(ist, 1, ord_type, algo);
                break;
            case FieldTypes::GF_2T:
            {
                poly.replace("a", "x");
                polst << poly.toStdString() + "\n";
                Token_stream<Field_ZP> tsl(polst, {2, 1}, ord_type);
                if (tsl.get().value.empty())
                    throw std::runtime_error{"Wrong primitive polynomial"};
                auto pol = expr(tsl, false);
                uint64_t prim = 0;
                for (auto& term : pol.terms) {
                    prim |= (1 << (int)term.mon.vars["x"]);
                }
                result = rInput(ist, {prim, 1}, ord_type, algo);
                break;
            }
            case FieldTypes::ZpT:
                result = rInput<Field_ZP>(ist, {poly.toULongLong(), 1}, ord_type, algo);
                break;
        }
        writeOutput(result);
    } catch (std::exception& ex) {
        showErrorMessage(QString::fromStdString(ex.what()));
    }

}


void Interface::writeOutput(QStringList& output) {
    QString html = "<html><head></head><body>";
    QStringList left;
    QStringList rest;
    for (auto& equation: output) {
        html += "<p>";
        equation.replace("^{", "<span style=\"vertical-align: super; color: red;\">");
        equation.replace("_{", "<span style=\"vertical-align: sub; color: green;\">");
        equation.replace("}", "</span>");

        html += equation;
        html += "</p>";
    }

    html += "</body></html>";
    ui->outputLabel->setTextFormat(Qt::RichText);
    ui->outputLabel->setText(html);

    //std::cout << html.toLocal8Bit().data() << std::endl;
}
