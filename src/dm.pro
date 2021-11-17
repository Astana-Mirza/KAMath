QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    interface.cpp \
    polynom/monom.cpp \
    fields/galois_fast.cpp \
    nonlinear_2n.cpp \
    fields/extra_fields.cpp \
    parser/field_reader.cpp \
    fields/modulo_prime.cpp

HEADERS += \
    interface.h \
    polynom/monom.h \
    polynom/polynom.h \
    fields/galois_fast.h \
    fields/modulo_prime.h \
    nonlinear_2n.h \
    fields/extra_fields.h \
    parser/field_reader.h \
    polynom/monom.h fields/galois_fast.h fields/extra_fields.h polynom/polynom.h grobner_alg/buchberger.h grobner_alg/faugere.h parser/field_reader.h parser/parser.h

FORMS += \
    interface.ui

OBJECTS_DIR += build
MOC_DIR += build
UI_DIR += build

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
