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

INCLUDEPATH += include extern/libkam/include
unix:LIBS   += -L$${PWD}/extern/libkam/build -lkam

make_libkam.target   = $${PWD}/extern/libkam/build/libkam.so
make_libkam.commands = cd $${PWD}/extern/libkam && \
		       make && \
		       cd ../../
PRE_TARGETDEPS      += $${PWD}/extern/libkam/build/libkam.so
QMAKE_EXTRA_TARGETS += make_libkam


SOURCES += \
     src/main.cpp \
     src/interface.cpp \

HEADERS += \
     include/interface.h \
     include/nonlinear_2n.h \
     include/parser/set_io.h \
     include/parser/parser.h

FORMS += \
     interface.ui

OBJECTS_DIR += build
MOC_DIR += build
UI_DIR += build

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
