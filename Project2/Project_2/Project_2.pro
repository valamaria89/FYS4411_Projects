TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    system.cpp \
    rbm.cpp \
    gradientdescent.cpp

HEADERS += \
    matrix.h \
    system.h \
    system.h \
    rbm.h \
    gradientdescent.h
