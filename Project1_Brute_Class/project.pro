TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    methods.cpp \
    position.cpp \
    system.cpp
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -llapack -lblas -larmadillo
HEADERS += \
    methods.h \
    matrix.h \
    position.h \
    system.h
