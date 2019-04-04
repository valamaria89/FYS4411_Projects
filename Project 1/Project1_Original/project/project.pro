TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    exercises.cpp \
    methods.cpp \
    methods_interaction.cpp
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -llapack -lblas -larmadillo
HEADERS += \
    exercises.h \
    ../../../Github_FYS4150/Computational-Physics-FYS3150/Projects/Project4/mag_sys_no_MPI/methods.h \
    methods.h \
    matrix.h \
    methods_interaction.h
