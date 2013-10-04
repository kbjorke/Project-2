TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    jacobi_algorithm.cpp \
    harmonic_oscillator_3d.cpp \
    output_functions.cpp

LIBS += -larmadillo -lblas -llapack -lunittest++

HEADERS += \
    jacobi_algorithm.h \
    harmonic_oscillator_3d.h \
    output_functions.h
