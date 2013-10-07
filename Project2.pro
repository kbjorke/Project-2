TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    jacobi_algorithm.cpp \
    harmonic_oscillator_3d.cpp \
    output_functions.cpp \
    lib.cpp

LIBS += -larmadillo -lblas -llapack -lunittest++ -lrt

HEADERS += \
    jacobi_algorithm.h \
    harmonic_oscillator_3d.h \
    output_functions.h \
    lib.h
