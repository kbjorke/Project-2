TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    jacobi_algorithm.cpp

LIBS += -larmadillo -lblas -llapack -lunittest++

HEADERS += \
    jacobi_algorithm.h
