#-------------------------------------------------
# test
# author: pdobrowolski
#-------------------------------------------------
TEMPLATE = app
TARGET   = test
CONFIG  += console
CONFIG  -= qt

DESTDIR = ../bin

#
# compiler
#
ARCH = "../3rdparty"

# architecture
contains(QMAKE_HOST.arch, x86_64):{
    ARCH = "$$ARCH/amd64"
} else {
    ARCH = "$$ARCH/i386"
}

# system
win32*:{
    ARCH = "$$ARCH/windows"
}
linux*:{
    ARCH = "$$ARCH/linux"
}

message(architecture: $$ARCH)

# lidia
LIDIA = "$$ARCH/lidia-2.3.0"
INCLUDEPATH += "$$LIDIA/include"
LIBS += -L"$$LIDIA/lib"

# libqi
LIBQI = "$$ARCH/libqi-0.9.33"
INCLUDEPATH += "$$LIBQI/include"
LIBS += -L"$$LIBQI/lib"

# realroot
REALROOT = "$$ARCH/realroot"
INCLUDEPATH += "$$REALROOT/include"
LIBS += -L"$$REALROOT/lib"

# optimization
CONFIG(release, debug|release) {
    QMAKE_CXXFLAGS += -DCGAL_NDEBUG
    QMAKE_CXXFLAGS -= -D_DEBUG
    QMAKE_CXXFLAGS += -DNDEBUG
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3 -march=core2
}

# parallel
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp

# CS
INCLUDEPATH += "../include"
LIBS += -L"../bin"

win32-g++:{
    QMAKE_CXXFLAGS += $$GXX_WARNING_FLAGS

    # boost needs that
    QMAKE_CXXFLAGS += -DBOOST_THREAD_USE_LIB

    # CGAL needs that
    QMAKE_CXXFLAGS += -frounding-math

    # boost
    BOOST = "$$ARCH/boost-1.49"
    INCLUDEPATH += "$$BOOST/include"
    LIBS += -L"$$BOOST/lib"

    # GMP
    GMP = "$$ARCH/gmp-5.0.5"
    INCLUDEPATH += "$$GMP/include"
    LIBS += -L"$$GMP/lib"

    # MPFR
    MPFR = "$$ARCH/mpfr-3.1.0"
    INCLUDEPATH += "$$MPFR/include"
    LIBS += -L"$$MPFR/lib"

    # MPFI
    MPFI = "$$ARCH/mpfi-1.5.1"
    INCLUDEPATH += "$$MPFI/include"
    LIBS += -L"$$MPFI/lib"

    # CGAL
    CGAL = "$$ARCH/CGAL-4.0"
    INCLUDEPATH += "$$CGAL/include"
    LIBS += -L"$$CGAL/lib"

    # APR
    APR = "$$ARCH/apr-1.4.6"
    INCLUDEPATH += "$$APR/include"
    LIBS += -L"$$APR/lib"

    # APRUTIL
    APRUTIL = "$$ARCH/apr-util-1.4.1"
    INCLUDEPATH += "$$APRUTIL/include"
    LIBS += -L"$$APRUTIL/lib"

    # LOG4CXX
    LOG4CXX = "$$ARCH/log4cxx-0.10.0"
    INCLUDEPATH += "$$LOG4CXX/include"
    LIBS += -L"$$LOG4CXX/lib"

    # GLEW
    GLEW = "$$ARCH/glew-1.7.0"
    INCLUDEPATH += "$$GLEW/include"
    LIBS += -L"$$GLEW/lib"

    # GLE
    GLE = "$$ARCH/gle-3.1.0"
    INCLUDEPATH += "$$GLE/include"
    LIBS += -L"$$GLE/lib"
}
linux*:{
    QMAKE_CXXFLAGS += $$GXX_WARNING_FLAGS

    # CGAL needs that
    QMAKE_CXXFLAGS += -frounding-math
}

#
# libs
#
#LIBS += -lqi
#LIBS += -lLiDIA
#LIBS += -lmpfi
LIBS += -lmpfr
#LIBS += -lgmp
LIBS += -lCGAL
#LIBS += -lrealroot
LIBS += -llog4cxx
LIBS += -lcs

win32-g++:{
    LIBS += -lboost_thread-mgw46-mt-1_49
    LIBS += -lapr-1
    LIBS += -laprutil-1
    LIBS += -lws2_32
    LIBS += -lexpat
    LIBS += -liconv
    LIBS += -lmswsock
}
linux*:{
    LIBS += -lboost_thread
}

# libcs
INCLUDEPATH += ../_deprecated/libcs/include
INCLUDEPATH += ..

# project path
INCLUDEPATH += src

# source
SOURCES += main.cpp \
#   test_various.cpp \
#   test_exact_perf.cpp \
#   test_cell_perf.cpp \
    test_inexact_perf.cpp

HEADERS += \
#   test_various.h \
#   test_exact_perf.h \
#   test_cell_perf.h \
    test_inexact_perf.h
