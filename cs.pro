#-------------------------------------------------
# libcs
# author: pdobrowolski
#-------------------------------------------------
QT       -= core gui qt
TARGET    = cs
TEMPLATE  = lib
DESTDIR   = ../bin

DEFINES  += LIB_LIBRARY

GXX_WARNING_FLAGS = \
#   -Werror                      \ # Make all warnings into errors.
    -Wall                        \ # Generate all warnings
    -Wextra                      \ # Generate even more extra warnings
#   -pedantic                    \ # Accept only pedantic code
#   -Weffc++                     \ # Accept only effective C++ code
    -Wwrite-strings              \ # Do not accept writing to contant string memory
    -Winit-self                  \ # Do not accept initializing variable with itself
    -Wcast-align                 \ # Do not accept misaligning with casting
    -Wcast-qual                  \ # Do not accept removing qualifiers with casting
#   -Wold-style-cast             \ # Do not accept old style casting
    -Wpointer-arith              \ # Warn about void pointer arthimetic
    -Wstrict-aliasing            \ # Ensure strict aliasing
    -Wuninitialized              \ # Do not accept uninitialized variables
#   -Wmissing-declarations       \ # Warn about global and non-accesible functions
    -Woverloaded-virtual         \ # Warn about incidental overiding non-virtual base methods
    -Wnon-virtual-dtor           \ # Warn about non-virtual destructor
#   -Wctor-dtor-privacy          \ # Warn about useless and non-constructible classes
#   -Wlong-long                  \ # Do not allow using long long
    -Wunreachable-code           \ # Warn about unreachable code
#   -Wfloat-equal                \ # Do not accept comparing floating points with equal operator
    -Wabi                        \ # Warn about possible ABI problems
#   -Wswitch-enum                \ # Check switch enumeration
    -Wformat=2                   \ # Check printf formatting
#   -Wundef                      \ # Warn if an undefined identifier is evaluated in an @if directive.
#   -Wshadow                     \ # Warn whenever a local variable shadows another local variable, parameter or global variable or whenever a built-in function is shadowed
#   -Wconversion                 \ # Warn for implicit conversions that may alter a value
    -Wlogical-op                 \ # Warn about suspicious uses of logical operators in expressions
#   -Waggregate-return           \ # Warn if any functions that return structures or unions are defined or called.
    -Wmissing-field-initializers \ # Warn if a structure's initializer has some fields missing.
    -Wredundant-decls            \ # Warn if anything is declared more than once in the same scope, even in cases where multiple declaration is valid and changes nothing.
    -Wmissing-include-dirs       \ # Warn if a user-supplied include directory does not exist.
#   -Wswitch-default             \ # Warn whenever a switch statement does not have a default case.
    -Wsync-nand                  \ # Warn when __sync_fetch_and_nand and __sync_nand_and_fetch built-in functions are used. These functions changed semantics in GCC 4.4.
    -Wunused                     \ # All the above -Wunused options combined.
#   -Wstrict-overflow=5          \ # Also warn about cases where the compiler reduces the magnitude of a constant involved in a comparison.
#   -Wunsafe-loop-optimizations  \ # Warn if the loop cannot be optimized because the compiler could not assume anything on the bounds of the loop indices. With -funsafe-loop-optimizations warn if the compiler made such assumptions.
    -Wmissing-format-attribute   \ # Warn about function pointers which might be candidates for format attributes.
#   -Wpadded                     \ # Warn if padding is included in a structure, either to align an element of the structure or to align the whole structure.
#   -Winline                     \ # Warn if a function can not be inlined and it was declared as inline.
    -Wdisabled-optimization      \ # Warn if a requested optimization pass is disabled.

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
}
linux*:{
    QMAKE_CXXFLAGS += $$GXX_WARNING_FLAGS

    # CGAL needs that
    QMAKE_CXXFLAGS += -frounding-math
}

#
# libs
#
LIBS += -lLiDIA
LIBS += -lqi
LIBS += -lCGAL
LIBS += -lmpfi
LIBS += -lmpfr
LIBS += -lgmp
LIBS += -lrealroot
LIBS += -llog4cxx

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

# sources
SOURCES += \
    lib/Uspensky_Hanrot.c \
    lib/Spin_cell_graph_3.cpp \
    lib/Spin_3.cpp \
    lib/Scene_loader.cpp \
    lib/Random.cpp \
    lib/Math_utils.cpp \
    lib/Hom_root.cpp \
    lib/Coordinate.cpp \
    lib/Contfrac.cpp \
    lib/Benchmark.cpp \
    lib/Build.cpp \
    lib/Config.cpp

# headers
HEADERS += \
    3rdparty/libqi.h \
    lib/Uspensky_Hanrot.h \
    lib/Uniform_random_spin_3.ipp \
    lib/Spin_straight_sample_generator_3.ipp \
    lib/Spin_quadric_tree_generator.ipp \
    lib/Spin_quadric_tree_3.ipp \
    lib/Spin_quadric_mesh_3.ipp \
    lib/Spin_quadric_3.ipp \
    lib/Spin_qsip_point.ipp \
    lib/Spin_qsip_mesh_3.ipp \
    lib/Spin_qsip_3.ipp \
    lib/Spin_qsic_3.ipp \
    lib/Spin_QC_sample_generator_3.ipp \
    lib/Spin_kernel_3.ipp \
    lib/Spin_configuration_space_3.ipp \
    lib/Spin_cell_graph_3.ipp \
    lib/Spin_3.ipp \
    lib/Smooth_triangle_3.ipp \
    lib/Random_spin_circle_3.ipp \
    lib/Quartic.ipp \
    lib/Quadratic.ipp \
    lib/Predicate_tt_3.ipp \
    lib/Predicate_s_3.ipp \
    lib/Predicate_h_3.ipp \
    lib/Predicate_g_3.ipp \
    lib/Predicate_bb_3.ipp \
    lib/Polynomial_sign_at.ipp \
    lib/Polynomial_interval.ipp \
    lib/Math_utils.ipp \
    lib/Isolate_polynomial_roots.ipp \
    lib/Hom_root.ipp \
    lib/Cubic.ipp \
    lib/Convert_polynomial.ipp \
    lib/Cell_n.ipp \
    lib/Bigint.ipp \
    lib/Benchmark.ipp \
    lib/Uniform_random_spin_3.h \
    lib/Spin_straight_sample_generator_3.h \
    lib/Spin_quadric_tree_generator.h \
    lib/Spin_quadric_tree_3.h \
    lib/Spin_quadric_mesh_3.h \
    lib/Spin_quadric_3.h \
    lib/Spin_qsip_point.h \
    lib/Spin_qsip_mesh_3.h \
    lib/Spin_qsip_3.h \
    lib/Spin_qsic_mesh_3.h \
    lib/Spin_qsic_3.h \
    lib/Spin_QC_sample_generator_3.h \
    lib/Spin_kernel_3.h \
    lib/Spin_configuration_space_3.h \
    lib/Spin_cell_graph_3.h \
    lib/Spin_3.h \
    lib/Sphere_triangle_intersection_3.h \
    lib/Smooth_triangle_3.h \
    lib/Random_spin_circle_3.h \
    lib/Random.h \
    lib/Quartic.h \
    lib/Quadratic.h \
    lib/Predicate_tt_3.h \
    lib/Predicate_s_3.h \
    lib/Predicate_h_3.h \
    lib/Predicate_g_3.h \
    lib/Predicate_bb_3.h \
    lib/Polynomial_sign_at.h \
    lib/Polynomial_interval.h \
    lib/Math_utils.h \
    lib/Isolate_polynomial_roots.h \
    lib/Hom_root.h \
    lib/Extended.h \
    lib/Cubic.h \
    lib/Coordinate.h \
    lib/Convert_polynomial.h \
    lib/Contfrac.h \
    lib/Cell_n.h \
    lib/Bigint.h \
    lib/Benchmark.h \
    lib/Sphere_triangle_intersection_3.ipp \
    lib/Spin_qsic_mesh_3.ipp \
    lib/Index_tree_n.ipp \
    lib/Index_tree_n.h \
    lib/Predicate_list_generator.h \
    lib/Predicate_tt_3_list_generator.ipp \
    lib/Predicate_tt_3_list_generator.h \
    lib/Predicate_bb_3_list_generator.ipp \
    lib/Predicate_bb_3_list_generator.h \
    lib/Ball_3.ipp \
    lib/Ball_3.h \
    lib/Loader_sphere_tree.ipp \
    lib/Loader_sphere_tree.h \
    lib/Loader_scene.h \
    lib/Loader_scene.ipp \
    lib/Build.h \
    lib/Spin_raster_graph_3.ipp \
    lib/Spin_raster_graph_3.h \
    lib/Voxel_3.h \
    lib/Voxel_3.ipp \
    lib/Spin_exact_graph_3.ipp \
    lib/Spin_exact_graph_3.h \
    lib/Config.h
