TEMPLATE = app

VERSION = 1.0.7

TARGET = Ofeli

QT += widgets
!isEmpty(QT.printsupport.name): QT += printsupport

simulator: warning(This example might not fully work on Simulator platform)


HEADERS += activecontour.hpp \
    filters.hpp \
    hausdorff_distance.hpp \
    ac_withoutedges.hpp \
    geodesic_ac.hpp \
    ac_withoutedges_yuv.hpp \
    linked_list.hpp \
    imageviewer.hpp \
    pixmapwidget.hpp \

SOURCES += activecontour.cpp \
    hausdorff_distance.cpp \
    filters.cpp \
    ac_withoutedges.cpp \
    ac_withoutedges_yuv.cpp \
    geodesic_ac.cpp \
    main.cpp \
    imageviewer.cpp \
    pixmapwidget.cpp \
    linked_list.tpp

# statically linked version of Boost library
DEFINES = BOOST_THREAD_USE_LIB

INCLUDEPATH += $$PWD/boost_1_52_0 \

TRANSLATIONS = $$PWD/Ofeli_fr.ts
CODECFORTR = UTF-8
CODECFORSRC = UTF-8

RESOURCES = $$PWD/Ofeli.qrc

macx : ICON = $$PWD/Ofeli.icns

win32 : RC_FILE = $$PWD/Ofeli.rc

QMAKE_CXXFLAGS_RELEASE += -O3
