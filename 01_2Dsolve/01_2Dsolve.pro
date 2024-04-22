TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

msvc{
    QMAKE_CFLAGS += \utf-8
    QMAKE_CXXFLAGS += \utf-8
}

SOURCES += \
        main.cpp \
        mashdata.cpp \
        solution.cpp

HEADERS += \
    mashdata.h \
    solution.h

DISTFILES += \
    beam1.msh \
    test5.msh

RESOURCES += \
    mash.qrc


