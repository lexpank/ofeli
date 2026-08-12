cmake_minimum_required(VERSION 3.16)

project(Ofeli VERSION 1.0.8 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Qt6 6.11.1 EXACT REQUIRED COMPONENTS Core Widgets PrintSupport LinguistTools)
qt_standard_project_setup()

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(OFELI_CORE_SOURCES
    activecontour.cpp
    ac_withoutedges.cpp
    ac_withoutedges_yuv.cpp
    filters.cpp
    geodesic_ac.cpp
    hausdorff_distance.cpp
)

set(OFELI_CORE_HEADERS
    activecontour.hpp
    ac_withoutedges.hpp
    ac_withoutedges_yuv.hpp
    filters.hpp
    geodesic_ac.hpp
    hausdorff_distance.hpp
    linked_list.hpp
    linked_list.tpp
)

add_library(ofeli_core STATIC ${OFELI_CORE_SOURCES} ${OFELI_CORE_HEADERS})
target_include_directories(ofeli_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(ofeli_core PUBLIC cxx_std_17)

set(OFELI_RESOURCES)
qt_add_resources(OFELI_RESOURCES Ofeli.qrc)

qt_add_executable(Ofeli
    main.cpp
    imageviewer.cpp
    imageviewer.hpp
    pixmapwidget.cpp
    pixmapwidget.hpp
    ${OFELI_RESOURCES}
)
target_link_libraries(Ofeli PRIVATE ofeli_core Qt6::Widgets Qt6::PrintSupport)
target_compile_definitions(Ofeli PRIVATE OFELI_VERSION="${PROJECT_VERSION}")
set_target_properties(Ofeli PROPERTIES WIN32_EXECUTABLE TRUE)

if(APPLE)
    target_sources(Ofeli PRIVATE Ofeli.icns)
    set_source_files_properties(Ofeli.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
    set_target_properties(Ofeli PROPERTIES MACOSX_BUNDLE TRUE MACOSX_BUNDLE_ICON_FILE Ofeli.icns)
elseif(WIN32)
    target_sources(Ofeli PRIVATE Ofeli.rc)
endif()

qt_add_translations(Ofeli
    SOURCE_TARGETS Ofeli
    TS_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Ofeli_fr.ts"
)

include(CTest)
if(BUILD_TESTING)
    add_executable(ofeli_core_smoke tests/core_smoke.cpp)
    target_link_libraries(ofeli_core_smoke PRIVATE ofeli_core)
    add_test(NAME ofeli_core_smoke COMMAND ofeli_core_smoke)
endif()
