cmake_minimum_required(VERSION 3.16)

project(Ofeli VERSION 1.0.8 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(OFELI_BUILD_GUI "Build the Qt desktop application" ON)
option(OFELI_BUILD_CLI "Build the Qt command-line application" ON)
option(OFELI_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

if(OFELI_ENABLE_SANITIZERS AND NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    message(FATAL_ERROR "OFELI_ENABLE_SANITIZERS requires GCC or Clang")
endif()

if(OFELI_BUILD_GUI OR OFELI_BUILD_CLI)
    set(OFELI_QT_COMPONENTS Core Gui)
    if(OFELI_BUILD_GUI)
        list(APPEND OFELI_QT_COMPONENTS Widgets PrintSupport LinguistTools)
    endif()
    find_package(Qt6 6.5 REQUIRED COMPONENTS ${OFELI_QT_COMPONENTS})
    qt_standard_project_setup()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

function(ofeli_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    if(OFELI_ENABLE_SANITIZERS)
        target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
endfunction()

set(OFELI_CORE_SOURCES
    src/core/activecontour.cpp
    src/core/ac_withoutedges.cpp
    src/core/ac_withoutedges_yuv.cpp
    src/core/filters.cpp
    src/core/geodesic_ac.cpp
    src/core/hausdorff_distance.cpp
    src/core/segmentation.cpp
)

set(OFELI_CORE_HEADERS
    include/ofeli/activecontour.hpp
    include/ofeli/ac_withoutedges.hpp
    include/ofeli/ac_withoutedges_yuv.hpp
    include/ofeli/filters.hpp
    include/ofeli/geodesic_ac.hpp
    include/ofeli/hausdorff_distance.hpp
    include/ofeli/linked_list.hpp
    include/ofeli/linked_list.tpp
    include/ofeli/segmentation.hpp
)

add_library(ofeli_core STATIC ${OFELI_CORE_SOURCES} ${OFELI_CORE_HEADERS})
set_target_properties(ofeli_core PROPERTIES EXPORT_NAME Core)
target_include_directories(ofeli_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(ofeli_core PUBLIC cxx_std_17)
ofeli_enable_warnings(ofeli_core)

if(OFELI_BUILD_GUI OR OFELI_BUILD_CLI)
    add_library(ofeli_qt_image_io STATIC src/qt/qt_image_io.cpp src/qt/qt_image_io.hpp)
    target_link_libraries(ofeli_qt_image_io PUBLIC ofeli_core Qt6::Core Qt6::Gui)
    target_include_directories(ofeli_qt_image_io PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/qt>
    )
    ofeli_enable_warnings(ofeli_qt_image_io)
endif()

if(OFELI_BUILD_GUI)
    set(OFELI_RESOURCES)
    qt_add_resources(OFELI_RESOURCES apps/gui/Ofeli.qrc)

    qt_add_executable(Ofeli
        apps/gui/main.cpp
        apps/gui/imageviewer.cpp
        apps/gui/imageviewer.hpp
        apps/gui/pixmapwidget.cpp
        apps/gui/pixmapwidget.hpp
        ${OFELI_RESOURCES}
    )
    target_link_libraries(Ofeli PRIVATE ofeli_core ofeli_qt_image_io Qt6::Widgets Qt6::PrintSupport)
    target_compile_definitions(Ofeli PRIVATE OFELI_VERSION="${PROJECT_VERSION}")
    set_target_properties(Ofeli PROPERTIES WIN32_EXECUTABLE TRUE)
    ofeli_enable_warnings(Ofeli)

    if(APPLE)
        target_sources(Ofeli PRIVATE apps/gui/Ofeli.icns)
        set_source_files_properties(apps/gui/Ofeli.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
        set_target_properties(Ofeli PROPERTIES MACOSX_BUNDLE TRUE MACOSX_BUNDLE_ICON_FILE Ofeli.icns)
    elseif(WIN32)
        target_sources(Ofeli PRIVATE apps/gui/Ofeli.rc)
    endif()

    qt_add_translations(Ofeli
        SOURCE_TARGETS Ofeli
        TS_FILES "${CMAKE_CURRENT_SOURCE_DIR}/apps/gui/Ofeli_fr.ts"
    )
endif()

if(OFELI_BUILD_CLI)
    qt_add_executable(ofeli_cli apps/cli/cli_main.cpp)
    target_link_libraries(ofeli_cli PRIVATE ofeli_core ofeli_qt_image_io Qt6::Core Qt6::Gui)
    target_compile_definitions(ofeli_cli PRIVATE OFELI_VERSION="${PROJECT_VERSION}")
    set_target_properties(ofeli_cli PROPERTIES OUTPUT_NAME ofeli-cli MACOSX_BUNDLE FALSE)
    ofeli_enable_warnings(ofeli_cli)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(ofeli_core_smoke tests/core_smoke.cpp)
    target_link_libraries(ofeli_core_smoke PRIVATE ofeli_core)
    ofeli_enable_warnings(ofeli_core_smoke)
    add_test(NAME ofeli_core_smoke COMMAND ofeli_core_smoke)

    add_executable(ofeli_segmentation_smoke tests/segmentation_smoke.cpp)
    target_link_libraries(ofeli_segmentation_smoke PRIVATE ofeli_core)
    ofeli_enable_warnings(ofeli_segmentation_smoke)
    add_test(NAME ofeli_segmentation_smoke COMMAND ofeli_segmentation_smoke)

    if(OFELI_BUILD_GUI OR OFELI_BUILD_CLI)
        add_executable(ofeli_qt_image_io_smoke tests/qt_image_io_smoke.cpp)
        target_link_libraries(ofeli_qt_image_io_smoke PRIVATE ofeli_qt_image_io)
        ofeli_enable_warnings(ofeli_qt_image_io_smoke)
        add_test(NAME ofeli_qt_image_io_smoke COMMAND ofeli_qt_image_io_smoke)
    endif()

    if(OFELI_BUILD_CLI)
        add_test(NAME ofeli_cli_smoke
                 COMMAND ${CMAKE_COMMAND}
                         -DCLI=$<TARGET_FILE:ofeli_cli>
                         -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/cli-smoke
                         -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_smoke.cmake)

        if(NOT OFELI_ENABLE_SANITIZERS)
            add_test(NAME ofeli_install_consumer
                     COMMAND ${CMAKE_COMMAND}
                             -DOFELI_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                             -DOFELI_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
                             -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/install_consumer.cmake)
        endif()
    endif()
endif()

install(TARGETS ofeli_core
        EXPORT OfeliTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(OFELI_BUILD_CLI)
    install(TARGETS ofeli_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(OFELI_BUILD_GUI)
    install(TARGETS Ofeli
            BUNDLE DESTINATION .
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

set(OFELI_PACKAGE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/Ofeli")
install(EXPORT OfeliTargets
        FILE OfeliTargets.cmake
        NAMESPACE Ofeli::
        DESTINATION ${OFELI_PACKAGE_INSTALL_DIR})

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/OfeliConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OfeliConfig.cmake"
    INSTALL_DESTINATION ${OFELI_PACKAGE_INSTALL_DIR})
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/OfeliConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/OfeliConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/OfeliConfigVersion.cmake"
        DESTINATION ${OFELI_PACKAGE_INSTALL_DIR})
