cmake_minimum_required(VERSION 3.16)

project(ofeli_install_consumer LANGUAGES CXX)

find_package(Ofeli CONFIG REQUIRED)

add_executable(ofeli_install_consumer main.cpp)
target_link_libraries(ofeli_install_consumer PRIVATE Ofeli::Core)
