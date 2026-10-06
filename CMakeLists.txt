cmake_minimum_required(VERSION 3.14)
project(lmeb VERSION 0.1.0
        DESCRIPTION "L-band Microwave Emission of the Biosphere (L-MEB) model, Wigneron et al. 2007"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(LMEB_BUILD_APPS  "Build command-line tools" ON)
option(LMEB_BUILD_TESTS "Build unit tests" ON)

add_library(lmeb
  src/dielectric.cpp
  src/soil.cpp
  src/vegetation.cpp
  src/atmosphere.cpp
  src/model.cpp
  src/presets.cpp
  src/retrieval.cpp
)
add_library(lmeb::lmeb ALIAS lmeb)
target_include_directories(lmeb PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_options(lmeb PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>)

if(LMEB_BUILD_APPS)
  add_executable(lmeb_forward apps/lmeb_forward.cpp apps/cli.hpp)
  target_link_libraries(lmeb_forward PRIVATE lmeb)

  add_executable(lmeb_retrieve apps/lmeb_retrieve.cpp apps/cli.hpp)
  target_link_libraries(lmeb_retrieve PRIVATE lmeb)
endif()

if(LMEB_BUILD_TESTS)
  enable_testing()
  add_executable(lmeb_tests tests/test_lmeb.cpp)
  target_link_libraries(lmeb_tests PRIVATE lmeb)
  add_test(NAME lmeb_tests COMMAND lmeb_tests)
endif()

include(GNUInstallDirs)
install(TARGETS lmeb EXPORT lmebTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/lmeb DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT lmebTargets NAMESPACE lmeb:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lmeb)
