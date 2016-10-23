cmake_minimum_required (VERSION 2.6)
project (fectest)
#set(CMAKE_VERBOSE_MAKEFILE ON)
add_definitions( -O0 -g -rdynamic -fno-inline -Wstrict-prototypes -std=gnu99 -Wall)
#add_definitions( -Os -Wstrict-prototypes -Wstrict-prototypes -std=gnu99 -Wall)

set(CMAKE_SYSTEM_NAME Linux)

SET(CMAKE_BUILD_TYPE release)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(fec_test src/fec_test.c)
target_link_libraries(fec_test)

add_executable(fec_test3 src/fec_test3.c)
target_link_libraries(fec_test3)
