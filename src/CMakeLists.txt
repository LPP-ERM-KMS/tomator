cmake_minimum_required (VERSION 3.13.0)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -Wall -Wextra")
project (Tomator1D)

file(GLOB_RECURSE FUNCS "Funcs/*.cpp" "Funcs/*.hpp")
file(GLOB_RECURSE VARS "Vars/*.cpp" "Vars/*.hpp" "Vars/*.c")

#Required
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
message(STATUS "Found Eigen version: ${Eigen3_VERSION}")
add_executable(${PROJECT_NAME} Tomator1D.cpp ${FUNCS} ${VARS})
target_compile_options(${PROJECT_NAME} PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
