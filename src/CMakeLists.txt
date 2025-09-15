cmake_minimum_required (VERSION 3.13.0)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")
project (Tomator1D)

file(GLOB_RECURSE EIGEN "Eigen/*.cpp" "Eigen/*.hpp")
file(GLOB_RECURSE FUNCS "Funcs/*.cpp" "Funcs/*.hpp")
file(GLOB_RECURSE VARS "Vars/*.cpp" "Vars/*.hpp" "Vars/*.c")


#Required
find_package(OpenMP REQUIRED)
add_executable(${PROJECT_NAME} Tomator1D.cpp ${EIGEN} ${FUNCS} ${VARS})
target_compile_options(${PROJECT_NAME} PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})

#If using IC
find_package(PythonLibs REQUIRED)
find_package(Python COMPONENTS Development REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${PYTHON_LIBRARIES})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
