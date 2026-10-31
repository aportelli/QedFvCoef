FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v3.0.2
)
if(EXISTS "${CMAKE_SOURCE_DIR}/.venv/bin/python")
    set(Python3_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python" CACHE FILEPATH "Python from local venv")
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
message(STATUS "Using Python3 executable: ${Python3_EXECUTABLE}")
FetchContent_MakeAvailable(pybind11)
pybind11_add_module(pyqedfv qedfv.cpp)
target_include_directories(pyqedfv PRIVATE ${PROJECT_SOURCE_DIR}/lib)
target_link_libraries(pyqedfv PUBLIC qedfv)
target_compile_definitions(pyqedfv PUBLIC VERSION_INFO=${CMAKE_PROJECT_VERSION})
