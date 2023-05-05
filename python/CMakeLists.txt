FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.10.4
)

FetchContent_MakeAvailable(pybind11)
pybind11_add_module(pyqedfv qedfv.cpp)
target_include_directories(pyqedfv PRIVATE ${PROJECT_SOURCE_DIR}/lib)
target_link_libraries(pyqedfv PUBLIC qedfv)
target_compile_definitions(pyqedfv PUBLIC VERSION_INFO=${CMAKE_PROJECT_VERSION})
