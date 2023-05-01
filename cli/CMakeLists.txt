file(DOWNLOAD https://raw.githubusercontent.com/aportelli/OptParser/main/OptParser.hpp 
  ${CMAKE_BINARY_DIR}/deps/include/OptParser.hpp)
add_library(optparser INTERFACE)
target_include_directories(optparser INTERFACE ${CMAKE_BINARY_DIR}/deps/include)

add_executable(qed-fv-coef coef.cpp)
target_include_directories(qed-fv-coef PRIVATE ${PROJECT_SOURCE_DIR}/lib)
target_link_libraries(qed-fv-coef qedfvcoef OpenMP::OpenMP_CXX optparser)
install(TARGETS qed-fv-coef DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(qed-fv-coef-spherescan sphere.cpp)
target_include_directories(qed-fv-coef-spherescan PRIVATE ${PROJECT_SOURCE_DIR}/lib)
target_link_libraries(qed-fv-coef-spherescan qedfvcoef OpenMP::OpenMP_CXX optparser)
install(TARGETS qed-fv-coef-spherescan DESTINATION ${CMAKE_INSTALL_BINDIR})
