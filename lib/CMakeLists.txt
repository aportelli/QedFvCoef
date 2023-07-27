add_library(qedfv SHARED)
target_sources(qedfv
  PRIVATE qedfv/coef.cpp qedfv/latticesum.cpp
  PUBLIC 
    FILE_SET public_headers
    TYPE HEADERS
    FILES
      qedfv/coef.hpp
      qedfv/global.hpp
      qedfv/latticesum.hpp
      qedfv/vectorutils.hpp)
target_link_libraries(qedfv PUBLIC OpenMP::OpenMP_CXX GSL::gsl GSL::gslcblas)
set_target_properties(qedfv PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION 1)
target_include_directories(
  qedfv
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS qedfv
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  FILE_SET public_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qedfv)
