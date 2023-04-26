add_library(qedfvcoef SHARED qedfvcoef.c)
set_target_properties(qedfvcoef PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(qedfvcoef PROPERTIES SOVERSION 1)
set_target_properties(qedfvcoef PROPERTIES PUBLIC_HEADER qedfvcoef.h)
target_link_libraries(qedfvcoef OpenMP::OpenMP_C GSL::gsl GSL::gslcblas)
install(TARGETS qedfvcoef
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
