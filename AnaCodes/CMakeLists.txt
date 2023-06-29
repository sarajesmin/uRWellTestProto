set(INCLUDES
  ${ROOT_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/build/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4)
include_directories(${INCLUDES})
link_directories(${CMAKE_SOURCE_DIR}/build/lib)
add_library(hipo4lib STATIC IMPORTED)
set_property( TARGET hipo4lib PROPERTY
  IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/libHipo4.a")
set(CMAKE_VERBOSE_MAKEFILE ON)


################################  SkimZeroSuppression  #############################################
message("includes is  " ${ROOT_INCLUDE_DIR})
add_executable(SkimZeroSuppression.exe SkimZeroSuppression.cc)
target_link_libraries(SkimZeroSuppression.exe Hipo4 -lpthread -lz -lm -ldl -lssl -lcrypto uRWellLib ${ROOT_LIBRARIES} )
target_include_directories(SkimZeroSuppression.exe PUBLIC ${INCLUDES} )
install( TARGETS SkimZeroSuppression.exe
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
install( TARGETS uRWellLib DESTINATION ${CMAKE_INSTALL_LIBDIR} )
#################################################################################################
