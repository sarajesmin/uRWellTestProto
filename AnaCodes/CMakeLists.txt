set(INCLUDES
  ${ROOT_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/build/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4)
include_directories(${INCLUDES})
link_directories(${CMAKE_SOURCE_DIR}/build/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)


################################  SkimZeroSuppression  #############################################
add_executable(SkimZeroSuppression.exe SkimZeroSuppression.cc)
target_link_libraries(SkimZeroSuppression.exe Hipo4 uRWellLib ${ROOT_LIBRARIES} )
#################################################################################################

################################  CheckDecoding.cc  #############################################
add_executable(CheckDecoding.exe CheckDecoding.cc)
target_link_libraries(CheckDecoding.exe Hipo4 uRWellLib ${ROOT_LIBRARIES} )
#################################################################################################

################################  AnaClustering.cc  #############################################
add_executable(AnaClustering.exe AnaClustering.cc)
target_link_libraries(AnaClustering.exe Hipo4 uRWellLib ${ROOT_LIBRARIES} )
#################################################################################################


################################ Installing Executable targets #############################
install( TARGETS SkimZeroSuppression.exe CheckDecoding.exe AnaClustering.exe
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )

################################ Installing Non target files ###############################
install( FILES DrawPedestals.cc DrawPlotsWithClustering.cc
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


################################ Installing needed directores ##############################
install( DIRECTORY
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/PedFiles)
install( DIRECTORY
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/Figs )


