#*****************************************************************//
#    Albany 3.0:  Copyright 2016 Sandia Corporation               //
#    This Software is released under the BSD license detailed     //
#    in the file "license.txt" in the top-level Albany directory  //
#*****************************************************************//

SET(SOURCES
  Albany_DealiiDiscretization.cpp
  Albany_DealiiMeshStruct.cpp
  )
SET(HEADERS
  Albany_DealiiDiscretization.hpp
  Albany_DealiiMeshStruct.hpp
  )

INCLUDE_DIRECTORIES (
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${Albany_SOURCE_DIR}/src
)
add_library(albanyDealii ${Albany_LIBRARY_TYPE} ${SOURCES} ${HEADERS})
target_link_libraries(albanyDealii ${Trilinos_LIBRARIES})

#set_target_properties(albanyDealii PROPERTIES PUBLIC_HEADER "${HEADERS}")

IF (INSTALL_ALBANY)
  install(TARGETS albanyDealii EXPORT albany-export
    LIBRARY DESTINATION "${LIB_INSTALL_DIR}/"
    ARCHIVE DESTINATION "${LIB_INSTALL_DIR}/"
    PUBLIC_HEADER DESTINATION "${INCLUDE_INSTALL_DIR}")
ENDIF()
