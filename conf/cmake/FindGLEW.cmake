# - Try to find GLEW
# Once done this will define
#
#  GLEW_FOUND - system has GLEW
#  GLEW_INCLUDE_DIR - the GLEW include directory
#  GLEW_LIBRARIES - Link these to use GLEW
#  GLEW_DEFINITIONS - Compiler switches required for using GLEW
#

FIND_PATH(GLEW_INCLUDE_DIR NAMES GL/glew.h
  PATHS
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/libraries/glew/include
  ENV CPATH
  /usr/include
  /usr/local/include
  /opt/local/include
  NO_DEFAULT_PATH
)

FIND_LIBRARY(GLEW_LIBRARIES NAMES GLEW GLEW32
  PATHS
  ${PROJECT_BINARY_DIR}/lib
  ${PROJECT_SOURCE_DIR}/lib
  ${PROJECT_SOURCE_DIR}/libraries/glew/lib
  ENV LD_LIBRARY_PATH
  ENV LIBRARY_PATH
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  NO_DEFAULT_PATH
)

IF(GLEW_INCLUDE_DIR AND GLEW_LIBRARIES)
   SET(GLEW_FOUND TRUE)
ENDIF(GLEW_INCLUDE_DIR AND GLEW_LIBRARIES)

# show the GLEW_INCLUDE_DIR and GLEW_LIBRARIES variables only in the advanced view
IF(GLEW_FOUND)
  MARK_AS_ADVANCED(GLEW_INCLUDE_DIR GLEW_LIBRARIES )
ENDIF(GLEW_FOUND)
