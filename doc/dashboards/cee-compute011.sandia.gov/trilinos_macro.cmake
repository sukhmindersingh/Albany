macro(do_trilinos CONFIGURE_OPTIONS BTYPE ILOC)

  message ("ctest state: BUILD_${BTYPE}")

  #
  # Configure the Trilinos build
  #


# Clean up build area
  IF (CLEAN_BUILD)
    IF(EXISTS "${CTEST_BINARY_DIRECTORY}/${BTYPE}" )
      FILE(REMOVE_RECURSE "${CTEST_BINARY_DIRECTORY}/${BTYPE}")
    ENDIF()
  ENDIF()

  if (NOT EXISTS "${CTEST_BINARY_DIRECTORY}/${BTYPE}")
    file (MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/${BTYPE})
  endif (NOT EXISTS "${CTEST_BINARY_DIRECTORY}/${BTYPE}")

  CTEST_CONFIGURE(
    BUILD "${CTEST_BINARY_DIRECTORY}/${BTYPE}"
    SOURCE "${CTEST_SOURCE_DIRECTORY}/Trilinos"
    OPTIONS "${CONFIGURE_OPTIONS}"
    RETURN_VALUE HAD_ERROR
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Configure
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message ("Cannot submit Trilinos configure results!")
    endif (S_HAD_ERROR)
  endif (CTEST_DO_SUBMIT)

  if (HAD_ERROR)
# No sense in going on if Trilinos will not config!
    message (FATAL_ERROR "Cannot configure Trilinos build!")
  endif (HAD_ERROR)

  #
  # Trilinos
  #
  # Build the rest of Trilinos and install everything
  #


  #set (CTEST_BUILD_TARGET all)
  set (CTEST_BUILD_TARGET install)

# Clean up Install area
  IF (CLEAN_BUILD)
    IF(EXISTS "${ILOC}" )
      FILE(REMOVE_RECURSE "${ILOC}")
    ENDIF()
  ENDIF()

  MESSAGE("\nBuilding target: '${CTEST_BUILD_TARGET}' ...\n")

  CTEST_BUILD(
    BUILD "${CTEST_BINARY_DIRECTORY}/${BTYPE}"
    RETURN_VALUE  HAD_ERROR
    NUMBER_ERRORS  BUILD_LIBS_NUM_ERRORS
    APPEND
    )

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Build
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message ("Cannot submit Trilinos build results!")
    endif (S_HAD_ERROR)

  endif (CTEST_DO_SUBMIT)

  if (HAD_ERROR)
# No sense in going on if Trilinos will not build!
    message (FATAL_ERROR "Cannot build Trilinos!")
  endif (HAD_ERROR)

  if (BUILD_LIBS_NUM_ERRORS GREATER 0)
# No sense in going on if Trilinos will not build!
    message (FATAL_ERROR "Encountered build errors in Trilinos build. Exiting!")
  endif (BUILD_LIBS_NUM_ERRORS GREATER 0)


# Run Trilinos tests 

  set (CTEST_TEST_TIMEOUT 600)
  CTEST_TEST(
    BUILD "${CTEST_BINARY_DIRECTORY}/${BTYPE}"
    #              PARALLEL_LEVEL "${CTEST_PARALLEL_LEVEL}"
    #              INCLUDE_LABEL "^${TRIBITS_PACKAGE}$"
    #NUMBER_FAILED  TEST_NUM_FAILED
    RETURN_VALUE  HAD_ERROR
    )

  if (HAD_ERROR)
   message("Some Trilinos tests failed.")
  endif (HAD_ERROR)

  if (CTEST_DO_SUBMIT)
    ctest_submit (PARTS Test
      RETURN_VALUE  S_HAD_ERROR
      )

    if (S_HAD_ERROR)
      message ("Cannot submit Trilinos test results!")
    endif (S_HAD_ERROR)
  endif (CTEST_DO_SUBMIT)


endmacro(do_trilinos CONFIGURE_OPTIONS BTYPE ILOC)
