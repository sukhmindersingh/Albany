if(NOT ALBANY_PARALLEL_ONLY)
  # Name the test with the directory name
  get_filename_component(testName ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  # Copy Input file from source to binary dir
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_steady.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_steady.yaml COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_be_with_source.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_be_with_source.yaml COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_fe_with_source.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_fe_with_source.yaml COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compute_errs_source_no_movie2avi.m
                 ${CMAKE_CURRENT_BINARY_DIR}/compute_errs_source_no_movie2avi.m COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_backward_euler.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_backward_euler.yaml COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_forward_euler.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_forward_euler.yaml COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_gerk.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_gerk.yaml COPYONLY)

  # Create the test with this name and standard executable
  add_test(${testName}_BackwardEuler ${Albany.exe} input_backward_euler.yaml)
  set_tests_properties(${testName}_BackwardEuler PROPERTIES LABELS "Demo;Tpetra;Forward")

  add_test(${testName}_BackwardEuler_with_Source ${Albany.exe} input_be_with_source.yaml)
  set_tests_properties(${testName}_BackwardEuler_with_Source PROPERTIES LABELS "Demo;Tpetra;Forward")

  add_test(${testName}_ForwardEuler_with_Source ${Albany.exe} input_fe_with_source.yaml)
  set_tests_properties(${testName}_ForwardEuler_with_Source PROPERTIES LABELS "Demo;Tpetra;Forward")

  add_test(${testName}_ForwardEuler ${Albany.exe} input_forward_euler.yaml)
  set_tests_properties(${testName}_ForwardEuler PROPERTIES LABELS "Demo;Tpetra;Forward")
endif()
