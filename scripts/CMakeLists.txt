add_test(
    NAME schema_examples
    COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/test_schema_examples.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
