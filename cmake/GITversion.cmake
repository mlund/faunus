# GIT
if (VERSION_STRING)
    set(GIT_LATEST_TAG ${VERSION_STRING})
else()
    find_package(Git)
    if (GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} log -1 --format="%h\ \(%cd\)" --date short
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        if (GIT_COMMIT_HASH)
            add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
        endif ()
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --tags --max-count=1
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags ${GIT_TAG_HASH}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_LATEST_TAG
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
    endif()
endif()
if (GIT_LATEST_TAG)
    add_definitions("-DGIT_LATEST_TAG=${GIT_LATEST_TAG}")
endif ()


