if (APPLE)
  add_test(NAME test_xcode_build COMMAND xcodebuild -project ${CMAKE_SOURCE_DIR}/xcode/faunus/faunus.xcodeproj -alltargets)
endif()
