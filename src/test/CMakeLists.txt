set( TESTS TestCppUtils )

add_executable( ${TESTS}
  tests_QuadTree.cpp
  tests_Geometry.cpp
  tests_MathUtility.cpp
  tests_Vec2.cpp
  tests_StringOps.cpp
  tests.cpp
  main.cpp
)

target_link_libraries( ${TESTS}
  cpputils
)

install( TARGETS ${TESTS} RUNTIME DESTINATION ${BIN} )
