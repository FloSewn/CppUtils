set( TESTS TestCppUtils )

add_executable( ${TESTS}
  tests_Geometry.cpp
  tests_Vec2.cpp
  tests.cpp
  main.cpp
)

target_link_libraries( ${TESTS}
  cpputils
)

install( TARGETS ${TESTS} RUNTIME DESTINATION ${BIN} )
