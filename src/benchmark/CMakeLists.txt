set( BENCHMARKS BenchmarkCppUtils )

add_executable( ${BENCHMARKS}
  benchmarks_RTreeND.cpp
  benchmarks_Container.cpp
  benchmarks_QuadTree.cpp
  benchmarks.cpp
  main.cpp
)

target_link_libraries( ${BENCHMARKS}
  cpputils
)

install( TARGETS ${BENCHMARKS} RUNTIME DESTINATION ${BIN} )
