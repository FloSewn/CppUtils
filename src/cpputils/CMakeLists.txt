#***********************************************************
# Module: cpputils
#***********************************************************
set( MODULE_UTIL cpputils )

add_library( ${MODULE_UTIL} INTERFACE )

target_include_directories( ${MODULE_UTIL}
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} )

target_link_libraries( ${MODULE_UTIL}
  INTERFACE m )
