
#find_package(fmt)
#find_package(spdlog)
#find_package(docopt)

# Generic test that uses conan libs
add_library(${PROJECT_NAME} INTERFACE)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_sources(${PROJECT_NAME}
  INTERFACE
    streamhist_bin.h
    streamhist_container.h
    streamhist_utils.h
    streamhist.h
  )

target_include_directories(${PROJECT_NAME}
  INTERFACE
    ${PROJECT_SOURCE_DIR}
  )


target_link_libraries(${PROJECT_NAME}
  INTERFACE 
    project_options
    project_warnings
#   docopt::docopt
#   fmt::fmt
#   spdlog::spdlog
  )
