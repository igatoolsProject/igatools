#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create the install target
macro(create_install_target)
  install(TARGETS   ${lib_name}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  install(DIRECTORY   ${PROJECT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
    PATTERN ".*" EXCLUDE
    PATTERN "*.in" EXCLUDE)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/base/config.h 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${lib_name}/base)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}Config.cmake 
    ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endmacro(create_install_target)