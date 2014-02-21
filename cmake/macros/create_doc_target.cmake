#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a target to generate the online documentation 
# using Doxygen
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
macro(create_doc_target)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in 
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc)
  add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/doc/html/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
  add_custom_target(doc COMMAND ${CMAKE_COMMAND} -P install_doc.cmake
    DEPENDS ${CMAKE_BINARY_DIR}/doc/html/index.html 
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMENT "The documentation is found at: ${CMAKE_INSTALL_PREFIX}/doc" 
    VERBATIM)
  file(WRITE ${CMAKE_BINARY_DIR}/install_doc.cmake
    "file( INSTALL ${CMAKE_BINARY_DIR}/doc/html
    DESTINATION ${CMAKE_INSTALL_PREFIX}/doc)\n"
    "file( REMOVE_RECURSE ${CMAKE_BINARY_DIR}/doc/html)")
endmacro(create_doc_target)