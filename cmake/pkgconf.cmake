# --- generate pkg-config .pc
set(pc_libs_private)
set(pc_req_private "ompi ompi-c orte zlib")

set(pc_req_public "p4est sc")
foreach(t p8est p6est)
  if(TARGET ${t})
    string(PREPEND pc_req_public "${t} ")
  endif()
endforeach()

set(pc_filename p4est-${git_version}.pc)
configure_file(${CMAKE_CURRENT_LIST_DIR}/pkgconf.pc.in ${pc_filename} @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${pc_filename} DESTINATION lib/pkgconfig)
