# the geturl command needs libcurl

find_package(CURL QUIET)
option(WITH_CURL "Enable libcurl support" ${CURL_FOUND})
if(WITH_CURL)
  find_package(CURL REQUIRED)
  target_compile_definitions(lammps PRIVATE -DLAMMPS_CURL)
  target_link_libraries(lammps PRIVATE CURL::libcurl)
endif()

