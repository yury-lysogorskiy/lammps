set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2023.11.25.fix.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

set(PACELIB_MD5 "b45de9a633f42ed65422567e3ce56f9f" CACHE STRING "MD5 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_MD5)
GetFallbackURL(PACELIB_URL PACELIB_FALLBACK)

# LOCAL_ML-PACE points to top-level dir with local lammps-user-pace repo,
# to make it easier to check local build without going through the public github releases
if(LOCAL_ML-PACE)
 set(lib-pace "${LOCAL_ML-PACE}")
else()
  # download library sources to build folder
  if(EXISTS ${CMAKE_BINARY_DIR}/libpace.tar.gz)
    file(MD5 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_MD5)
  endif()
  if(NOT "${DL_MD5}" STREQUAL "${PACELIB_MD5}")
    message(STATUS "Downloading ${PACELIB_URL}")
    file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz STATUS DL_STATUS SHOW_PROGRESS)
    file(MD5 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_MD5)
    if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_MD5}" STREQUAL "${PACELIB_MD5}"))
      message(WARNING "Download from primary URL ${PACELIB_URL} failed\nTrying fallback URL ${PACELIB_FALLBACK}")
      file(DOWNLOAD ${PACELIB_FALLBACK} ${CMAKE_BINARY_DIR}/libpace.tar.gz EXPECTED_HASH MD5=${PACELIB_MD5} SHOW_PROGRESS)
    endif()
  else()
    message(STATUS "Using already downloaded archive ${CMAKE_BINARY_DIR}/libpace.tar.gz")
  endif()


  # uncompress downloaded sources
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E remove_directory lammps-user-pace*
    COMMAND ${CMAKE_COMMAND} -E tar xzf libpace.tar.gz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  get_newest_file(${CMAKE_BINARY_DIR}/lammps-user-pace-* lib-pace)
endif()

add_subdirectory(${lib-pace} build-pace)
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

# get default python
execute_process(
        COMMAND which python
        OUTPUT_VARIABLE python_name_default
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(PACE_PYTHON_EXEC  ${python_name_default})
#message("PACE_PYTHON_EXEC=${PACE_PYTHON_EXEC}")
execute_process(
  COMMAND ${PACE_PYTHON_EXEC} -c "import os;import pkgutil;package = pkgutil.get_loader('tensorflow');print(os.path.dirname(package.get_filename()))"
  OUTPUT_VARIABLE TF_DISCOVER
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
#message("TF_DISCOVER=${TF_DISCOVER}")
string(STRIP "${TF_DISCOVER}" TF_DISCOVER)
set(TF_PATH ${TF_DISCOVER})

#message("TF_PATH=${TF_PATH}")
set(TF_LIB_FILE "${TF_PATH}/libtensorflow_cc.so.2")

if(EXISTS ${TF_LIB_FILE})
  message("-- TensorFlow library is FOUND at ${TF_LIB_FILE}")
  add_library(tensorflow SHARED IMPORTED)
  set_target_properties(tensorflow PROPERTIES
          IMPORTED_LOCATION ${TF_LIB_FILE}
          INTERFACE_INCLUDE_DIRECTORIES "${TF_PATH}/include"
  )
  ###############################
  # download cppflow
  if(NOT EXISTS ${CMAKE_BINARY_DIR}/cppflow-2.0.0)
      if(NOT EXISTS ${CMAKE_BINARY_DIR}/libcppflow.tar.gz)
          set(CPPFLOW_URL "https://github.com/serizba/cppflow/archive/refs/tags/v2.0.0.tar.gz" CACHE STRING "URL for cppflow")
          message(STATUS "Downloading ${CPPFLOW_URL}")
          file(DOWNLOAD ${CPPFLOW_URL} ${CMAKE_BINARY_DIR}/libcppflow.tar.gz STATUS DL_CPPFLOW_STATUS)
          # uncompress downloaded sources
          execute_process(
                  COMMAND ${CMAKE_COMMAND} -E remove_directory cppflow-*
                  COMMAND ${CMAKE_COMMAND} -E tar xzf libcppflow.tar.gz
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
          )
      else()
        message(STATUS "Using already downloaded archive  ${CMAKE_BINARY_DIR}/libcppflow.tar.gz")
      endif()
  else()
    message(STATUS "Using already existing CPPFLOW ${CMAKE_BINARY_DIR}/cppflow-2.0.0")
  endif()

  set(cppflow_path "${CMAKE_BINARY_DIR}/cppflow-2.0.0")
  add_library(cppflow INTERFACE)
  target_include_directories(cppflow
          INTERFACE
          ${tensorflow_INCLUDE_DIRS}
          $<BUILD_INTERFACE:${cppflow_path}/include>
  )

  target_compile_features(cppflow INTERFACE cxx_std_17)
  target_link_libraries(cppflow INTERFACE
          ${tensorflow_LIBRARIES}
  )

  set(PACE_TP ON)
  find_package(OpenMP)
else()
  message("-- TensorFlow library is NOT found at ${TF_LIB_FILE}")
endif()

if(CMAKE_PROJECT_NAME STREQUAL "lammps")
  target_link_libraries(lammps PRIVATE pace)
  if(DEFINED PACE_TP)
    add_definitions(-DPACE_TP)
    target_link_libraries(lammps PRIVATE tensorflow)
    target_link_libraries(lammps PRIVATE cppflow)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(lammps PUBLIC OpenMP::OpenMP_CXX)
    endif()
  endif()

endif()
