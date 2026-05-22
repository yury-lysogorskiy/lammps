# PACE library support for ML-PACE package

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2025.12.4.p1.tar.gz" CACHE STRING "URL for PACE evaluator library sources")
set(PACELIB_SHA256 "21e9d7ad2094eef0f19958d154866fc725fc6ccfa82ec3681ef2b006545ced96" CACHE STRING "SHA256 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL PACELIB_SHA256)
GetFallbackURL(PACELIB_URL PACELIB_FALLBACK)

if(LOCAL_ML-PACE)
  message(STATUS "Using LOCAL ML-PACE ${LOCAL_ML-PACE}")
  set(lib-pace "${LOCAL_ML-PACE}")
else()
  set(_pace_archive "${CMAKE_BINARY_DIR}/libpace.tar.gz")
  set(_pace_stamp   "${CMAKE_BINARY_DIR}/.pace_extracted_sha256")

  # Download archive if missing or hash mismatch
  if(EXISTS "${_pace_archive}")
    file(SHA256 "${_pace_archive}" _pace_dl_sha256)
  endif()
  if(NOT "${_pace_dl_sha256}" STREQUAL "${PACELIB_SHA256}")
    message(STATUS "Downloading ${PACELIB_URL}")
    file(DOWNLOAD "${PACELIB_URL}" "${_pace_archive}"
      EXPECTED_HASH SHA256=${PACELIB_SHA256}
      STATUS _dl_status
    )
    list(GET _dl_status 0 _dl_code)
    if(NOT _dl_code EQUAL 0)
      message(WARNING "Download from primary URL ${PACELIB_URL} failed. Trying fallback ${PACELIB_FALLBACK}")
      file(DOWNLOAD "${PACELIB_FALLBACK}" "${_pace_archive}"
        EXPECTED_HASH SHA256=${PACELIB_SHA256}
        SHOW_PROGRESS
      )
    endif()
    set(_pace_dl_sha256 "${PACELIB_SHA256}")
  else()
    message(STATUS "Using cached ${_pace_archive}")
  endif()

  # Extract only when archive changed or directory missing
  file(GLOB _pace_dirs "${CMAKE_BINARY_DIR}/lammps-user-pace*")
  if(EXISTS "${_pace_stamp}")
    file(READ "${_pace_stamp}" _pace_stamp_sha256)
    string(STRIP "${_pace_stamp_sha256}" _pace_stamp_sha256)
  endif()
  if(NOT "${_pace_stamp_sha256}" STREQUAL "${_pace_dl_sha256}" OR NOT _pace_dirs)
    file(GLOB _pace_dirs "${CMAKE_BINARY_DIR}/lammps-user-pace*")
    foreach(_d ${_pace_dirs})
      file(REMOVE_RECURSE "${_d}")
    endforeach()
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E tar xzf "${_pace_archive}"
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
      RESULT_VARIABLE _pace_tar_result
    )
    if(_pace_tar_result EQUAL 0)
      file(WRITE "${_pace_stamp}" "${_pace_dl_sha256}")
    else()
      message(FATAL_ERROR "Failed to extract ${_pace_archive} (exit code ${_pace_tar_result})")
    endif()
  endif()

  get_newest_file(${CMAKE_BINARY_DIR}/lammps-user-pace-* lib-pace)
endif()

add_subdirectory(${lib-pace} build-pace)
set_target_properties(pace PROPERTIES CXX_EXTENSIONS ON OUTPUT_NAME lammps_pace${LAMMPS_MACHINE})

if(NOT DEFINED NO_GRACE_TF)

  # --- Locate TensorFlow ---
  if(TF_LIB_FILE)
    message(STATUS "User-defined TF_LIB_FILE: ${TF_LIB_FILE}")
  else()
    # Try Python discovery (result cached across configure runs)
    if(NOT PACE_PYTHON_EXEC)
      find_program(PACE_PYTHON_EXEC NAMES python3 python)
    endif()
    message(STATUS "Python interpreter: ${PACE_PYTHON_EXEC}")
    if(NOT _tf_python_path)
      execute_process(
        COMMAND ${PACE_PYTHON_EXEC} -c
          "import os,pkgutil; p=pkgutil.get_loader('tensorflow'); print(os.path.dirname(p.get_filename()))"
        OUTPUT_VARIABLE _tf_python_path
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
      set(_tf_python_path "${_tf_python_path}" CACHE INTERNAL "TF directory discovered via Python" FORCE)
    endif()
    set(TF_PATH "${_tf_python_path}")

    if(APPLE)
      set(TF_LIB_FILE "${TF_PATH}/libtensorflow_cc.2.dylib")
    elseif(WIN32)
      set(TF_LIB_FILE "${TF_PATH}/tensorflow.dll")
    else()
      set(TF_LIB_FILE "${TF_PATH}/libtensorflow_cc.so.2")
    endif()
    set(TF_INCLUDE_PATH "${TF_PATH}/include")

    # Download TF if not found via Python
    if(NOT EXISTS "${TF_LIB_FILE}")
      if(WIN32)
        set(_tf_url     "https://storage.googleapis.com/tensorflow/versions/2.18.1/libtensorflow-cpu-windows-x86_64.zip")
        set(_tf_sha256  "28acdcea6c6b34828cf0e95e67802b0f3577d51bc2e8915de811b7aa0b04452d")
        set(_tf_archive "${CMAKE_BINARY_DIR}/libtensorflow.zip")
        set(_tf_extract ${CMAKE_COMMAND} -E tar xf)
      elseif(APPLE)
        set(_tf_url     "https://storage.googleapis.com/tensorflow/versions/2.18.0/libtensorflow-cpu-darwin-arm64.tar.gz")
        set(_tf_sha256  "462257d2792730dcb131fcf21bc826192ae5a2c418535f6347d051f10fc8be8a")
        set(_tf_archive "${CMAKE_BINARY_DIR}/libtensorflow.tar.gz")
        set(_tf_extract ${CMAKE_COMMAND} -E tar xzf)
      else()
        set(_tf_url     "https://storage.googleapis.com/tensorflow/versions/2.18.0/libtensorflow-gpu-linux-x86_64.tar.gz")
        set(_tf_sha256  "6ca25aae03548cf76f6f68f00bdf53ec39710f08cee23bf6419b9e6e27feca5c")
        set(_tf_archive "${CMAKE_BINARY_DIR}/libtensorflow.tar.gz")
        set(_tf_extract ${CMAKE_COMMAND} -E tar xzf)
      endif()

      if(NOT EXISTS "${_tf_archive}")
        message(STATUS "Downloading TensorFlow from ${_tf_url}")
        file(DOWNLOAD "${_tf_url}" "${_tf_archive}"
          SHOW_PROGRESS
          EXPECTED_HASH SHA256=${_tf_sha256}
          STATUS _tf_dl_status
        )
        list(GET _tf_dl_status 0 _tf_dl_code)
        list(GET _tf_dl_status 1 _tf_dl_msg)
        if(NOT _tf_dl_code EQUAL 0)
          message(FATAL_ERROR "Failed to download TensorFlow from ${_tf_url}. Error: ${_tf_dl_msg}")
        endif()
      else()
        message(STATUS "Using cached ${_tf_archive}")
      endif()

      # Extract only when archive changed or lib directory missing
      set(_tf_stamp "${CMAKE_BINARY_DIR}/.tf_extracted_sha256")
      if(EXISTS "${_tf_stamp}")
        file(READ "${_tf_stamp}" _tf_stamp_sha256)
        string(STRIP "${_tf_stamp_sha256}" _tf_stamp_sha256)
      endif()
      if(NOT "${_tf_stamp_sha256}" STREQUAL "${_tf_sha256}" OR NOT EXISTS "${CMAKE_BINARY_DIR}/lib")
        message(STATUS "Extracting TensorFlow to ${CMAKE_BINARY_DIR}...")
        execute_process(
          COMMAND ${_tf_extract} "${_tf_archive}"
          WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
          RESULT_VARIABLE _tf_tar_result
        )
        if(_tf_tar_result EQUAL 0)
          file(WRITE "${_tf_stamp}" "${_tf_sha256}")
        else()
          message(FATAL_ERROR "Failed to extract ${_tf_archive} (exit code ${_tf_tar_result})")
        endif()
      endif()

      set(TF_PATH "${CMAKE_BINARY_DIR}")
      if(WIN32)
        set(TF_LIB_FILE "${TF_PATH}/lib/tensorflow.dll")
        string(REPLACE ".dll" ".lib" TF_IMPORTS_LIB_FILE "${TF_LIB_FILE}")
      elseif(APPLE)
        set(TF_LIB_FILE "${TF_PATH}/lib/libtensorflow.2.dylib")
      else()
        set(TF_LIB_FILE "${TF_PATH}/lib/libtensorflow.so.2")
      endif()
      set(TF_INCLUDE_PATH "${TF_PATH}/include")
    endif()
  endif()

  # --- Import TensorFlow target ---
  if(EXISTS "${TF_LIB_FILE}")
    message(STATUS "TensorFlow library found at ${TF_LIB_FILE}")
    add_library(tensorflow SHARED IMPORTED)
    if(WIN32)
      set_target_properties(tensorflow PROPERTIES
        IMPORTED_LOCATION    "${TF_LIB_FILE}"
        IMPORTED_IMPLIB      "${TF_IMPORTS_LIB_FILE}"
        INTERFACE_INCLUDE_DIRECTORIES "${TF_INCLUDE_PATH}")
    else()
      get_filename_component(_tf_lib_dir "${TF_LIB_FILE}" DIRECTORY)
      # find_library can't match versioned .so.2 files, so check directly
      find_library(TF_FRAMEWORK_LIB NAMES tensorflow_framework
        HINTS "${_tf_lib_dir}" NO_DEFAULT_PATH)
      if(NOT TF_FRAMEWORK_LIB)
        file(GLOB _tf_fw_candidates "${_tf_lib_dir}/libtensorflow_framework.so*"
                                    "${_tf_lib_dir}/libtensorflow_framework.*dylib")
        list(GET _tf_fw_candidates 0 TF_FRAMEWORK_LIB)
      endif()
      if(TF_FRAMEWORK_LIB)
        message(STATUS "TensorFlow framework library found at ${TF_FRAMEWORK_LIB}")
        set(_tf_iface_libs "${TF_FRAMEWORK_LIB}")
      else()
        set(_tf_iface_libs "")
      endif()
      set_target_properties(tensorflow PROPERTIES
        IMPORTED_LOCATION    "${TF_LIB_FILE}"
        INTERFACE_INCLUDE_DIRECTORIES "${TF_INCLUDE_PATH}"
        INTERFACE_LINK_LIBRARIES "${_tf_iface_libs}")
    endif()

    # In newer TF versions (e.g. 2.20+), TF_Message and other symbols were moved
    # to libtensorflow_framework.so.2, which must be explicitly linked if present.
    if(APPLE)
      set(TF_FW_LIB_FILE "${TF_PATH}/libtensorflow_framework.2.dylib")
    elseif(WIN32)
      set(TF_FW_LIB_FILE "${TF_PATH}/tensorflow_framework.dll")
    else()
      set(TF_FW_LIB_FILE "${TF_PATH}/libtensorflow_framework.so.2")
    endif()

    if(EXISTS ${TF_FW_LIB_FILE})
      message("-- TensorFlow framework library is FOUND at ${TF_FW_LIB_FILE}")
      add_library(tensorflow_framework SHARED IMPORTED)
      set_target_properties(tensorflow_framework PROPERTIES
              IMPORTED_LOCATION "${TF_FW_LIB_FILE}")
    endif()

    # --- cppflow ---
    if(DEFINED CPPFLOW_PATH AND EXISTS "${CPPFLOW_PATH}")
      message(STATUS "Using provided local cppflow at: ${CPPFLOW_PATH}")
    else()
      set(CPPFLOW_VERSION "2.0.4aw")
      set(CPPFLOW_URL    "https://github.com/ACEworksGmbH/cppflow/archive/refs/tags/v${CPPFLOW_VERSION}.tar.gz" CACHE STRING "URL for cppflow")
      set(CPPFLOW_SHA256 "10a0956dd5acb3515a94ef8172261c77494e70d3116ba078b3aac10b8be199e5")
      set(_cppflow_archive "${CMAKE_BINARY_DIR}/libcppflow.tar.gz")
      set(_cppflow_stamp   "${CMAKE_BINARY_DIR}/.cppflow_extracted_sha256")

      if(NOT EXISTS "${_cppflow_archive}")
        message(STATUS "Downloading cppflow from ${CPPFLOW_URL}")
        file(DOWNLOAD "${CPPFLOW_URL}" "${_cppflow_archive}"
          EXPECTED_HASH SHA256=${CPPFLOW_SHA256}
          STATUS _cppflow_dl_status
        )
        list(GET _cppflow_dl_status 0 _cppflow_dl_code)
        list(GET _cppflow_dl_status 1 _cppflow_dl_msg)
        if(NOT _cppflow_dl_code EQUAL 0)
          message(FATAL_ERROR "Failed to download cppflow from ${CPPFLOW_URL}. Error: ${_cppflow_dl_msg}")
        endif()
      else()
        message(STATUS "Using cached ${_cppflow_archive}")
      endif()

      if(EXISTS "${_cppflow_stamp}")
        file(READ "${_cppflow_stamp}" _cppflow_stamp_sha256)
        string(STRIP "${_cppflow_stamp_sha256}" _cppflow_stamp_sha256)
      endif()
      set(_cppflow_dir "${CMAKE_BINARY_DIR}/cppflow-${CPPFLOW_VERSION}")
      if(NOT "${_cppflow_stamp_sha256}" STREQUAL "${CPPFLOW_SHA256}" OR NOT EXISTS "${_cppflow_dir}")
        file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/cppflow-${CPPFLOW_VERSION}")
        execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xzf "${_cppflow_archive}"
          WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
          RESULT_VARIABLE _cppflow_tar_result
        )
        if(_cppflow_tar_result EQUAL 0)
          file(WRITE "${_cppflow_stamp}" "${CPPFLOW_SHA256}")
        else()
          message(FATAL_ERROR "Failed to extract ${_cppflow_archive} (exit code ${_cppflow_tar_result})")
        endif()
      endif()

      set(CPPFLOW_PATH "${CMAKE_BINARY_DIR}/cppflow-${CPPFLOW_VERSION}")
    endif()

    add_library(cppflow INTERFACE)
    target_include_directories(cppflow INTERFACE $<BUILD_INTERFACE:${CPPFLOW_PATH}/include>)
    target_compile_features(cppflow INTERFACE cxx_std_17)
    target_link_libraries(cppflow INTERFACE tensorflow)

    set(PACE_TP ON)
    find_package(OpenMP)
  else()
    message(STATUS "TensorFlow not found; GRACE/TF pair styles will not be available")
  endif()

else()
  message(STATUS "NO_GRACE_TF set: skipping TensorFlow/GRACE compilation")
endif()

if(CMAKE_PROJECT_NAME STREQUAL "lammps")
  target_link_libraries(lammps PRIVATE pace)

  if(DEFINED PACE_TP)
    target_compile_definitions(lammps PRIVATE PACE_TP)
    target_link_libraries(lammps PRIVATE tensorflow)
    if(TARGET tensorflow_framework)
      target_link_libraries(lammps PRIVATE tensorflow_framework)
    endif()
    target_link_libraries(lammps PRIVATE cppflow)
    if(OpenMP_CXX_FOUND)
      target_link_libraries(lammps PUBLIC OpenMP::OpenMP_CXX)
    endif()
  endif()

  if(DEFINED NO_GRACE_TF)
    target_compile_definitions(lammps PRIVATE NO_GRACE_TF=1)
  endif()

  if(MLPACE_DO_NOT_DISABLE_TFLOAT32)
    target_compile_definitions(lammps PRIVATE MLPACE_DO_NOT_DISABLE_TFLOAT32)
  endif()

  if(WIN32)
    if(BUILD_SHARED_LIBS)
      target_compile_definitions(lammps PRIVATE YAML_CPP_DLL)
    else()
      target_compile_definitions(lammps PRIVATE YAML_CPP_STATIC_DEFINE)
    endif()
  endif()
endif()
