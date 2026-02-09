# PACE library support for ML-PACE package

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(PACELIB_URL "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/v.2025.12.4.patch1.tar.gz" CACHE STRING "URL for PACE evaluator library sources")

set(PACELIB_SHA256 "1267a4d9e6a3a5f9583af29d269a492f9fbafa89d4cdc4d0a497a53aaae734ac" CACHE STRING "SHA256 checksum of PACE evaluator library tarball")
mark_as_advanced(PACELIB_URL)
mark_as_advanced(PACELIB_SHA256)
GetFallbackURL(PACELIB_URL PACELIB_FALLBACK)

# LOCAL_ML-PACE points to top-level dir with local lammps-user-pace repo,
# to make it easier to check local build without going through the public github releases
if(LOCAL_ML-PACE)
  message(STATUS "Using LOCAL ML-PACE ${LOCAL_ML-PACE}")
  set(lib-pace "${LOCAL_ML-PACE}")
else()
  # download library sources to build folder
  if(EXISTS ${CMAKE_BINARY_DIR}/libpace.tar.gz)
    file(SHA256 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_SHA256)
  endif()
  if(NOT "${DL_SHA256}" STREQUAL "${PACELIB_SHA256}")
    message(STATUS "Downloading ${PACELIB_URL}")
    file(DOWNLOAD ${PACELIB_URL} ${CMAKE_BINARY_DIR}/libpace.tar.gz
            STATUS DL_STATUS
#            SHOW_PROGRESS
    )
    file(SHA256 ${CMAKE_BINARY_DIR}/libpace.tar.gz DL_SHA256)
    if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_SHA256}" STREQUAL "${PACELIB_SHA256}"))
      message(WARNING "Download from primary URL ${PACELIB_URL} failed\nTrying fallback URL ${PACELIB_FALLBACK}")
      file(DOWNLOAD ${PACELIB_FALLBACK} ${CMAKE_BINARY_DIR}/libpace.tar.gz EXPECTED_HASH SHA256=${PACELIB_SHA256} SHOW_PROGRESS)
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

# GRACE/TensorFlow is compiled by default
if(NOT DEFINED NO_GRACE_TF)
  # We will compile with TF support

  # Check, if TF_LIB_FILE is provided
  if(TF_LIB_FILE)
    message("User-defined TF_LIB_FILE is provided: ${TF_LIB_FILE}")
  else()
    # 1) try to find TensorFlow library  from Python installation (for older versions of TF)

    # get default python
    if(NOT PACE_PYTHON_EXEC)
      find_package(Python COMPONENTS Interpreter QUIET)
      set(PACE_PYTHON_EXEC ${Python_EXECUTABLE})
    endif()
    message("Python interpreter found: ${PACE_PYTHON_EXEC}")
    execute_process(
      COMMAND ${PACE_PYTHON_EXEC} -c "import os;import pkgutil;package = pkgutil.get_loader('tensorflow');print(os.path.dirname(package.get_filename()))"
      OUTPUT_VARIABLE TF_DISCOVER
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # message("TF_DISCOVER=${TF_DISCOVER}")
    string(STRIP "${TF_DISCOVER}" TF_DISCOVER)
    set(TF_PATH ${TF_DISCOVER})

    if(APPLE)
      set(TF_LIB_FILE "${TF_PATH}/libtensorflow_cc.2.dylib")
    elseif(WIN32)
      set(TF_LIB_FILE "${TF_PATH}/tensorflow.dll")
    else()
      set(TF_LIB_FILE "${TF_PATH}/libtensorflow_cc.so.2")
    endif()
    # setup include path
    set(TF_INCLUDE_PATH "${TF_PATH}/include")


    # 2) If not found, download it
    if(NOT EXISTS ${TF_LIB_FILE})
      # Define URLs for TensorFlow C++ library for different platforms
      set(TF_URL_WINDOWS "https://storage.googleapis.com/tensorflow/versions/2.18.1/libtensorflow-cpu-windows-x86_64.zip")
      set(TF_URL_LINUX   "https://storage.googleapis.com/tensorflow/versions/2.18.0/libtensorflow-gpu-linux-x86_64.tar.gz")
      set(TF_URL_MACOS   "https://storage.googleapis.com/tensorflow/versions/2.18.0/libtensorflow-cpu-darwin-arm64.tar.gz")

      # Define SHA256 Checksums
      set(TF_SHA256_WINDOWS "28acdcea6c6b34828cf0e95e67802b0f3577d51bc2e8915de811b7aa0b04452d")
      set(TF_SHA256_LINUX   "6ca25aae03548cf76f6f68f00bdf53ec39710f08cee23bf6419b9e6e27feca5c")
      set(TF_SHA256_MACOS   "462257d2792730dcb131fcf21bc826192ae5a2c418535f6347d051f10fc8be8a")

      set(TF_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/tensorflow-library-download")

      message(STATUS "TensorFlow library not found via Python discovery. Attempting to download.")

      if(WIN32)
        set(TF_URL ${TF_URL_WINDOWS})
        set(TF_SHA256 ${TF_SHA256_WINDOWS})
        set(TF_ARCHIVE "${CMAKE_BINARY_DIR}/libtensorflow.zip")
        set(EXTRACT_COMMAND ${CMAKE_COMMAND} -E tar xf)
      elseif(APPLE)
        set(TF_URL ${TF_URL_MACOS})
        set(TF_SHA256 ${TF_SHA256_MACOS})
        set(TF_ARCHIVE "${CMAKE_BINARY_DIR}/libtensorflow.tar.gz")
        set(EXTRACT_COMMAND ${CMAKE_COMMAND} -E tar xzf)
      else() # linux
        set(TF_URL ${TF_URL_LINUX})
        set(TF_SHA256 ${TF_SHA256_LINUX})
        set(TF_ARCHIVE "${CMAKE_BINARY_DIR}/libtensorflow.tar.gz")
        set(EXTRACT_COMMAND ${CMAKE_COMMAND} -E tar xzf)
      endif()

      # 2. Download if missing (or if just deleted above)
      if(NOT EXISTS ${TF_ARCHIVE})
        message(STATUS "Downloading TensorFlow C library from ${TF_URL}")

        # Added EXPECTED_HASH to automatically verify integrity upon download
        file(DOWNLOAD ${TF_URL} ${TF_ARCHIVE}
                SHOW_PROGRESS
                EXPECTED_HASH SHA256=${TF_SHA256}
                STATUS DL_STATUS
        )

        list(GET DL_STATUS 0 DL_CODE)
        list(GET DL_STATUS 1 DL_MSG)

        if(NOT DL_CODE EQUAL 0)
          message(FATAL_ERROR "Failed to download TensorFlow from ${TF_URL}. Error: ${DL_MSG}")
        endif()
      else()
        message(STATUS "Using already downloaded archive ${TF_ARCHIVE} (Hash verified)")
      endif()

      message(STATUS "Clean folder for archive ${TF_DOWNLOAD_DIR}...")
      message(STATUS "Extracting TensorFlow library archive ${TF_ARCHIVE} to current working dir ${CMAKE_BINARY_DIR}...")
      execute_process(
        COMMAND ${EXTRACT_COMMAND} ${TF_ARCHIVE}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
      set(TF_PATH ${CMAKE_BINARY_DIR})

      # setup library path
      if(WIN32)
        set(TF_LIB_FILE "${TF_PATH}/lib/tensorflow.dll") # Path inside downloaded archive
        string(REPLACE ".dll" ".lib" TF_IMPORTS_LIB_FILE "${TF_LIB_FILE}")
      elseif(APPLE)
        set(TF_LIB_FILE "${TF_PATH}/lib/libtensorflow.2.dylib")
      else() # linux
        set(TF_LIB_FILE "${TF_PATH}/lib/libtensorflow.so.2")
      endif()

      # setup include path
      set(TF_INCLUDE_PATH "${TF_PATH}/include")
    endif()
  endif()


  # 3) Finally, import library or fail
  if(EXISTS ${TF_LIB_FILE})
    message("-- TensorFlow library is FOUND at ${TF_LIB_FILE}")
    add_library(tensorflow SHARED IMPORTED)
    if(WIN32)
      set_target_properties(tensorflow PROPERTIES
              IMPORTED_LOCATION "${TF_LIB_FILE}"
              IMPORTED_IMPLIB "${TF_IMPORTS_LIB_FILE}"
              INTERFACE_INCLUDE_DIRECTORIES "${TF_INCLUDE_PATH}")
    else()
      # This handles both Linux and macOS correctly
      set_target_properties(tensorflow PROPERTIES
              IMPORTED_LOCATION "${TF_LIB_FILE}"
              INTERFACE_INCLUDE_DIRECTORIES "${TF_INCLUDE_PATH}")
    endif()

    ###############################

    # -------------------------------------------------------------------------
    # Logic: Use local path if provided, otherwise download and extract
    # -------------------------------------------------------------------------

    # Check if cppflow_path is set and refers to a valid directory
    if(DEFINED CPPFLOW_PATH AND EXISTS "${CPPFLOW_PATH}")
      message(STATUS "Using provided local cppflow at: ${CPPFLOW_PATH}")
    else()
      message(STATUS "Local CPPFLOW_PATH not found or not set. Proceeding with download...")
      # download cppflow
      set(CPPFLOW_VERSION "2.0.3aw")
      set(CPPFLOW_URL "https://github.com/ACEworksGmbH/cppflow/archive/refs/tags/v${CPPFLOW_VERSION}.tar.gz" CACHE STRING "URL for cppflow")
      set(CPPFLOW_SHA256 "f1144030aa6d6ed8f1a843f6e5fb5ae4b8e25383620096e2d086f8c27c0a6ef0")

      set(CPPFLOW_ARCHIVE "${CMAKE_BINARY_DIR}/libcppflow.tar.gz")

      # --- START OF YOUR DOWNLOAD LOGIC ---

      # 1. Verify existing file integrity
      if(EXISTS ${CPPFLOW_ARCHIVE})
        file(SHA256 ${CPPFLOW_ARCHIVE} CURRENT_CPPFLOW_SHA256)
        if(NOT CURRENT_CPPFLOW_SHA256 STREQUAL CPPFLOW_SHA256)
          message(WARNING "Existing cppflow archive hash mismatch.\nExpected: ${CPPFLOW_SHA256}\nActual:   ${CURRENT_CPPFLOW_SHA256}\nDeleting and re-downloading...")
          file(REMOVE ${CPPFLOW_ARCHIVE})
        endif()
      endif()

      # 2. Download with hash check
      if(NOT EXISTS ${CPPFLOW_ARCHIVE})
        message(STATUS "Downloading ${CPPFLOW_URL}")
        file(DOWNLOAD ${CPPFLOW_URL} ${CPPFLOW_ARCHIVE}
#                SHOW_PROGRESS
                EXPECTED_HASH SHA256=${CPPFLOW_SHA256}
                STATUS DL_CPPFLOW_STATUS
        )

        list(GET DL_CPPFLOW_STATUS 0 DL_CPPFLOW_CODE)
        list(GET DL_CPPFLOW_STATUS 1 DL_CPPFLOW_MSG)

        if(NOT DL_CPPFLOW_CODE EQUAL 0)
          message(FATAL_ERROR "Failed to download cppflow from ${CPPFLOW_URL}. Error: ${DL_CPPFLOW_MSG}")
        endif()
      else()
        message(STATUS "Using already downloaded cppflow archive (Hash verified)")
      endif()

      # 3. Uncompress downloaded sources
      # Note: I updated the tar command to use ${CPPFLOW_ARCHIVE} variable for consistency
      execute_process(
              COMMAND ${CMAKE_COMMAND} -E remove_directory cppflow-${CPPFLOW_VERSION}
              COMMAND ${CMAKE_COMMAND} -E tar xzf ${CPPFLOW_ARCHIVE}
              WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )

      # Set the path to the auto-downloaded location
      set(CPPFLOW_PATH "${CMAKE_BINARY_DIR}/cppflow-${CPPFLOW_VERSION}")

      # --- END OF YOUR DOWNLOAD LOGIC ---

    endif()

    message("CPPFLOW_PATH=${CPPFLOW_PATH}")

    add_library(cppflow INTERFACE)
    target_include_directories(cppflow
            INTERFACE
            ${tensorflow_INCLUDE_DIRS}
            $<BUILD_INTERFACE:${CPPFLOW_PATH}/include>
    )

    target_compile_features(cppflow INTERFACE cxx_std_17)
    target_link_libraries(cppflow INTERFACE tensorflow)

    set(PACE_TP ON)
    find_package(OpenMP)
  else()
    message("-- TensorFlow library is NEITHER found at ${TF_LIB_FILE} NOR downloaded/extracted")
  endif()
else()
  message("-- NO GRACE/TensorFlow will be compiled (because flag NO_GRACE_TF is set)")
  add_definitions(-DNO_GRACE_TF=1)
endif() # if(NOT DEFINED NO_GRACE_TF)

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

  if(WIN32)
    if(BUILD_SHARED_LIBS)
      # Building lammps AND pace/yaml-cpp as shared libs (.dll)
      # Tell lammps (pair_grace.cpp) to IMPORT symbols from the DLL.
      message(STATUS "ML-PACE: Configuring 'lammps' for shared library import (YAML_CPP_DLL)")
      target_compile_definitions(lammps PRIVATE YAML_CPP_DLL)
    else()
      # Building lammps AND pace/yaml-cpp as static libs (.lib)
      # Tell lammps (pair_grace.cpp) it's a STATIC lib.
      message(STATUS "ML-PACE: Configuring 'lammps' for static library link (YAML_CPP_STATIC_DEFINE)")
      target_compile_definitions(lammps PRIVATE YAML_CPP_STATIC_DEFINE)
    endif()
  endif()

endif()
