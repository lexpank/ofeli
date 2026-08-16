set(_install_dir "${OFELI_BINARY_DIR}/install-consumer")
set(_consumer_build "${_install_dir}/consumer-build")
set(_input "${_install_dir}/input.pgm")
set(_mask "${_install_dir}/mask.png")

file(REMOVE_RECURSE "${_install_dir}")

execute_process(
    COMMAND "${CMAKE_COMMAND}" --install "${OFELI_BINARY_DIR}" --prefix "${_install_dir}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error
)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "install failed:\n${_output}\n${_error}")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${OFELI_SOURCE_DIR}/tests/install_consumer"
            -B "${_consumer_build}" -DCMAKE_PREFIX_PATH=${_install_dir}
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error
)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "consumer configure failed:\n${_output}\n${_error}")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${_consumer_build}" --parallel
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error
)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "consumer build failed:\n${_output}\n${_error}")
endif()

if(WIN32)
    set(_consumer_exe "${_consumer_build}/ofeli_install_consumer.exe")
    set(_cli_exe "${_install_dir}/bin/ofeli-cli.exe")
else()
    set(_consumer_exe "${_consumer_build}/ofeli_install_consumer")
    set(_cli_exe "${_install_dir}/bin/ofeli-cli")
endif()

execute_process(
    COMMAND "${_consumer_exe}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error
)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "installed consumer failed:\n${_output}\n${_error}")
endif()

file(WRITE "${_input}" "P2\n4 4\n255\n0 32 64 96 32 64 96 128 64 96 128 160 96 128 160 192\n")
execute_process(
    COMMAND "${_cli_exe}" segment "${_input}" --mask "${_mask}" --force
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _error
)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "installed CLI failed:\n${_output}\n${_error}")
endif()

if(NOT EXISTS "${_mask}")
    message(FATAL_ERROR "installed CLI did not create ${_mask}")
endif()
file(SIZE "${_mask}" _mask_size)
if(_mask_size LESS 1)
    message(FATAL_ERROR "installed CLI created an empty mask")
endif()
