if(NOT DEFINED CLI OR NOT EXISTS "${CLI}")
    message(FATAL_ERROR "CLI executable is missing")
endif()
if(NOT DEFINED WORK_DIR)
    message(FATAL_ERROR "WORK_DIR is required")
endif()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(WRITE "${WORK_DIR}/gray.pgm" [=[P2
4 4
255
0 32 64 96
32 64 96 128
64 96 128 160
96 128 160 192
]=])
file(WRITE "${WORK_DIR}/rgb.ppm" [=[P3
4 4
255
255 0 0 0 255 0 0 0 255 255 255 255
32 64 96 64 96 128 96 128 160 128 160 192
16 32 48 32 48 64 48 64 80 64 80 96
96 80 64 80 64 48 64 48 32 48 32 16
]=])
file(WRITE "${WORK_DIR}/initial-mask.pgm" [=[P2
4 4
255
0 255 0 255
255 0 255 0
0 255 0 255
255 0 255 0
]=])
file(WRITE "${WORK_DIR}/corrupt.pgm" "not an image\n")

function(run_cli name expected)
    execute_process(
        COMMAND "${CLI}" ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE)
    if(NOT "${result}" STREQUAL "${expected}")
        message(FATAL_ERROR "${name}: expected exit ${expected}, got ${result}; stdout='${output}' stderr='${error}'")
    endif()
    set(cli_output_${name} "${output}" PARENT_SCOPE)
    set(cli_error_${name} "${error}" PARENT_SCOPE)
endfunction()

function(require_file path)
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "missing output ${path}")
    endif()
    file(SIZE "${path}" size)
    if(size LESS 1)
        message(FATAL_ERROR "empty output ${path}")
    endif()
endfunction()

run_cli(help 0 --help)
if(NOT cli_output_help MATCHES "Usage:")
    message(FATAL_ERROR "help output has no usage line")
endif()
run_cli(version 0 --version)
if(NOT cli_output_version MATCHES "1\\.0\\.8")
    message(FATAL_ERROR "version output is incorrect: ${cli_output_version}")
endif()
run_cli(missing 2 segment)
run_cli(extra 2 segment "${WORK_DIR}/gray.pgm" --mask "${WORK_DIR}/extra.png" extra)
run_cli(bad_model 2 segment "${WORK_DIR}/gray.pgm" --mask "${WORK_DIR}/bad.png" --model invalid)
run_cli(missing_input 1 segment "${WORK_DIR}/missing.pgm" --mask "${WORK_DIR}/missing.png")
run_cli(corrupt_input 1 segment "${WORK_DIR}/corrupt.pgm" --mask "${WORK_DIR}/corrupt.png")

run_cli(gray 0 segment "${WORK_DIR}/gray.pgm" --mask "${WORK_DIR}/gray-mask.png")
require_file("${WORK_DIR}/gray-mask.png")
if(NOT cli_output_gray MATCHES "^model=chan-vese width=4 height=4 iterations=[0-9]+ stop=(converged|iteration-limit|cancelled) elapsed_ms=[0-9]+$")
    message(FATAL_ERROR "unexpected grayscale success line: ${cli_output_gray}")
endif()
run_cli(overwrite 2 segment "${WORK_DIR}/gray.pgm" --mask "${WORK_DIR}/gray-mask.png")

run_cli(initial_mask 0 segment "${WORK_DIR}/gray.pgm" --mask "${WORK_DIR}/initial-mask-output.png" --initial-mask "${WORK_DIR}/initial-mask.pgm")
require_file("${WORK_DIR}/initial-mask-output.png")
run_cli(rgb 0 segment "${WORK_DIR}/rgb.ppm" --mask "${WORK_DIR}/rgb-mask.png" --overlay "${WORK_DIR}/rgb-overlay.png")
require_file("${WORK_DIR}/rgb-mask.png")
require_file("${WORK_DIR}/rgb-overlay.png")
run_cli(geodesic 0 segment "${WORK_DIR}/rgb.ppm" --mask "${WORK_DIR}/geodesic-mask.png" --model geodesic)
require_file("${WORK_DIR}/geodesic-mask.png")
