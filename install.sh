#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
BUILD_DIR="$DIR/build"
VERBOSE=false
project=FLINT

function usage() {
    cat <<EOF

Install script for $project

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -v, --verbose             Enable verbose output

Commands:
  build                     Perform a full build
    --compilers=<name>      Set compilers suite (intel,gnu)
    --include-orion=<path>  Build a local ORION checkout instead of a release
    --include-oslo=<path>   Build a local OSLO checkout instead of a release
    --update-deps           Select again the newest releases of ORION and OSLO (deps.lock)
    --use-cantera           Use Cantera (Sundials required)
    --use-sundials          Use Sundials (via OSLO)
    --use-tecio             Use TecIO (via ORION)

EOF
    exit 1
}


log() {
    if [ "$VERBOSE" = true ]; then
        # Bold and dim gray (ANSI escape: bold + color 90)
        echo -e "\033[1;90m$1\033[0m"
    fi
}

error() {
    # Bold red + [ERROR] tag, output to stderr
    echo -e "\033[1;31m[ERROR] $1\033[0m" >&2
}

task() {
    # Bold yellow + ==> tag, output to stdout
    echo -e "\033[1;38;5;186m==> $1\033[0m"
}


# Create default CMakePresets.json if it doesn't exist
function write_presets() {
  FC=$(grep '^CMAKE_Fortran_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CC=$(grep '^CMAKE_C_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CXX=$(grep '^CMAKE_CXX_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  # the local checkouts, when given
  local sources=""
  if [[ -n $ORION_PATH ]]; then
    sources+="        \"FETCHCONTENT_SOURCE_DIR_ORION\": \"${ORION_PATH}\","$'\n'
  fi
  if [[ -n $OSLO_PATH ]]; then
    sources+="        \"FETCHCONTENT_SOURCE_DIR_OSLO\": \"${OSLO_PATH}\","$'\n'
  fi
  
  cat <<EOF > CMakePresets.json
{
  "version": 3,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 23
  },
  "configurePresets": [
    {
      "name": "default",
      "description": "Default preset",
      "binaryDir": "\${sourceDir}/build",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "${BUILD_TYPE}",
${sources}        "ORION_TECIO_DIR": "${DEPS_DIR}/tecio",
        "OSLO_SUNDIALS_DIR": "${SUNDIALS_DIR}",
        "CMAKE_Fortran_COMPILER": "${FC}",
        "CMAKE_CXX_COMPILER": "${CXX}",
        "CMAKE_C_COMPILER": "${CC}",
        "USE_CANTERA": "${USE_CANTERA}",
        "USE_SUNDIALS": "${USE_SUNDIALS}",
        "USE_TECIO": "${USE_TECIO}"
      }
    }
  ]
}
EOF
}


# Default global values
COMMAND=""
COMPILERS=""
# ORION and OSLO: the releases of deps.lock, or the newest ones that FLINT accepts (see CMakeLists.txt)
ORION_PATH=""
OSLO_PATH=""
UPDATE_DEPS=false
# TecIO and SUNDIALS, built by ORION and OSLO, here: the build folder is removed at each build
DEPS_DIR="$DIR/build-deps"
BUILD_TYPE=RELEASE
USE_SUNDIALS=false
USE_CANTERA=false
USE_TECIO=false

# Define allowed options for each command using regular arrays
CMD=("build")
CMD_OPTIONS_build=("--compilers --include-orion --include-oslo --update-deps --use-sundials --use-cantera --use-tecio")

# Parse global options
while getopts "v-:" opt; do
    case "$opt" in
        -)
            case "$OPTARG" in
                verbose) VERBOSE=true ;;
                help) usage ;;
                *) error "Unknown global option '--$OPTARG'"; usage ;;
            esac
            ;;
        v) VERBOSE=true ;;
        ?) error "Unknown global option '-$OPTARG'"; usage ;;
    esac
done
shift $((OPTIND -1))

# Ensure a command was provided
if [[ $# -eq 0 ]]; then
  error "No command provided!"
  usage
fi

COMMAND="$1"
# Check if the command is valid
if [[ ! " ${CMD[@]} " =~ " ${COMMAND} " ]]; then
  error "Unknown command '$COMMAND'"
  usage
fi
shift

# Parse command-specific options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --compilers=*)
            [[ "$COMMAND" == "build" ]] || { error " --compilers is only valid for 'build' command"; exit 1; }
            COMPILERS="${1#*=}"
            ;;
        --include-orion=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-orion is only valid for 'build' command"; exit 1; }
            # absolute: CMake takes the folder from the build directory
            ORION_PATH=$(cd "${1#*=}" && pwd) || { error "No folder ${1#*=}"; exit 1; }
            ;;
        --include-oslo=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-oslo is only valid for 'build' command"; exit 1; }
            # absolute: CMake takes the folder from the build directory
            OSLO_PATH=$(cd "${1#*=}" && pwd) || { error "No folder ${1#*=}"; exit 1; }
            ;;
        --update-deps)
            [[ "$COMMAND" == "build" ]] || { error " --update-deps is only valid for 'build' command"; exit 1; }
            UPDATE_DEPS=true
            ;;
        --use-cantera)
            [[ "$COMMAND" == "build" ]] || { error " --use-cantera is only valid for 'build' command"; exit 1; }
            USE_CANTERA=true
            USE_SUNDIALS=true
            ;;
        --use-sundials)
            [[ "$COMMAND" == "build" ]] || { error " --use-sundials is only valid for 'build' command"; exit 1; }
            USE_SUNDIALS=true
            ;;
        --use-tecio)
            [[ "$COMMAND" == "build" ]] || { error " --use-tecio is only valid for 'build' command"; exit 1; }
            USE_TECIO=true
            ;;
        *)
            eval "opts=(\"\${CMD_OPTIONS_${COMMAND}[@]}\")"
            error "Unknown option '$1' for command '$COMMAND'. Valid options: ${opts[@]}"
            exit 1
            ;;
    esac
    shift
done


# Execute the selected command
case "$COMMAND" in
    build)
        task "Building $project"

        if [[ $COMPILERS == "intel" ]]; then 
            export FC="ifx"
            export CXX="icpx"
            export CC="icx"
        elif [[ $COMPILERS == "gnu" ]]; then 
            export FC="gfortran"
            export CXX="g++"
            export CC="gcc"
        fi
        log "Build dir: $BUILD_DIR"
        log "Build type: $BUILD_TYPE"
        log "ORION: ${ORION_PATH:-release}"
        log "OSLO: ${OSLO_PATH:-release}"
        log "Update ORION and OSLO: $UPDATE_DEPS"
        log "Use Cantera: $USE_CANTERA"
        log "Use Sundials: $USE_SUNDIALS"
        log "Use TecIO: $USE_TECIO"
        if [[ -z "${FC+x}" || -z "${CC+x}" || -z "${CXX+x}" ]]; then
          log "Compilers not set. CMake will decide."
        else
          log "Compilers: FC=$FC, CC=$CC", CXX=$CXX
        fi
        # one SUNDIALS per compilers suite (OSLO rebuilds an installation made with other compilers)
        SUNDIALS_DIR="$DEPS_DIR/sundials-install-${COMPILERS:-default}"
        CMAKE_ARGS=(-DUSE_CANTERA=$USE_CANTERA -DUSE_SUNDIALS=$USE_SUNDIALS -DCMAKE_BUILD_TYPE=$BUILD_TYPE
                    -DUSE_OPENMP=OFF -DUSE_MPI=OFF -DUSE_TECIO=$USE_TECIO
                    -DORION_TECIO_DIR="$DEPS_DIR/tecio" -DOSLO_SUNDIALS_DIR="$SUNDIALS_DIR")
        if [[ -n $ORION_PATH ]]; then
          CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_ORION="$ORION_PATH")
        fi
        if [[ -n $OSLO_PATH ]]; then
          CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_OSLO="$OSLO_PATH")
        fi
        if [[ $UPDATE_DEPS == true ]]; then
          CMAKE_ARGS+=(-DSEMVER_UPDATE=ALL)
        fi
        rm -rf $BUILD_DIR
        cmake -B $BUILD_DIR "${CMAKE_ARGS[@]}" || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"

        task "Write CMakePresets.json"
        write_presets
        log "[OK] CMakePresets.json created"
        ;;
    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac