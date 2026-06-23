cmake_minimum_required(VERSION 3.15...3.27)
project(pyigl_one_shot)

# >=17 for return value optimization
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

include(FetchContent)

# nanobind (same pin as libigl-python-bindings)
FetchContent_Declare(
  nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG        v2.7.0
)
FetchContent_MakeAvailable(nanobind)

# Build igl::core from the parent libigl checkout (this repo, which contains
# winding_number_one_shot_2d). Added as a subdirectory so libigl is NOT the
# top-level project, which disables tutorials/tests/glfw/etc. automatically and
# leaves only the header-only igl::core target.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/.." "${CMAKE_CURRENT_BINARY_DIR}/libigl")

# --- collect bindings (mirrors the libigl-python-bindings convention) --------
# Every src/<name>.cpp that defines `void bind_<name>(nb::module_ &m)` is
# auto-declared and auto-invoked in module.cpp via the generated .in files.
file(GLOB sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

set(BINDING_SOURCES ${sources})
list(FILTER BINDING_SOURCES EXCLUDE REGEX ".*/module\\.cpp$")

set(BINDING_DECLARATIONS "")
set(BINDING_INVOCATIONS "")
foreach(source_file ${BINDING_SOURCES})
  get_filename_component(filename ${source_file} NAME_WE)
  set(BINDING_DECLARATIONS "${BINDING_DECLARATIONS}extern void bind_${filename}(nb::module_ &m);\n")
  set(BINDING_INVOCATIONS  "${BINDING_INVOCATIONS}    bind_${filename}(m);\n")
endforeach()

set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/include")
file(MAKE_DIRECTORY "${generated_dir}")
file(WRITE "${generated_dir}/BINDING_DECLARATIONS.in" "${BINDING_DECLARATIONS}")
file(WRITE "${generated_dir}/BINDING_INVOCATIONS.in" "${BINDING_INVOCATIONS}")

nanobind_add_module(pyigl_core STABLE_ABI ${sources})
target_link_libraries(pyigl_core PRIVATE igl::core)
target_include_directories(pyigl_core PRIVATE "${generated_dir}")
target_include_directories(pyigl_core PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Drop the compiled extension into an `igl` package dir next to this file so it
# can be imported as `import igl` (see README for the local-shadowing caveat).
set(output_dir "${CMAKE_CURRENT_SOURCE_DIR}/igl")
file(MAKE_DIRECTORY "${output_dir}")
file(WRITE "${output_dir}/__init__.py" "from .pyigl_core import *\n")
set_target_properties(pyigl_core PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${output_dir}"
  RUNTIME_OUTPUT_DIRECTORY "${output_dir}")
