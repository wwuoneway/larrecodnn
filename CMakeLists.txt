# ======================================================================
#  larrecodnn main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larrecodnn/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.12...3.14 FATAL_ERROR)

project(larrecodnn LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules $ENV{ART_DIR}/Modules)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( larreco )
find_ups_product( nutools )
find_ups_product( nug4 )
find_ups_product( nurandom )
find_ups_product( art )
find_ups_product(art_root_io)
find_ups_product( cetbuildtools )
find_ups_product( postgresql )
find_ups_product( range )
find_ups_product( eigen )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# cannot build tensorflow with our clang
# at this time, trtis_clients is only available for python 3
if( DEFINED ENV{TENSORFLOW_DIR} )
  # find tensorflow library
  find_ups_product( protobuf )
  find_ups_product( tensorflow )
  cet_find_library(TENSORFLOW NAMES tensorflow-core PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
  cet_find_library(PROTOBUF NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
  # This is a hack for a broken tensorflow link, which
  # requires -lrt for glibc earlier than 2.17 (i.e. SLF6)
  # We take advantage of OSTYPE defined by cetbuildtools
  if ( "${OSTYPE}" MATCHES "slf6" )
    set(FOR_TENSORFLOW "-lrt")
  endif ()
endif ()
if( DEFINED ENV{TRTIS_CLIENTS_DIR} )
  find_ups_product(trtis_clients)
endif ()

# source
add_subdirectory(larrecodnn)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
