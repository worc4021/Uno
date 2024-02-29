project(bqpd)
enable_language(CXX Fortran)

set(FORTRANSYMBOLS "")
include(FortranCInterface)
include(GNUInstallDirs)

if (HAVE_BQPD)
    add_library(bqpd bqpd/auxil.f bqpd/bqpd.f bqpd/sparseA.f bqpd/sparseL.f bqpd/util.f)
    target_include_directories(bqpd PRIVATE include)
    add_library(bqpd::bqpd ALIAS bqpd)
    target_sources(bqpd PUBLIC 
                            FILE_SET public_headers
                            TYPE HEADERS
                            BASE_DIRS include
                            FILES include/bqpd.hpp include/fc_help.h)

    set_property(TARGET bqpd PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

    list(APPEND FORTRANSYMBOLS bqpd wsc kktalphac epsc repc refactorc)
    install(TARGETS bqpd FILE_SET public_headers)
endif()

if (HAVE_MA27)
    add_library(ma27 ma27/common/deps.f ma27/ma27/ma27d.f ma27/ma28/ma28d.f ma27/mc19/mc19d.f)
    add_library(hsl::ma27 ALIAS ma27)
    target_sources(ma27  PUBLIC 
                        FILE_SET public_headers
                        TYPE HEADERS
                        BASE_DIRS include
                        FILES include/ma27.hpp include/fc_help.h)
    set_property(TARGET ma27 PROPERTY
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    list(APPEND FORTRANSYMBOLS ma27id ma27ad ma27bd ma27cd)
    install(TARGETS ma27 FILE_SET public_headers)
endif()


FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(${CMAKE_CURRENT_SOURCE_DIR}/include/fc_help.h 
                        SYMBOL_NAMESPACE FC_
                        SYMBOLS ${FORTRANSYMBOLS}
)


