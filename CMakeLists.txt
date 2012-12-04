project(Lorg)
cmake_minimum_required(VERSION 2.8)

# The version number.
set (Lorg_VERSION_MAJOR 1)
set (Lorg_VERSION_MINOR 0)


# profiling
option(PROFILING "Enable gprof gmon.out generation" off)
if (PROFILING)
        add_definitions(-pg)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif (PROFILING)

# Boost nécessaire
find_package(Boost 1.46 COMPONENTS regex program_options system filesystem thread unit_test_framework REQUIRED)
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )
# librt nécessaire sous Linux pour Boost
if (CMAKE_SYSTEM_NAME STREQUAL Linux)
    set(RTLIB rt)
endif()


# multithreading
option(MULTITHREAD "Enable multithreading where available" on)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/")

#set(TBB_INSTALL_DIR "" CACHE STRING  "" CACHE PATH "Path to the libtbb installation folder")
#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${TBBLIB_PATH})

set(TB2GRAM 0)
if (MULTITHREAD)
    find_package(TBB)
    if (TBB_FOUND)
        add_definitions(-DUSE_THREADS=1)
        add_definitions(-DUSE_TBB)
        set(TBB_LIBS ${TBB_LIBRARIES})
        set(WITH_TB2GRAM 1)
        set(TBBINCLUDES ${TBB_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "tbb libraries not found on your system. multithreaded version unavailable.\n"
        "  - use -DTBB_INSTALL_DIR=<folder containing libtbb.so> -DTBB_INCLUDE_DIR=<folder containing tbb/tbb.h>\n"
        "  - or  -DMULTITHREAD=NO")
    endif()
else()
    if (TBB_FOUND)
        message(WARNING "As the MULTITHREAD option is not set, the tb2gram executable will not be produced (tb2gram is the grammar training tool).")
    else()
        message(WARNING "Install TBB and enable the MULTITHREAD option if you want to compile tb2gram (the grammar training tool)")
    endif()
endif()
    

# tcmalloc
option(WANT_TCMALLOC "Require use of tcmalloc" off)
if(WANT_TCMALLOC)
    find_package(GooglePerfTools COMPONENTS tcmalloc REQUIRED)
    if (NOT GOOGLE_PERFTOOLS_FOUND)
        message(FATAL_ERROR "tcmalloc not found on your system")
    endif()
    set(TCMALLOC_LIB ${TCMALLOC_LIBRARY})
endif()



# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LORG_License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${Lorg_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${Lorg_VERSION_MINOR}")
# four \ needed to escape special meaning of .
set(CPACK_SOURCE_IGNORE_FILES
     "~$"
     ".*/\\\\."
     "^${PROJECT_SOURCE_DIR}/install"
     "^${PROJECT_SOURCE_DIR}/test"
     "^${PROJECT_SOURCE_DIR}/build")
include (CPack)





# compilation options

if (CMAKE_COMPILER_IS_GNUCXX)
    set(C11_CXX_FLAGS "-std=c++0x")
    set(C11_CXX_FLAGS_RELEASE "-O6 -Wall -Wextra -pipe -DNDEBUG")
    #set(CMAKE_CXX_FLAGS_RELEASE -Ofast -Wall -Wextra -pipe
    #set(CMAKE_CXX_FLAGS_RELEASE -O3 -Wall -pipe -fomit-frame-pointer -march=core2  -mtune=core2")
    #set(CMAKE_CXX_FLAGS_RELEASE -O3 -Wall -pipe -DNDEBUG -fomit-frame-pointer")
    set(C11_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
else()
    set(C11_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(C11_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    set(C11_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()
    
set(CMAKE_CXX_FLAGS "${C11_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${C11_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "${C11_CXX_FLAGS_RELEASE}")

#message(WARNING "Compiler is GNUCXX = " ${CMAKE_COMPILER_IS_GNUCXX})
#message(WARNING "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
#message(WARNING "CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})


# targets are in src and src/training subdirectories
add_subdirectory(src)
add_subdirectory(test)
    
    