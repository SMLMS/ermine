###########
# Begin CMakeLists.txt
###########

cmake_minimum_required(VERSION 3.2)
project(ermine)
enable_testing()

###########
# Define Compiler
###########

set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -fPIC -std=c++11 -O2 -fopenmp")

MESSAGE("compiler flags: " ${CMAKE_CXX_FLAGS})
##########
# Define required Boost libraries
##########

set (REQUIRED_BOOST_LIBS
	system
	filesystem
	program_options
	random
)

set (MIN_BOOST_VERSION 1.72.0)

message ("Searching for Boost...\n")

find_package(
	Boost
	${MIN_BOOST_VERSION} EXACT
	REQUIRED
	COMPONENTS ${REQUIRED_BOOST_LIBS}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BOOST_CXX_FLAGS}")
MESSAGE("")


##########
# Define required Root libraries
##########
#set(DROOT_DIR $ROOTSYS/cmake)
# tell CMake where to find ROOt
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set (MIN_ROOT_VERSION 5)

message("Searching for Root...\n")
find_package(
	ROOT 
	${MIN_ROOT_VERSION}
	CONFIG
	REQUIRED
	#${MIN_ROOT_VERSION} REQUIRED
)

include(${ROOT_USE_FILE})

MESSAGE("")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


##########
# Define required HDF5 libraries
##########
message("Searching for HDF5...\n")
find_package(
	HDF5
	REQUIRED
	COMPONENTS CXX HL
)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HDF5_CXX_DEFINITIONS}")

##########
# Define further compiler flags
##########

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -fPIC -std=c++11 -O2 -fopenmp")

##########
# Define main.cpp
##########
set(ERMINE_MAIN_SOURCE
	"${PROJECT_SOURCE_DIR}/main.cpp"
)

##########
# Add the executable and link to external libraries
#########

set(EXECUTABLENAME
	ermine_cli
)

###########
# Link smlms Library
##########
add_subdirectory(
	"${PROJECT_SOURCE_DIR}/source"
	"${PROJECT_SOURCE_DIR}/header"
)

##########
# Define include directories
##########
##INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

include_directories(
	${PROJECT_SOURCE_DIR}
	"${PROJECT_SOURCE_DIR}/header/"
	"${PROJECT_SOURCE_DIR}/source/"
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIRS}
	${HDF5_CXX_INCLUDE_DIRS}
	)

link_directories(
	"${PROJECT_SOURCE_DIR}/source"
	"${PROJECT_SOURCE_DIR}/header"
)

add_executable(
	${EXECUTABLENAME}
	${ERMINE_MAIN_SOURCE}
)

target_link_libraries(
	${EXECUTABLENAME}
	smlms
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	#${HDF5_CXX_LIBRARIES}
	${HDF5_CXX_HL_LIBRARIES}
)

MESSAGE("stuff..\n")
MESSAGE("Boost_INCLUDE_DIRS: " ${Boost_INCLUDE_DIRS})
MESSAGE("ROOT_INCLUDE_DIRS: " ${ROOT_INCLUDE_DIRS})
MESSAGE("HDF5_INCLUDE_DIRS: " ${HDF5_INCLUDE_DIRS})
MESSAGE("Boost_LIBRARIES: " ${Boost_LIBRARIES})
MESSAGE("ROOT_LIBRARIES: " ${ROOT_LIBRARIES})
MESSAGE("HDF5_CXX_LIBRARIES: " ${HDF5_CXX_LIBRARIES})


##########
# install the executable
##########

#INSTALL (TARGETS
#	${EXECUTABLENAME}
#	DESTINATION /opt/${EXECUTABLENAME}
#)

##########
# pack the executable
##########
#set(ERMINE_VERSION_MAJOR 19)
#set(ERMINE_VERSION_MINOR 02)
#include (InstallRequiredSystemLibraries)
#set(CPACK_GENERATOR "RPM")
#set(CPACK_PACKAGE_NAME "${EXECUTABLENAME}")
#set(CPACK_PACKAGE_DESCRIPTION "ermine command line interface")
#set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "estimation of reaction-rates from markov-based investigation on nanoscopy experiments")
#set(CPACK_PACKAGE_VERSION "${ERMINE_VERSION_MAJOR}.${ERMINE_VERSION_MINOR}")
#set(CPACK_PACKAGE_VERSION_MAJOR "${ERMINE_VERSION_MAJOR}")
#set(CPACK_PACKAGE_VERSION_MINOR "${ERMINE_VERSION_MINOR}")
#set(CPACK_PACKAGE_CONTACT "Sebastian Malkusch")
#set(CPACK_PACKAGE_VENDOR "single molecule biophysics")
#set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
#set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
#set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
#set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")

#set(CPACK_COMPONENTS_ALL Libraries ApplicationData)
#include(CPack)
 

##########
# end of CMakeLists.txt
##########
