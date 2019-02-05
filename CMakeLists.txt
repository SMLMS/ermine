cmake_minimum_required(VERSION 3.2)

project(ermine)


enable_testing()

# Define required Boost libraries
set (REQUIRED_BOOST_LIBS
	system
	filesystem
	program_options
	random
)

set (MIN_BOOST_VERSION1.5.7)

message ("Searching for Boost...\n")

find_package(
	Boost
	COMPONENTS ${REQUIRED_BOOST_LIBS}
	${MIN_BOOST_VERDION} REQUIRED
)
MESSAGE("")


# Define required Root libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set (MIN_ROOT_VERSION 5)

message("Searching for Root...\n")
find_package(
	ROOT 
	${MIN_ROOT_VERSION} REQUIRED
)
include(${ROOT_USE_FILE})
MESSAGE("")

message("Searching for HDF5...\n")
find_package(
	HDF5
	REQUIRED
	COMPONENTS CXX HL
)
MESSAGE("")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} ${HDF5_DEFINITIONS} -Wall -fPIC -std=c++17 -O2 -fopenmp")

include_directories(
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_SOURCE_DIR}/header/"
	"${PROJECT_SOURCE_DIR}/source/"
	"${Boost_INCLUDE_DIRS}"
	"${ROOT_INCLUDE_DIRS}"
	"${HDF5_INCLUDE_DIR}"
)

set(ERMINE_MAINSOURCE
	"${PROJECT_SOURCE_DIR}/main.cpp"
)

set(EXECUTABLENAME
	ermine_cli
)

#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBS OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

add_subdirectory(
	"${PROJECT_SOURCE_DIR}/source"
)


link_directories(
	${PROJECT_SOURCE_DIR}/source
	${PROJECT_SOURCE_DIR}/header
)

add_executable(
	${EXECUTABLENAME}
	${ERMINE_MAINSOURCE}
)

target_link_libraries(
	${EXECUTABLENAME}
	smlms
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${HDF5_CXX_LIBRARIES}
)

INSTALL (TARGETS
	${EXECUTABLENAME}
	DESTINATION /opt/${EXECUTABLENAME}
)
