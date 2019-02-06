###########
# Begin CMakeLists.txt
###########

cmake_minimum_required(VERSION 3.2)
project(ermine)
enable_testing()
set(onHive false) #true or false


###########
# Define Compiler
###########

if (onHive)
	set(CMAKE_CXX_COMPILER g++-6)
else()
	set(CMAKE_CXX_COMPILER g++-7)
endif(onHive)


##########
# Define required Boost libraries
##########

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BOOST_CXX_FLAGS}")
MESSAGE("")


##########
# Define required Root libraries
##########

if(onHive)
	EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
	list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
	set (MIN_ROOT_VERSION 5)

	message("Searching for Root...\n")
	find_package(
		ROOT 
		${MIN_ROOT_VERSION} REQUIRED
	)
	include(${ROOT_USE_FILE})
	MESSAGE("")
endif(onHive)
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
MESSAGE("")

if(onHive)
	EXECUTE_PROCESS(COMMAND h5c++ -show OUTPUT_VARIABLE HDF5_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)
	string(SUBSTRING ${HDF5_COMPILER} 23 -1 HDF5_CXX_FLAGS)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HDF5_CXX_FLAGS}")
else()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HDF5_DEFINITIONS}")
endif(onHive)


##########
# Define further compiler flags
##########

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -fPIC -std=c++17 -O2 -fopenmp")


##########
# Define include directories
##########

if(onHive)
	include_directories(
		${PROJECT_SOURCE_DIR}
		${PROJECT_SOURCE_DIR}/header/
		${PROJECT_SOURCE_DIR}/source/
		${Boost_INCLUDE_DIRS}
		${HDF5_INCLUDE_DIR}
	)
else()
	include_directories(
		"${PROJECT_SOURCE_DIR}"
		"${PROJECT_SOURCE_DIR}/header/"
		"${PROJECT_SOURCE_DIR}/source/"
		"${Boost_INCLUDE_DIRS}"
		"${ROOT_INCLUDE_DIRS}"
		"${HDF5_INCLUDE_DIR}"
)
endif(onHive)


##########
# Define main.cpp
##########

if(onHive)
	set(ERMINE_MAINSOURCE
		${PROJECT_SOURCE_DIR}/main.cpp
	)
else()
	set(ERMINE_MAINSOURCE
		"${PROJECT_SOURCE_DIR}/main.cpp"
	)
endif(onHive)


###########
# Link smlms Library
##########

if(onHive)
	add_subdirectory(
		${PROJECT_SOURCE_DIR}/source
	)
else()
	add_subdirectory(
		"${PROJECT_SOURCE_DIR}/source"
	)
endif(onHive)

link_directories(
	${PROJECT_SOURCE_DIR}/source
	${PROJECT_SOURCE_DIR}/header
)


##########
# Add the executable and link to external libraries
#########

set(EXECUTABLENAME
	ermine_cli
)

add_executable(
	${EXECUTABLENAME}
	${ERMINE_MAINSOURCE}
)

if(onHive)
	EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
	EXECUTE_PROCESS(COMMAND root-config --glibs OUTPUT_VARIABLE ROOT_GLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_EXE_LINKER_FLAGS "${ROOT_LIBS} ${ROOT_GLIBS}")
	target_link_libraries(
		${EXECUTABLENAME}
		smlms
		${Boost_LIBRARIES}
		${HDF5_CXX_LIBRARIES}
)
else()
	target_link_libraries(
		${EXECUTABLENAME}
		smlms
		${Boost_LIBRARIES}
		${ROOT_LIBRARIES}
		${HDF5_CXX_LIBRARIES}
)
endif(onHive)


##########
# install the executable
##########

INSTALL (TARGETS
	${EXECUTABLENAME}
	DESTINATION /opt/${EXECUTABLENAME}
)


##########
# end of CMakeLists.txt
##########
