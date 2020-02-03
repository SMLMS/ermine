##########
# begin of CMakeLists.txt
##########


##########
# Define Sources of SMLMS
##########

include_directories(
	${PROJECT_SOURCE_DIR}
)

file(GLOB SMLMS_HEADERS
	"${PROJECT_SOURCE_DIR}/header/*.hpp"
)

file(GLOB SMLMS_SOURCES
	"${PROJECT_SOURCE_DIR}/source/*.cpp"
)


##########
# Add SMLMS Library
##########

set(LIBRARYNAME
	smlms
)

add_library(
	${LIBRARYNAME} STATIC
	${SMLMS_HEADERS}
	${SMLMS_SOURCES}
)


##########
# end of CMakeLists.txt
##########
