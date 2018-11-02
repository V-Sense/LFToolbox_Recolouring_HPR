CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
SET(PROJ_NAME cpm)
PROJECT (${PROJ_NAME})

ADD_DEFINITIONS(-std=c++11)
ADD_DEFINITIONS(-DWITH_OPENMP)

find_package(OpenCV COMPONENTS core highgui imgproc features2d REQUIRED)
find_package(LAPACK REQUIRED)

SET(EPIC_DIR EpicFlow_v1.00_CVPR2015)

INCLUDE_DIRECTORIES(
    ./
    TinyVLib
    ${EPIC_DIR}
    ${OpenCV_INCLUDE_DIRS}
    )

ADD_LIBRARY(epic_lib 
    ${EPIC_DIR}/epic.cpp
    ${EPIC_DIR}/epic_aux.cpp
    ${EPIC_DIR}/image.c
    ${EPIC_DIR}/io.c
    ${EPIC_DIR}/variational.c
    ${EPIC_DIR}/variational_aux.c
    ${EPIC_DIR}/solver.c
)

TARGET_LINK_LIBRARIES(epic_lib
    jpeg
    png
    -O2 blas lapack
    )

ADD_EXECUTABLE(${PROJ_NAME} main.cpp dt.cpp OpticalFlow.cpp)

ADD_EXECUTABLE(match2flow match2flow.cpp)
#ADD_EXECUTABLE(sift sift.cpp)
ADD_LIBRARY(match_eval_dll SHARED match_eval_dll.cpp)

TARGET_LINK_LIBRARIES(${PROJ_NAME}
    ${OpenCV_LIBRARIES}
    )

TARGET_LINK_LIBRARIES(match2flow
    ${OpenCV_LIBRARIES}
    epic_lib
    )

TARGET_LINK_LIBRARIES(match_eval_dll
    ${OpenCV_LIBRARIES}
    )

#TARGET_LINK_LIBRARIES(sift
#    ${OpenCV_LIBRARIES}
#    )
