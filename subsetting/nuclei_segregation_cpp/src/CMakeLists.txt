cmake_minimum_required( VERSION 2.8 )
project( NucleiSegregation )

find_package( OpenCV REQUIRED )

add_executable( ${CMAKE_SOURCE_DIR}/NucleiSegregation nuclei_segregation.cpp params.c image.cpp thresholds.cpp fg_bg.cpp watershed.cpp )

target_link_libraries( ${CMAKE_SOURCE_DIR}/NucleiSegregation ${OpenCV_LIBS} )