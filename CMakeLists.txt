cmake_minimum_required(VERSION 2.8)
project( RGB2HSI )
find_package( OpenCV REQUIRED )
add_executable( rgb2hsi rgb2hsi.cpp )
target_link_libraries( rgb2hsi ${OpenCV_LIBS} )

