cmake_minimum_required(VERSION 3.10)
project(troyrelu)

option(USE_BFV "Use BFV version of ReLU" ON)

if(USE_BFV)
    set(SOURCE_FILE ReLU_bfv.cu)
else()
    set(SOURCE_FILE ReLU_ckks.cu)
endif()

add_executable(troyrelu ${SOURCE_FILE})

target_link_libraries(troyrelu PRIVATE troy)
