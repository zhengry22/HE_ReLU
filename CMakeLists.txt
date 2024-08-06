add_executable(troyrelu)

target_sources(troyrelu PRIVATE
    ReLU_bfv.cu
)

target_link_libraries(troyrelu PRIVATE troy)
