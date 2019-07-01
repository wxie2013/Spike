file(GLOB_RECURSE MAIN_FILES *.C)
file(GLOB_RECURSE CUDA_FILES *.cu)

set(WHOLE_ARCHIVE_FLAG "-Wl,--whole-archive")
set(NO_WHOLE_ARCHIVE_FLAG "-Wl,--no-whole-archive")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

####################################
## handle the main files first ###
set(dict_list)
foreach(main_file ${MAIN_FILES}) 

  get_filename_component(barename ${main_file} NAME_WE) # get the name without the extension

  ROOT_GENERATE_DICTIONARY(G__${barename} ${barename}.h LINKDEF ${barename}LinkDef.h)

  #message("+++ list of files +++++")  # print out 
  #message(${barename})

  LIST(APPEND dict_list G__${barename}.cxx)
endforeach()

add_library(production SHARED
  ${MAIN_FILES}
  ${dict_list}
)

target_link_libraries(production ${ROOT_LIBRARIES})

#install(FILES  ${PROJECT_SOURCE_DIR}/Build/Source/libproduction_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES  ${CMAKE_CURRENT_BINARY_DIR}/libproduction_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_PATH})


###############################################3
### now handle the cuda files ##
#set(CUDA_NVCC_FLAGS "-std=c++14")
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch=sm_37")
#function(add_cuda_exe_lib name files libraries is_lib)
#    set (obj_list)
#    foreach(file ${files})
#        get_filename_component(barename ${file} NAME)
#        add_custom_command(
#            OUTPUT ${barename}.o
#            DEPENDS ${file}
#            COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --blue "Building NVCC Device object ${file}"
#            COMMAND ${CUDA_NVCC_EXECUTABLE} ${CUDA_NVCC_FLAGS} -I${PROJECT_SOURCE_DIR} -Xcompiler '-fPIC' -dc "${file}" -o ${barename}.o
#            COMMENT "Building ${file}"
#            VERBATIM
#        )
#
#        LIST(APPEND obj_list ${barename}.o)
#    endforeach()
#
#    set (lib_list)
#    #LIST(APPEND lib_list "-lcudadevrt")
#    #LIST(APPEND lib_list "-lcudart")
#    foreach(library_name ${libraries})
#        LIST(APPEND lib_list "-l${library_name}")
#    endforeach()
#
#    set (flags ${CUDA_NVCC_FLAGS})
#    if (is_lib)
#        LIST(APPEND flags "-dlink")
#        set (obj_name "${CMAKE_CURRENT_BINARY_DIR}/lib${name}")
#    else()
#        set (obj_name "${CMAKE_CURRENT_BINARY_DIR}/${name}")
#    endif()
#
#    add_custom_target(${name} ALL
#        COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red "Linking CXX executable ${name}"
#        COMMAND ${CUDA_NVCC_EXECUTABLE} ${flags} -Xcompiler '-fPIC' ${lib_list} -o ${obj_name}.dlink.o ${obj_list} -lcudadevrt
#        #COMMAND g++ -shared -o ${obj_name}.so ${obj_list} ${obj_name}.dlink.o ${CUDART}
#        COMMAND ${CUDA_NVCC_EXECUTABLE} ${CUDA_NVCC_FLAGS} ${lib_list} -lib -o ${obj_name}.a ${obj_list} ${obj_name}.dlink.o -lcudadevrt
#        #COMMAND ${CUDA_NVCC_EXECUTABLE} -G -g ${CUDA_NVCC_FLAGS} -Xcompiler '-fPIC' -I${PROJECT_SOURCE_DIR} -lib -o ${obj_name}.a ${files}
#
#        DEPENDS ${obj_list}
#        COMMENT "Linking ${name}"
#    )
#endfunction()
#
#function(add_cuda_exe name files libraries)
#    add_cuda_exe_lib(${name} "${files}" "${libraries}" OFF)
#endfunction()
#
#function(add_cuda_lib name files libraries)
#    add_cuda_exe_lib(${name} "${files}" "${libraries}" ON)
#endfunction()
#
#add_cuda_lib(productionCUDA
#  "${CUDA_FILES}"
#  "")
#
#add_dependencies(production
#  productionCUDA
#)

