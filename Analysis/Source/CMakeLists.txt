file(GLOB_RECURSE MAIN_FILES *.C)
file(GLOB_RECURSE CUDA_FILES *.cu)

set(WHOLE_ARCHIVE_FLAG "-Wl,--whole-archive")
set(NO_WHOLE_ARCHIVE_FLAG "-Wl,--no-whole-archive")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

####################################
## handle the main files first ###
set(dict_list)
foreach(main_file ${MAIN_FILES}) 
  get_filename_component(barename ${main_file} NAME_WE)
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

