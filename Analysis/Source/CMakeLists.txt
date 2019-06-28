## loop over each source file in this directory 
foreach(source
    production
    )

   # The following folder will be included
   if(BUILD_WITH_ROOT)

      ROOT_GENERATE_DICTIONARY(G__${source} ${source}.h LINKDEF ${source}LinkDef.h)

      #---Create a shared library with geneated dictionary
      add_library( ${source} SHARED ${source}.C G__${source}.cxx)
      target_link_libraries(${source} ${ROOT_LIBRARIES})

   endif()

   if (BUILD_WITH_CUDA)
      target_link_libraries(${source} ${CUDA_LIBRARIES})
  endif()

endforeach()
