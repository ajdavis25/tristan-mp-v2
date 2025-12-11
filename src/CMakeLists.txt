set(executable ${PROJECT_NAME}.xc)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
)

add_executable(${executable} ${SOURCES} ${userfile})
set_target_properties(${executable} PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_include_directories(${executable} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/algorithms)
target_include_directories(${executable} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/fieldsolvers)
target_include_directories(${executable} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/vay)
target_include_directories(${executable} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/boris)
target_include_directories(${executable} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tools)

target_link_libraries(${executable} PRIVATE MPI::MPI_Fortran)

if (${hdf5} STREQUAL "ON")
  if (TARGET HDF5::HDF5_Fortran)
    target_link_libraries(${executable} PRIVATE HDF5::HDF5_Fortran)
  else()
    target_link_libraries(${executable} PRIVATE HDF5::HDF5)
  endif()
endif()
