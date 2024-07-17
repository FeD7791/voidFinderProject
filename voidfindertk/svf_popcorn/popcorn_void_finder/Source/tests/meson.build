gtest_dep = dependency('gtest', main : false, required : false)

testprog = executable('test_mpi',
    [ 'test_comm.cpp', 'test_mpi_main.cpp' ],
    dependencies : [ gtest_dep, mpi, openmp ],
    link_with : common_lib
)

test('MPI tests', testprog)
