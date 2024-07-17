project('popcorn', 'cpp',
        meson_version: '>=0.56.0',
        default_options: [ 'cpp_std=c++17',
                           'buildtype=release',
                           'warning_level=2'
                         ]
)

mpi = dependency('mpi', language: 'cpp')
openmp = dependency('openmp', language: 'cpp')

subdir('io_lib')
subdir('configuration_lib')

common_sources = [
  'allvars.cpp',
  'arvo_variables_tools.cpp',
  'arvo_tools.cpp',
  'arvo_functions.cpp',
  'finder.cpp',
  'grid.cpp',
  'objects.cpp',
  'popcorn.cpp',
  't.cpp',
]
common_lib = static_library('common',
  common_sources,
  dependencies: openmp
)

svf = executable('svf',
  [ 'svf_run.cpp' , 'popcorn.cpp', 'comm.cpp' ],
  dependencies: [mpi, openmp],
  link_with: [ common_lib, io_lib, var_lib ]
)

popcorn = executable('popcorn',
  [ 'svf.cpp', 'popcorn_run.cpp', 'comm.cpp' ],
  dependencies: [mpi, openmp],
  link_with: [ common_lib, io_lib, var_lib ]
)

clean_duplicates = executable('clean_duplicates',
  [ 'clean_duplicates_run.cpp', 'clean_duplicates.cpp' ],
  link_with: [ common_lib, io_lib, var_lib ]
)

compute_intersecs = executable('compute_intersecs',
  [ 'compute_intersecs_run.cpp', 'compute_intersecs.cpp' ],
  link_with: [ common_lib, io_lib, var_lib ],
  dependencies: openmp
)

popcorn_app = executable('popcorn_app',
  [ 'popcorn_app.cpp', 'compute_intersecs.cpp', 'clean_duplicates.cpp', 'svf.cpp' ],
  dependencies: openmp,
  link_with: [ common_lib, io_lib, var_lib ]
)

subdir('tests')
