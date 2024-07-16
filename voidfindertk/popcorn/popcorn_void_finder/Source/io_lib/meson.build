io_sources = [
  'io.cpp',
  'io_tracers_ascii.cpp',
  'io_tracers_gadget.cpp',
  'io_tracers_hdf5.cpp',
  'io_tracers_stream.cpp'
]

hdf5 = dependency('hdf5', language: 'cpp')

io_lib = static_library('io',
  io_sources,
  dependencies: [hdf5, openmp]
)
