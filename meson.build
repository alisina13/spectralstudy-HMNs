project('hmn-percolation', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3'])

ff = meson.get_compiler('cpp')
#libgomp = ff.find_library('gomp')

#dependency('openmp')

dependency('mpi', language: 'cpp')



exe = executable('percolation',['FileReader.cpp','main.cpp'],
  install : true)
  
test('basic', exe)
