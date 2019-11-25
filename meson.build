project(
	'hmn-percolation',
	'cpp',
	version: '0.9',
	meson_version: '>=0.42',
	default_options: ['cpp_std=c++14','warning_level=3']
)

mpi_dep = dependency('mpi', language: 'cpp')

inc = include_directories('include')

src = []
subdir('include')
subdir('src')

framework = get_option('framework')

executable_prefix = 'hmn-percolation'

cpp = meson.get_compiler('cpp')


if not get_option('buildtype').contains('debug')
	add_global_arguments('-DNDEBUG', language: 'cpp')
else
	# debug Blitz++
	add_global_arguments('-O3', language : 'cpp')
	# add_global_arguments('-march=native', language : 'cpp')
	add_global_arguments('-DBZ_DEBUG', language: 'cpp')
endif

deps = []

if get_option('blas') == 'openblas'
	deps += cpp.find_library('blas')
endif

if get_option('lapack') == 'netlib'
	lapackdep = dependency('lapack-netlib', required: false)
	lapacklib = cpp.find_library('lapack')
	if lapackdep.found()
		deps += lapackdep
	else
		deps += lapacklib
	endif
	deps += cpp.find_library('lapacke')
	deps += cpp.find_library('gfortran')
endif


libmkl = declare_dependency(
	link_args: [
		'-lmkl_intel_lp64',
		'-lmkl_intel_thread',
		'-lmkl_core',
		'-liomp5',
		'-lpthread',
		'-lm',
		'-ldl',
	],
)

deps += libmkl


if framework == 'openmp' or framework == 'opencl'
	add_global_arguments('-DARG_OPENMP=1', language: 'cpp')
	if cpp.get_id() == 'gcc'
		add_global_arguments('-fopenmp', language: 'cpp')
		add_global_link_arguments('-fopenmp', language: 'cpp')
	elif cpp.get_id() == 'intel'
		add_global_arguments('-qopenmp', language: 'cpp')
		add_global_link_arguments('-qopenmp', language: 'cpp')
	endif
	# make Blitz++ reference counters thread-safe
	add_global_arguments('-DBZ_THREADSAFE', language: 'cpp')
endif
if framework == 'opencl'
	add_global_arguments('-DARG_OPENCL=1', language: 'cpp')
	deps += [
		dependency('OpenCL'),
		dependency('threads'),
		dependency('clFFT'),
	]
endif

if framework == 'none'
	add_global_arguments('-DARG_NONE=1', language: 'cpp')
endif


config = configuration_data()
config.set('mt_config_file', 'mt.dat')

if cpp.has_header('lapacke.h')
	config.set('ARMA_LAPACKE_H', true)
elif cpp.has_header('lapacke/lapacke.h')
	config.set('ARMA_LAPACKE_LAPACKE_H', true)
endif
if cpp.has_header('cblas.h')
	config.set('ARMA_CBLAS_H', true)
elif cpp.has_header('openblas/cblas.h')
	config.set('ARMA_OPENBLAS_CBLAS_H', true)
endif


cpp_args = []
cpp_args += deps

hm_args = '' #'-mkl=parallel'

exe = executable(
	executable_prefix,
	sources: src,
	include_directories: inc,
	dependencies: deps,
	install_dir : './',
	install: false
)

#exe = executable(executable_prefix, [src], install : true)
