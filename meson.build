project('vc-satreduce', 'cpp', license : 'MIT')

add_global_arguments(language : 'cpp')

threads = dependency('threads')
boost_dep = dependency('boost', modules : ['program_options', 'system'])

# We require CaDiCaL to be installed.
cadical_dep = dependency('cadical')

executable('vc-bnb',
           'src/algorithms.cpp',
           'src/bnb.cpp',
           'src/graph.cpp',
           'src/input-reader.cpp',
           'src/main.cpp',
           'src/reductions.cpp',
           dependencies : [
               boost_dep,
               cadical_dep,
               threads
           ],
           include_directories: [
               include_directories('include/'),
           ],
           cpp_args: ['-std=c++17'],
           install: true
)
