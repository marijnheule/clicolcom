project('gc', ['c', 'cpp'],
	default_options : ['c_std=c11', 'cpp_std=c++1z'])
compiler = meson.get_compiler('cpp')

basic_lib = static_library('basic', ['src/intstack.cpp',
						 'src/partition.cpp',
				     'src/options.cpp',
				     'src/vertices_vec.cpp',
						 'src/interval_list.cpp',
				     'src/sparse_dynamic_graph.cpp'])
basic_dep = declare_dependency(link_with:basic_lib,
			       include_directories: include_directories('src'))

prop_lib = static_library('prop', ['src/prop.cpp',
				   'src/rewriter.cpp',
				   'src/mycielski.cpp',
				   'src/statistics.cpp'])
prop_dep = declare_dependency(link_with: prop_lib,
			      include_directories: include_directories('src'))

minicsp_inc = include_directories('minicsp')
minicsp_lib = static_library('minicsp', ['minicsp/core/solver.cpp',
					 'minicsp/core/cons.cpp',
					 'minicsp/core/utils.cpp'])
minicsp_dep = declare_dependency(link_with: minicsp_lib,
				 include_directories: minicsp_inc)

segundo_lib = static_library('segundo',
		      ['sota/Segundo/DSATUR/bornes.cpp',
		       'sota/Segundo/DSATUR/dsatur_algo.cpp',
		       'sota/Segundo/DSATUR/graphe.cpp',
		       'sota/Segundo/DSATUR/cliquer/cliquer.cpp',
		       'sota/Segundo/DSATUR/cliquer/graph.cpp',
		       'sota/Segundo/DSATUR/cliquer/reorder.cpp',
		      ])
segundo_dep = declare_dependency(link_with: segundo_lib)

domega_lib = static_library('domega',
		      ['sota/dOmega/LS/src/Buss.cpp',
		       'sota/dOmega/LS/src/Clique.cpp',
		       'sota/dOmega/LS/src/Graph.cpp',
				   'sota/dOmega/LS/src/NemhauserTrotter.cpp',
					 'sota/dOmega/LS/src/VertexCover.cpp',
		      ])
domega_dep = declare_dependency(link_with: domega_lib)

thread_dep = dependency('threads')

executable('gc', ['src/gc.cpp'],
	   dependencies: [minicsp_dep, basic_dep, prop_dep,
			  segundo_dep, domega_dep, thread_dep])

executable('segundo', ['sota/Segundo/DSATUR/main.cpp'],
	   dependencies: [segundo_dep, basic_dep])

use_gecode = get_option('use_gecode')
if use_gecode
  gecode_root = get_option('gecode')
  gecode_libs = gecode_root+'/lib'
  gecode_inc = include_directories(gecode_root+'/include')
  gecodedriver_lib = compiler.find_library('gecodedriver', dirs:gecode_root+'/lib')
  gecodesearch_lib = compiler.find_library('gecodesearch', dirs:gecode_root+'/lib')
  gecodeminimodel_lib = compiler.find_library('gecodeminimodel', dirs:gecode_root+'/lib')
  gecodeint_lib = compiler.find_library('gecodeint', dirs:gecode_root+'/lib')
  gecodeset_lib = compiler.find_library('gecodeset', dirs:gecode_root+'/lib')
  gecodefloat_lib = compiler.find_library('gecodefloat', dirs:gecode_root+'/lib')
  gecodekernel_lib = compiler.find_library('gecodekernel', dirs:gecode_root+'/lib')
  gecodesupport_lib = compiler.find_library('gecodesupport', dirs:gecode_root+'/lib')
  
  gecode_dep = declare_dependency(dependencies: [gecodedriver_lib,
						 gecodesearch_lib,
						 gecodeminimodel_lib,
						 gecodeint_lib,
						 gecodeset_lib,
						 gecodefloat_lib,
						 gecodekernel_lib,
						 gecodesupport_lib],
				  include_directories: gecode_inc)
  
  executable('gecode-gc', ['sota/graph-color.cpp'],
	     dependencies: [gecode_dep, basic_dep])
endif

executable('minicolor-converter',
	   ['sota/converter/main.c',
	    'sota/converter/converter.c'])

gzip_lib = compiler.find_library('z')
gzip_dep = declare_dependency(dependencies: gzip_lib)
executable('minicolor',
	   ['sota/MiniMerge/core/Main.cpp',
	    'sota/MiniMerge/core/Solver.cpp'],
	   include_directories: include_directories(['sota/MiniMerge/mtl']),
	   dependencies: gzip_dep)

use_cplex = get_option('use_cplex')
if use_cplex
  cplex_root = get_option('cplex')
  cplex_libdir = cplex_root+'cplex/lib/x86-64_linux/static_pic'
  concert_libdir = cplex_root+'concert/lib/x86-64_linux/static_pic'
  cplex_incdir = cplex_root+'cplex/include'
  concert_incdir = cplex_root+'concert/include'

  cplex_lib = compiler.find_library('cplex', dirs:cplex_libdir)
  ilocplex_lib = compiler.find_library('ilocplex', dirs:cplex_libdir)
  concert_lib = compiler.find_library('concert', dirs:concert_libdir)

  cplex_inc = include_directories([cplex_incdir,
				   concert_incdir])
  threads_dep = dependency('threads')
  cplex_dep = declare_dependency(dependencies: [cplex_lib,
						ilocplex_lib,
						concert_lib,
						threads_dep],
				 include_directories: cplex_inc)

  executable('cplex-color',
	     ['sota/cplex.cpp'],
	     dependencies: [cplex_dep, basic_dep],
	     cpp_args: ['-DIL_STD'])
endif


executable('dsatur',
	   ['sota/trick/dsatur.c'])

executable('color6',
	   ['sota/color6/color6.c'])
		 
