from distutils.core import setup, Extension

mfinder_module = Extension('_mfinder',
                           sources=['/'.join(['mfinder',f]) for f in ['clustering.c',
                                                                      'globals.c',
                                                                      'grassberger.c',
                                                                      'hash.c',
                                                                      'list.c',
                                                                      'mat.c',
                                                                      'metropolis.c',
                                                                      'motif_ids.c',
                                                                      'output.c',
                                                                      'permutation.c',
                                                                      'prob.c',
                                                                      'random.c',
                                                                      'results.c',
                                                                      'role.c',
                                                                      'stubs.c',
                                                                      'switches.c',
                                                                      'wrapper.c',
                                                                      'mfinder.i',
                                                                      ]
                                    ],
                            define_macros=[('UNIX', None),
                                           #("NDEBUG", None), # appears to be enabled by default by distutils
                                           ],
                            extra_compile_args = ["-O3"],
                           )


setup(
    name = "pymfinder",
    version = "0.23",
    description = "Python wrapper for mfinder 1.2",
    author = "Daniel B. Stouffer",
    author_email = "daniel@stoufferlab.org",
    url = 'http://github.com/stoufferlab/pymfinder',
    packages = ['pymfinder',
                'pymfinder.mfinder',
                'pymfinder.test',
                #'pymfinder.data',                
                ],
    package_dir = {'pymfinder' : '',
                   'pymfinder.mfinder' : 'mfinder',
                   'pymfinder.test' : 'test',
                   #'pymfinder.data' : 'data',
                   },
    package_data = {'pymfinder.test' : ['test.net',]
                    },
    ext_package = 'pymfinder.mfinder',
    ext_modules = [mfinder_module,]

)
