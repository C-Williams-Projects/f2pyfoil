project('xfoil-wrap', 'c', 'fortran', version: '0.1.0')

py = import('python').find_installation(pure: false)

incdir_numpy = run_command(py,
    ['-c', 'import numpy; print(numpy.get_include())'],
    check: true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check: true
).stdout().strip()

f2py_source = custom_target('xfoilmodule',
    input: 'src/api.f',
    output: ['xfoilmodule.c', 'xfoil-f2pywrappers.f'],
    command: [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'xfoil', '--lower'],
)

py.extension_module('xfoil',
    sources: [
        f2py_source,
        'src/naca.f',
        'src/ntcalc.f',
        'src/spline.f',
        'src/userio.f',
        'src/xbl.f',
        'src/xblsys.f',
        'src/xfoil.f',
        'src/xgdes.f',
        'src/xgeom.f',
        'src/xmdes.f',
        'src/xoper.f',
        'src/xpanel.f',
        'src/xqdes.f',
        'src/xsolve.f',
        'src/xutils.f',
        'src/osmap.f',
        'src/cst.f',
        'src/parsec.f',
        'src/api.f',
        incdir_f2py / 'fortranobject.c',
    ],
    include_directories: include_directories('src'),
    c_args: ['-I' + incdir_numpy, '-I' + incdir_f2py],
    fortran_args: ['-static', '-O', '-fdefault-real-8', '-finit-real=zero'],
    link_args: ['-static', '-lgfortran', '-lquadmath', '-lm'],
    install: true,
)

install_data(
    ['src/osmap.dat'],
    install_dir: py.get_install_dir(),
)
