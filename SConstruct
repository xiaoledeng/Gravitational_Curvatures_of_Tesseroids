# Modified by Xiao-Le Deng (Email: xldeng@whu.edu.cn)
# No.1 Date : 2017/01/05
# Mo.2 Date : 2017/04/17
# No.3 Date : 2018/01/24
# No.4 Date : 2020/09/29 revised for python3

# Build the Tesseroids programs
import sys
from datetime import date
import os
import fnmatch

# Make version.c by replacing stuff in version.template
from print_version import version
year = '{:d}'.format(date.today().year)
print('Replacing year ({}) and version ({}) in version.c...'.format(year,
                                                                    version))
with open('src/lib/version.template') as f:
    template = f.readlines()
with open('src/lib/version.c', 'w') as f:
    for line in template:
        f.write(line.replace('$VERSION', version).replace('$YEAR', year))


# get the mode flag from the command line
mode = ARGUMENTS.get('mode', 'default')
if not (mode in ['default', 'check', 'bin32', 'win32']):
   print("Error: unknown mode '%s'" % mode)
   Exit(1)
print('**** Compiling in ' + mode + ' mode...')

if sys.platform == 'win32':
    env = Environment(
        CPPPATH='src/lib')
elif mode == 'check':
    env = Environment(
        CFLAGS='-ansi -pedantic-errors -Wall -ggdb',
        LIBS=['m'],
        CPPPATH='src/lib')
elif mode == 'win32':
    env = Environment(
        CFLAGS='-O3',
        LIBS=['m'],
        CPPPATH='src/lib')
    env.Tool('crossmingw', toolpath=['scons-tools'])
elif mode == 'bin32':
    env = Environment(
        CFLAGS='-O3 -m32',
        LINKFLAGS='-m32',
        LIBS=['m'],
        CPPPATH='src/lib')
else:
    env = Environment(
        CFLAGS='-O3',
        LIBS=['m'],
        CPPPATH='src/lib')

# Build the tessg* programs
tesssrc = Split("""
    src/lib/logger.c
    src/lib/version.c
    src/lib/grav_tess.c
    src/lib/glq.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    src/lib/tessg_main.c
    """)
fields = ['pot', 'gx', 'gy', 'gz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz','gxxx','gxxy','gxxz','gyyx','gyyy','gyyz','gzzx','gzzy','gzzz','gxyz']
for f in fields:
    sources = ['src/tess%s.c' % (f)] + tesssrc
    env.Program('bin/tess%s' % (f), source=sources)

# Build the prismg* programs
tesssrc = Split("""
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    src/lib/prismg_main.c
    src/lib/grav_prism.c
    """)
fields = ['pot', 'gx', 'gy', 'gz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
for f in fields:
    sources = ['src/prism%s.c' % (f)] + tesssrc
    env.Program('bin/prism%s' % (f), source=sources)

# Build prismpots, prismgs, and prismggts
env.Program('bin/prismpots', source=Split("""
    src/prismpots.c
    src/lib/grav_prism_sph.c
    src/lib/grav_prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))
env.Program('bin/prismgs', source=Split("""
    src/prismgs.c
    src/lib/grav_prism_sph.c
    src/lib/grav_prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))
env.Program('bin/prismggts', source=Split("""
    src/prismggts.c
    src/lib/grav_prism_sph.c
    src/lib/grav_prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))

# Build tess2prism
env.Program('bin/tess2prism', source=Split("""
    src/tess2prism.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/geometry.c
    src/lib/parsers.c
    """))
# Build tessdefaults
env.Program('bin/tessdefaults', source=Split("""
    src/tessdefaults.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/constants.c
    src/lib/glq.c
    src/lib/geometry.c
    """))
# Build tessgrd
env.Program('bin/tessgrd', source=Split("""
    src/tessgrd.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/constants.c
    """))
# Build tessmass
env.Program('bin/tessmass', source=Split("""
    src/tessmass.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/geometry.c
    src/lib/constants.c
    """))
# Build tessmodgen
env.Program('bin/tessmodgen', source=Split("""
    src/tessmodgen.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/geometry.c
    src/lib/constants.c
    """))
# Build tesslayers
env.Program('bin/tesslayers', source=Split("""
    src/tesslayers.c
    src/lib/logger.c
    src/lib/version.c
    src/lib/parsers.c
    src/lib/geometry.c
    src/lib/constants.c
    """))

# Build the test runner
sources = ['test/test_all.c']
sources.extend(Glob("src/lib/*.c"))
tesstest = env.Program('tesstest', source=sources)

# Clean files
clean = '*.exe *.pyc *~'.split()
for root, folders, files in os.walk('.'):
    if '.git' not in root:
        for target in clean:
            for fname in fnmatch.filter(files, target):
                Clean('.', os.path.join(root, fname))
Clean('.', 'src/lib/version.c')
Clean('.', '__pycache__')
