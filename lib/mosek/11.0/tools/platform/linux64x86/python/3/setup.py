from setuptools import Extension
from setuptools import setup
import logging
import pathlib
import platform
import setuptools
import setuptools.command.build_ext
import setuptools.command.install
import shutil
import subprocess
import sys,os,re

class InstallationError(Exception): pass

major,minor,_,_,_ = sys.version_info
setupdir = pathlib.Path(__file__).resolve().parent

python_versions = [(3, 9), (3, 10), (3, 11), (3, 12), (3, 13)]
if (major,minor) not in python_versions: raise InstallationError("Unsupported python version")


class build_ext(setuptools.command.build_ext.build_ext):
    """
    Extend the default `build_ext` command replacing the extension
    building functionality with one that simply copies a pre-built
    extension module.
    """
    def build_extension(self,ext):
        tgtdir = pathlib.Path(self.build_lib).joinpath(*ext.name.split('.')[:-1])

        try: os.makedirs(tgtdir)
        except OSError: pass
        for s in ext.sources:
            logging.info("copying %s -> %s" % (s,tgtdir))
            shutil.copy(s,tgtdir)

class install(setuptools.command.install.install):
    """
    Extend the default install command, adding an additional operation
    that installs the dynamic MOSEK libraries.
    """
    libdir   = ['../../bin']
    instlibs = [('libtbb.so.12.12', 'libtbb.so.12'), ('libmosek64.so.11.0', 'libmosek64.so.11.0')]
    
    def findlib(self,lib):
        for p in self.libdir:
            f = pathlib.Path(p).joinpath(lib)
            if f.exists():
                return f
        raise InstallationError(f"Library not found: {lib}")
    
    def install_libs(self): 
        mskdir = pathlib.Path(self.install_lib).joinpath('mosek')
        for lib,tgtname in [ (self.findlib(lib),tgtname) for (lib,tgtname) in self.instlibs ]:
            logging.info(f"copying {lib} -> {mskdir}")
            shutil.copy(lib,mskdir.joinpath(tgtname),follow_symlinks=False)
            
    def run(self):
        super().run()
        self.execute(self.install_libs, (), msg="Installing native libraries")

os.chdir(setupdir)
setup(name =        'Mosek',
      version =     '11.0.26',
      packages =    ['mosek', 'mosek.fusion', 'mosek.fusion.impl'],
      ext_modules = [ Extension("mosek._msk",sources=[r"mosek/_msk.cpython-39-x86_64-linux-gnu.so",r"mosek/_msk.cpython-310-x86_64-linux-gnu.so",r"mosek/_msk.cpython-311-x86_64-linux-gnu.so",r"mosek/_msk.cpython-312-x86_64-linux-gnu.so",r"mosek/_msk.cpython-313-x86_64-linux-gnu.so"]),Extension("mosek.fusion.impl.fragments",sources=[r"mosek/fusion/impl/fragments.cpython-39-x86_64-linux-gnu.so",r"mosek/fusion/impl/fragments.cpython-310-x86_64-linux-gnu.so",r"mosek/fusion/impl/fragments.cpython-311-x86_64-linux-gnu.so",r"mosek/fusion/impl/fragments.cpython-312-x86_64-linux-gnu.so",r"mosek/fusion/impl/fragments.cpython-313-x86_64-linux-gnu.so"]) ],
      cmdclass =    { "build_ext" : build_ext, "install" : install })
