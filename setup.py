'''
Created on Jun 7, 2012

@author: Steven
'''

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

version_info = [0,1,0]
version = '.'.join([str(i) for i in version_info])
def generate_version_py():
    fid = open("__version.py",'w')
    try:
        fid.write("version_info = %s\n" % version_info)
        fid.write("version = %s\n" % version)
    finally:
        fid.close()
        
generate_version_py()

config = Configuration()
config.add_subpackage("PyGS")
config.add_subpackage("fort")
config.add_extension('fort.DFT',['fort/DFT.f90'],extra_f90_compile_args='-Wtabs')


if __name__=="__main__":
    
    setup(name = 'PyGS',
          version = version,
          author = 'Steven Murray',
          author_email = 'steven.jeanette.m@gmail.com',
          description = 'Interactive program to deal with Galaxy Surveys.',
          url = 'doesnt.have.one.yet.com',
          **config.todict()
    )