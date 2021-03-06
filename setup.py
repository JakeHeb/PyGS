'''
Created on Jun 7, 2012

@author: Steven
'''

from numpy.distutils.core import setup,Extension


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

fort = Extension('PyGS.fort.DFT',['PyGS/fort/DFT.f90'],extra_f90_compile_args=['-Wtabs'],f2py_options=['--quiet'])


if __name__=="__main__":
    
    setup(name = 'PyGS',
          version = version,
          author = 'Steven Murray',
          author_email = 'steven.jeanette.m@gmail.com',
          description = 'Interactive program to deal with Galaxy Surveys.',
          url = 'doesnt.have.one.yet.com',
          ext_modules = [fort],
          packages = ['PyGS','PyGS.fort']
    )