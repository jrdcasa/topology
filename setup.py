import logging
import sys
import os
import glob
import site
import subprocess
import shutil
import tempfile
import urllib.request
import tarfile
from datetime import datetime
from setuptools import setup, Extension
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler

# Formatter for the logger
class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""
    FORMATS = {
        logging.ERROR: "\n\tERROR: %(asctime)s: %(msg)s",
        logging.WARNING: "\n\tWARNING: %(msg)s",
        logging.DEBUG: "%(asctime)s: %(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        date_fmt = '%d-%m-%Y %d %H:%M:%S'
        formatter = logging.Formatter(log_fmt, date_fmt)
        return formatter.format(record)


# Install packages from pip ==============================================================
def install_with_pip(pack, vers=None, log=None, namepkg=None):

    # Update pip
    p = subprocess.Popen([sys.executable, "-m", "pip", "install", "--upgrade", "pip"],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()

    # sys.executable gives the path of the python interpreter
    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if vers is None:
        m = "{}: ** {}: Installing {}".format(now, namepkg, pack)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}".format(pack)])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}".format(pack)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()
    else:
        m = "{}: ** {}: Installing {}=={}".format(now, namepkg, pack, vers)
        print(m) if log is None else log.info(m)
        # subprocess.call([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers), " &>install.log"])
        p = subprocess.Popen([sys.executable, "-m", "pip", "install", "{0}=={1}".format(pack, vers)],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.communicate()


# Install indigox-bond software ======================================================================================
def install_indigox_bond(log=None, namepkg=None):

    """
    Installing the indigo-bond library if is not present in the python enviorement.
    """

    import git

    giturl = 'https://github.com/allison-group/indigo-bondorder.git'
    install_dir = 'thirdparty/indigo-bondorder'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    m = "\n\t\t ============== COMPILING & INSTALLING INDIGO-BONDORDER PACKAGE ==============\n\n"

    try:
        import indigox as ix
        m += "{}: ** {}: indigo-bondorder is already installed in your system. {}".format(now, namepkg, giturl)
        print(m) if log is None else log.info(m)
    except ModuleNotFoundError:
        m += "{}: ** {}: indigo-bondorder is not installed in your system\n".format(now, namepkg)
        m += "{}: ** {}: Installing from git... {}\n".format(now, namepkg, giturl)
        print(m) if log is None else log.info(m)

        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            m = "================= ERROR INSTALL ================\n"
            m += "** {}: Cannot find CMake executable needed \n".format(namepkg)
            m += "** {}: for indigo-bondorder compilation.\n".format(namepkg)
            m += "** {}: Instal CMake in your Linux\n".format(namepkg)
            m += "** {}: The installation is aborted\n".format(namepkg)
            m += "================= ERROR INSTALL ================"
            print(m) if log is None else log.info(m)
            exit()

        # Look at thirdparty directory
        if os.path.isdir("thirdparty"):
            pass
        else:
            os.makedirs("thirdparty")

        # Share data in indigox
        envdir = None
        for ipath in site.getsitepackages():
            g = glob.glob(os.path.join(ipath))
            if g:
                envdir = g[0]
                break

        fullpathlib_cmake = os.path.abspath(install_dir)
        fullpathdata_cmake = os.path.abspath(envdir+"/indigox/share")
        # Check if exists a distribution of indigox in the thirdparty directory
        # git clone https://github.com/allison-group/indigo-bondorder.git
        if os.path.isdir("thirdparty/indigo-bondorder"):
            pass
        else:
            try:
                git.Repo.clone_from(giturl, install_dir)
            except git.GitCommandError:
                if not os.path.isdir(install_dir):
                    m = "================= ERROR INSTALL ================"
                    m += "** {}: The github repository for indigo-bondorder is not valid or not exists.!!!\n".format(namepkg)
                    m += "** {}: giturl     : {}\n".format(namepkg, giturl)
                    m += "** {}: install_dir: {}\n".format(namepkg, install_dir)
                    m += "** {}: Indigo-bondorder cannot be installed\n".format(namepkg)
                    m += "** {}: The installation is aborted\n".format(namepkg)
                    m += "================= ERROR INSTALL ================"
                    print(m) if log is None else log.info(m)
                    exit()
                else:
                    pass

            subprocess.call(["rm", "-rf", "thirdparty/indigo-bondorder/build"])
            subprocess.call(["mkdir", "thirdparty/indigo-bondorder/build"])
            os.chdir("thirdparty/indigo-bondorder/build")
            cmake_arguments = ["-DCMAKE_INSTALL_PREFIX={}".format(fullpathlib_cmake),
                               "-DCMAKE_INSTALL_DATAROOTDIR={}".format(fullpathdata_cmake)]
            subprocess.check_call(["cmake", "{}".format(fullpathlib_cmake)] + cmake_arguments)
            subprocess.call("make")
            subprocess.call(["make", "install"])
            os.chdir("../../")
            subprocess.call(["tar", "cvfz", "indigo-bondorder.tar.gz", "indigo-bondorder"])
            subprocess.call(["rm", "-rf", "./indigo-bondorder"])
            os.chdir("..")

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "** {}: {}\n".format(namepkg, now)
        m += "** {}: envdir={}\n".format(namepkg, envdir)
        m += "** {}: The *.so library has been installed in: {}/indigox/" \
             "pyindigox.cpython-36m-x86_64-linux-gnu.so\n".format(namepkg, envdir)
        m += "                                        {envdir}/indigox/__init__.py\n"
        m += "** {}: Share library for indigo in {}\n".format(namepkg, envdir+"/indigox/share")
        print(m) if log is None else log.info(m)

    try:
        import indigox as ix
        m = "\n{}: ** {}: indigo-bondorder is correctly imported. {}".format(now, namepkg, giturl)
        print(m) if log is None else log.info(m)
    except ModuleNotFoundError:
        m = "================= ERROR INSTALL ================\n"
        m += "{}: ** {}: indigo-bondorder libray cannot be imported as:\n".format(now, namepkg)
        m += "{}: ** {}: \timport indigox as ix\n".format(now, namepkg)
        m += "{}: ** {}: Something wrong during compilation.\n".format(now, namepkg)
        m += "================= ERROR INSTALL ================"
        print(m) if log is None else log.info(m)
        exit()

# Disabling-output-when-compiling-with-distutil =================================================
def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/
    #            7018879/disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            with open(fname, 'w') as fout:
                if include is not None:
                    fout.write('#include {0!s}\n'.format(include))
                fout.write('int main(void) {\n')
                fout.write('    {0!s};\n'.format(funcname))
                fout.write('}\n')
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir,
                                 extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except Exception:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)


# Does this compiler support OpenMP parallelization?""" ==============================================================
def detect_openmp():
    print("TOPOLOGY: Attempting to autodetect OpenMP support... ", end="")
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler.add_library('gomp')
    include = '<omp.h>'
    extra_postargs = ['-fopenmp']
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()', include=include,
                            extra_postargs=extra_postargs)
    if hasopenmp:
        print("POLYANAGRO: Compiler supports OpenMP")
    else:
        print("POLYANAGRO: Did not detect OpenMP support.")

    return hasopenmp


# Setup external extensions ==============================================================
def setup_external_extensions(debug_cflags=False, use_openmp=True):
    has_openmp = detect_openmp()

    # parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    mathlib = ['m']
    define_macros = []
    extra_compile_args = ['-std=c99', '-ffast-math', '-O3', '-funroll-loops', '-Wno-cpp']
    if debug_cflags:
        extra_compile_args.extend(['-Wall', '-pedantic'])
        define_macros.extend([('DEBUG', '1')])

    parallel_args = ['-fopenmp'] if has_openmp and use_openmp else []
    parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    parallel_macros = [('PARALLEL', None)] if has_openmp and use_openmp else []

    extensions_install = [
        Extension("ext_libc.c_distances", ["topology/ext_libc/c_distances.pyx"],
                  libraries=mathlib,
                  define_macros=define_macros,
                  extra_compile_args=extra_compile_args, ),
        Extension("ext_libc.c_distances_openmp",
                  ["topology/ext_libc/c_distances_openmp.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=parallel_args),
        Extension("ext_libc.c_unwrap_openmp", ["topology/ext_libc/c_unwrap_openmp.pyx"],
                  libraries=mathlib + parallel_libraries,
                  define_macros=define_macros + parallel_macros,
                  extra_compile_args=parallel_args + extra_compile_args,
                  extra_link_args=parallel_args),
    ]

    return extensions_install


# Install eigen library software ======================================================================================
def install_eigen(log=None, namepkg=None):

    """
    Installing the eigen library which is needed in openbabel.
    """

    import git

    # Version eigen 3.4.0 might not work. Error obtained compiler too old.
    # giturl = 'https://gitlab.com/libeigen/eigen.git'
    giturl = 'https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz'
    tar_download_file = os.path.basename(giturl)
    install_dir = 'thirdparty/eigen-3.3.9'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    m = "\n\t\t ============== COMPILING & INSTALLING EIGEN PACKAGE ==============\n\n"
    print(m) if log is None else log.info(m)

    if not os.path.isdir(install_dir+"/eigen_library/include"):
        m = "*****************************************************\n"
        m += "** {}: eigen is not installed in your system\n".format(namepkg)
        m += "http://eigen.tuxfamily.org/index.php?title=Main_Page\n"
        m += "** {}: Installing version 3.3.9 from git... {}".format(namepkg, giturl)
        print(m) if log is None else log.info(m)

    try:
        subprocess.check_output(['cmake', '--version'])
    except OSError:
        m = "================= ERROR INSTALL ================\n"
        m += "** {}: Cannot find CMake executable\n".format(namepkg)
        m += "** {}: The installation is aborted\n".format(namepkg)
        m += "================= ERROR INSTALL ================\n"
        print(m) if log is None else log.info(m)
        sys.exit()

    # Look at thirdparty directory
    if os.path.isdir("thirdparty"):
        pass
    else:
        os.makedirs("thirdparty")

    fullpath_cmake = os.path.abspath(install_dir)

    # Check if exists a distribution of indigox in the thirdparty directory
    if os.path.isdir(install_dir+"/eigen_library/include"):
        m = "{}: ** {}: eigen_library seems to be already compiled in your system. " \
            "{}".format(now, namepkg, install_dir+"/eigen_library/include")
        print(m) if log is None else log.info(m)
    else:
        # git clone https://gitlab.com/libeigen/eigen.git
        try:
            # git.Repo.clone_from(giturl, install_dir)
            urllib.request.urlretrieve(giturl, "thirdparty/"+tar_download_file)
            tar = tarfile.open("thirdparty/"+tar_download_file)
            tar.extractall(path="./thirdparty/")
            tar.close()
        except (urllib.error.HTTPError, FileNotFoundError) as e:
            if not os.path.isdir(install_dir):
                m = "================= ERROR INSTALL ================\n"
                m += "** {}: The github repository for eigen is not valid or not exists.!!!\n".format(namepkg)
                m += "** {}: giturl     : {}\n".format(namepkg, giturl)
                m += "** {}: install_dir: {}\n".format(namepkg, install_dir)
                m += "** {}: eigen cannot be installed\n".format(namepkg)
                m += "** {}: The installation is aborted\n".format(namepkg)
                m += "================= ERROR INSTALL ================"
                print(m) if log is None else log.info(m)
                exit()
            else:
                pass

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "{}: ** {}: Compiling eigen\n".format(now, namepkg)
        print(m) if log is None else log.info(m)
        subprocess.call(["rm", "-rf", install_dir+"/build"])
        subprocess.call(["mkdir",  install_dir+"/build"])
        subprocess.call(["mkdir", install_dir+"/eigen_library"])
        os.chdir(install_dir+"/build")
        cmake_arguments1 = ["-DCMAKE_INSTALL_PREFIX={}".format(fullpath_cmake+"/eigen_library")]
        subprocess.check_call(["cmake", "{}", "{}".format(fullpath_cmake)]+cmake_arguments1)
        subprocess.call("make")
        subprocess.call(["make", "install"])
        os.chdir("../../")
        os.chdir("..")

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "{}: ** {}: The eigen library has been installed in: thirdparty/eigen/eigen_library\n".format(now, namepkg)
        print(m) if log is None else log.info(m)


# Install indigox-bond software ======================================================================================
def install_openbabel(log=None, namepkg=None):

    """
    Installing the openbabel library if is not present in the python enviorement.
    """

    import git

    giturl = 'https://github.com/openbabel/openbabel.git'

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m = "\n\t\t ============== COMPILING & INSTALLING OPENBABEL PACKAGE ==============\n\n"
    print(m) if log is None else log.info(m)

    # Look at thirdparty directory
    if not os.path.isdir("thirdparty"):
        os.makedirs("thirdparty")

    install_dir = 'thirdparty'
    eigen_dir = 'thirdparty/eigen-3.3.9/eigen_library/include/eigen3'
    fullpath_cmake = os.path.abspath(install_dir)
    fullpath_eigen = os.path.abspath(eigen_dir)

    try:
        from openbabel import openbabel as ob
        m = "{}: ** {}: openbabel is already installed in your system. " \
            "{}".format(now, namepkg, fullpath_cmake + "/openbabel")
        print(m) if log is None else log.info(m)
        return
    except (ModuleNotFoundError, ImportError) as e:
        if os.path.isdir(os.path.join(fullpath_cmake, "openbabel")):
            subprocess.call(["rm", "-rf", os.path.join(fullpath_cmake, "openbabel")])
        m = "{}: Downloading: ... openbabel-3-1-1(Wait for one minute)\n".format(now)
        print(m) if log is None else log.info(m)
        git.Repo.clone_from(giturl, fullpath_cmake + "/openbabel")

    try:
        subprocess.check_output(['cmake', '--version'])
    except OSError:
        m = "================= ERROR INSTALL ================\n"
        m += "** {}: Cannot find CMake executable\n".format(namepkg)
        m += "** {}: The installation is aborted\n".format(namepkg)
        m += "================= ERROR INSTALL ================\n"
        print(m) if log is None else log.info(m)
        exit()

    # Check if swig is installed in the system. This is needed in order to build the python bindings
    error = subprocess.call(["which", "swig"])
    if error != 0:
        m = "================= ERROR INSTALL ================\n"
        m += "** {}: Cannot find Swig executable\n".format(namepkg)
        m += "** {}: Try to install swig in your system (Ubuntu: apt-get install swig)\n".format(namepkg)
        m += "** {}: The installation is aborted\n".format(namepkg)
        m += "================= ERROR INSTALL ================\n"
        print(m) if log is None else log.info(m)
        exit()

    home_directory = os.path.expanduser('~')
    if not os.path.isdir(fullpath_cmake+"/openbabel/build"):
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "{}: ** {}: Compiling openbabel\n".format(now, namepkg)
        print(m) if log is None else log.info(m)
        subprocess.call(["rm", "-rf", "thirdparty/openbabel/build"])
        subprocess.call(["mkdir", "thirdparty/openbabel/build"])
        os.chdir("thirdparty/openbabel/build")
        cmake_arguments1 = ["-DCMAKE_INSTALL_PREFIX={}".format(os.path.join(home_directory,".local/openbabel"))]
        cmake_arguments2 = ["-DPYTHON_BINDINGS=ON"]
        cmake_arguments3 = ["-DRUN_SWIG=ON"]
        cmake_arguments4 = ["-DEIGEN3_INCLUDE_DIR={}".format(fullpath_eigen)]
        subprocess.check_call(["cmake", "{}", "{}".format(fullpath_cmake+"/openbabel")] +
                              cmake_arguments1+cmake_arguments2+cmake_arguments3+cmake_arguments4)
        subprocess.call(["make", "-j4"])
        subprocess.call(["make", "install"])
        os.chdir("../../")
        os.chdir("..")

        # Copy the library to the root site-engines of the python distribution
        exe = "_openbabel.so"
        for root, dirs, files in os.walk(home_directory+"/.local/openbabel/"):
            for name in files:
                if name == exe:
                    ll = os.path.join(root, name)
                    exe_path = os.path.split(ll)[0]

    
        dir_env_python = site.getsitepackages()[0]
        #dir_openbabel_installed = os.path.join(home_directory, os.path.split(exe_path))
        dir_openbabel_installed =  exe_path
        subprocess.call(["rm", "-rf", "{}".format(os.path.join(dir_env_python, "openbabel"))])
        subprocess.call(["cp", "-r", "{}".format(dir_openbabel_installed), "{}".format(dir_env_python)])
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "{}: ** {}: Open babel copy from {} to {}\n".format(now, namepkg,
                                                                os.path.join(dir_env_python, "openbabel"), dir_env_python)
        print(m) if log is None else log.info(m)

        m += "HOME DIR      : {}\n".format(home_directory)
        m += "DIR ENV PYTHON: {}\n".format(dir_env_python)
        m += "EXE_PATH      : {}\n".format(exe_path)
        print(m) if log is None else log.info(m)


    try:
        from openbabel import openbabel as ob
        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "\n{}: ** {}: openbabel is correctly imported. {}".format(now, namepkg, giturl)
        print(m) if log is None else log.info(m)
    except (ModuleNotFoundError, ImportError) as e:
        m = "================= ERROR INSTALL ================\n"
        m += "{}: ** {}: openbabel libray cannot be imported as:\n".format(now, namepkg)
        m += "{}: ** {}: \tfrom openbabel import openbabel as ob\n".format(now, namepkg)
        m += "{}: ** {}: Something wrong during compilation.\n".format(now, namepkg)
        m += "Error: {}\n".format(e)
        m += "================= ERROR INSTALL ================"
        print(m) if log is None else log.info(m)
        exit()


# Requeriments to be manually installed  ===========================================================================
def check_requirements_outside(namepkg=None):


    # Check for swig (http://www.swig.org) ====================================
    if os.system("which swig"):
        m = "ERROR. Please install SWIG in your system (http://www.swig.org)\n"
        m += "ERROR. Ubuntu: apt get install swig"
        print(m) if logger is None else logger.info(m)
        exit()

    # Check for cmake =========================================================
    if os.system("which cmake"):
        m = "ERROR. Please install CMAKE in your system\n"
        m += "ERROR. Ubuntu: apt get install cmake"
        print(m) if logger is None else logger.info(m)
        exit()

    # Check for pygraphviz ====================================================
    try:
        import pygraphviz as pag
    except ImportError:
        m = "ERROR. Please install PYGRAPHVIZ in your system\n"
        m += "ERROR. Ubuntu: sudo apt get python3-dev\n"
        m += "ERROR. Ubuntu: sudo apt get libgraphviz-dev\n"
        m += "ERROR. Ubuntu: pip3 install wheel\n"
        m += "ERROR. Ubuntu: pip3 install pygraphviz\n"
        print(m) if logger is None else logger.info(m)
        exit()

    # Wheel can be needed for install ==========================================
    try:
        import wheel
    except ModuleNotFoundError:
        nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        m = "Install wheel in your system\n"
        m = "================= ERROR INSTALL ================\n"
        m += "{}: ** {}: wheel libray cannot be imported as:\n".format(nowm, "TOPOLOGY")
        m += "{}: ** {}: \timport wheel\n".format(nowm, "TOPOLOGY")
        m += "{}: ** {}: Install wheel in your system:\n".format(nowm, "TOPOLOGY")
        m += "{}: ** {}: \tpython -m pip install wheel\n".format(nowm, "TOPOLOGY")
        m += "================= ERROR INSTALL ================"
        print(m) if logger is None else logger.info(m)
        exit()


# Check if libraries can be imported  ===========================================================================
def last_import_check(log=None, namepkg="TOPOLOGY"):

    m1 = "\n\t\t ******************* SUMMARY *******************\n"

    # Pip packages ============================================================
    with open('requirements.txt') as f:
        required = f.read().splitlines()
    for ipack in required:
        try:
            pkg, version = ipack.split(">=")[0:2]
            if pkg[0] == "#":
                continue
        except ValueError:
            pkg = ipack
            if pkg[0] == "#":
                continue

        try:
            if pkg == "GitPython":
                pkg = "git"
            elif pkg == "rdkit-pypi":
                pkg = "rdkit"
            elif pkg == "Sphinx":
                pkg = "sphinx"
            __import__(pkg)
        except ImportError:
            m1 += "\t\t ERROR: Package {} cannot be imported.\n".format(pkg)
            m1 += "\t\t ERROR: The installation is unsuccesfully!!!!!.\n"
            print(m1) if log is None else log.info(m1)
            exit()
    m1 += "\t\t Pip packages in requirements.txt file have been succesfully imported.\n"

    # Indigox package ==========================================================
    try:
        import indigox as ix
        m1 += "\t\t Indigox has been succesfully imported.\n"
    except ImportError:
        m1 += "\t\t ERROR: Package {} cannot be imported.\n".format("indigox")
        m1 += "\t\t ERROR: The installation is unsuccesfully!!!!!.\n"
        print(m1) if log is None else log.info(m1)
        exit()

    # openbabel package ==========================================================
    try:
        import openbabel as ob
        m1 += "\t\t Openbabel has been succesfully imported.\n"
        home_directory = os.path.expanduser('~')
        m1 += "\t\t Openbabel has been installed in {}.\n".format(os.path.join(home_directory,".local/openbabel"))
    except ImportError:
        m1 += "\t\t ERROR: Package {} cannot be imported.\n".format("openbabel")
        m1 += "\t\t ERROR: The installation is unsuccesfully!!!!!.\n"
        print(m1) if log is None else log.info(m1)
        exit()

    # topology package ==========================================================
    try:
        import topology
        import topology.readmol
        m1 += "\t\t Topology has been succesfully imported.\n"
    except ImportError:
        m1 += "\t\t ERROR: Package {} cannot be imported.\n".format("Topology")
        m1 += "\t\t ERROR: The installation is unsuccesfully!!!!!.\n"
        print(m1) if log is None else log.info(m1)
        exit()

    m1 += "\t\t ******************* SUMMARY *******************"
    print(m1) if log is None else log.info(m1)


# Main setup
if __name__ == '__main__':

    # Creating the logger to install.log file ==================================
    logger = logging.getLogger(name="INSTALL_LOG")
    logger.setLevel(logging.DEBUG)
    h1 = logging.FileHandler("install.log", 'w')
    h1.setFormatter(CustomFormatter())
    # Output also in the screen
    logger.addHandler(h1)
    f1 = logging.StreamHandler()
    f1.setFormatter(CustomFormatter())
    logger.addHandler(f1)

    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Starting installation!!!! at {}\n\n".format(nowm)
    m1 += "\t\tPython version {}.{}\n".format(sys.version_info.major, sys.version_info.minor)
    print(m1) if logger is None else logger.info(m1)

    # External requirements ===================================================
    check_requirements_outside(namepkg="TOPOLOGY")

    # Install requirements ===================================
    m1 = "\t\t ============== SYS PATH ==============\n"
    for item in sys.path:
        m1 += item + "\n"
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 += "\n\t\t ============== INSTALLING PIP PACKAGES ==============\n"
    print(m1) if logger is None else logger.info(m1)

    with open('requirements.txt') as f:
        required = f.read().splitlines()
    for ipack in required:
        try:
            pkg, version = ipack.split(">=")[0:2]
            if pkg[0] == "#":
                continue
            install_with_pip(pkg, vers=version, log=logger, namepkg="TOPOLOGY")
        except ValueError:
            pkg = ipack
            if pkg[0] == "#":
                continue
            install_with_pip(pkg, log=logger, namepkg="TOPOLOGY")

    # Install third-party software ===========================
    import git
    install_indigox_bond(log=logger, namepkg="TOPOLOGY")
    install_eigen(log=logger, namepkg="TOPOLOGY")
    install_openbabel(log=logger, namepkg="TOPOLOGY")

    # Setup TOPOLOGY ===========================================
    from Cython.Build import cythonize
    import numpy
    # Extensions
    extensions = setup_external_extensions()
    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t ============== RUNNING SETUP FROM SETUPTOOLS {} ==============\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
    setup(
        ext_modules=cythonize(extensions,
                              compiler_directives={'language_level': sys.version_info[0]}),
        include_dirs=[numpy.get_include()]
    )

    last_import_check(log=logger, namepkg="TOPOLOGY")



    nowm = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    m1 = "\n\t\t Installation Done!!!! at {}\n\n".format(nowm)
    print(m1) if logger is None else logger.info(m1)
