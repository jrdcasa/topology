import os
import glob
import site
import shutil
import numpy
import subprocess
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

# Read version from version.py
version = {}
try:
    with open("version.py") as f:
        exec(f.read(), version)
except FileNotFoundError:
    version["__version__"] = "N/A"


class BuildIndigox(build_ext):
    """
    Custom build_ext command to clone and compile indigox on Linux.
    Features:
    - Handles existing source directories
    - Checks for required C++ dependencies (Eigen3, Boost)
    - Captures and prints CMake and build output in real time
    - Searches for the built library in multiple possible locations
    - Copies the library into the Python package
    """
    def run(self):
        # --- Check for dependencies ---
        required_cmds = ["cmake", "git"]
        for cmd in required_cmds:
            if not shutil.which(cmd):
                raise RuntimeError(f"** setup.py **: Required command '{cmd}' not found. Please install it first.")

        # Check for Eigen3
        try:
            subprocess.check_output(["pkg-config", "--exists", "eigen3"])
        except subprocess.CalledProcessError:
            raise RuntimeError("** setup.py **: Eigen3 not found. Install it, e.g., 'sudo apt install libeigen3-dev'")

        # Check for Boost
        try:
            subprocess.check_output(["pkg-config", "--exists", "boost"])
        except subprocess.CalledProcessError:
            print("** setup.py **: Warning: Boost not detected via pkg-config. Ensure Boost is installed if required.")

        # --- Prepare build directories ---
        build_dir = os.path.abspath("build/indigox")
        os.makedirs(build_dir, exist_ok=True)

        repo_url = "https://github.com/allison-group/indigo-bondorder.git"
        src_dir = os.path.join(build_dir, "src")

        # Clone repo if it doesn't exist
        if not os.path.exists(src_dir):
            print("Cloning indigox repository...")
            try:
               # subprocess.check_call(["git", "clone", "--recursive", repo_url, src_dir])
                subprocess.check_call(["git", "clone", repo_url, src_dir])
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"** setup.py **: Git clone failed: {e}\n"
                    "Check your internet connection or remove the build/indigox/src folder and try again."
                )
        else:
            print(f"** setup.py **: Indigox source already exists at {src_dir}, skipping clone.")

        # Create CMake build directory
        build_subdir = os.path.join(build_dir, "build")
        os.makedirs(build_subdir, exist_ok=True)

        # Share data in indigox
        envdir = None
        for ipath in site.getsitepackages():
            g = glob.glob(os.path.join(ipath))
            if g:
                envdir = g[0]
                break

        fullpathdata_cmake = os.path.abspath(envdir+"/topology/indigox_lib")
        print("")
        print("** setup.py **: IndigoX data are located at: {}".format(fullpathdata_cmake))
        print("** setup.py **: DCMAKE_INSTALL_DATAROOTDIR = {}".format(fullpathdata_cmake))
        print("** setup.py **: Build_subdir: {}".format(build_subdir))
        print("")

        # Configure with CMake
        print("Configuring indigox with CMake...")

        cmake_arguments = [
            "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_INSTALL_DATAROOTDIR={}".format(fullpathdata_cmake),
        ]

        cmake_proc = subprocess.Popen(
            ["cmake", "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_INSTALL_PREFIX={}".format(fullpathdata_cmake),
             "-DCMAKE_INSTALL_DATAROOTDIR={}".format(fullpathdata_cmake), src_dir],
            cwd=build_subdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        for line in cmake_proc.stdout:
            print(line, end="")
        cmake_proc.wait()
        if cmake_proc.returncode != 0:
            raise RuntimeError("** setup.py **: CMake configuration failed. Check dependencies and logs above.")

        # Build with all available cores
        print("Building indigox...")
        build_proc = subprocess.Popen(
            ["cmake", "--build", ".", "--config", "Release", "-j", str(os.cpu_count())],
            cwd=build_subdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        for line in build_proc.stdout:
            print(line, end="")
        build_proc.wait()
        if build_proc.returncode != 0:
            raise RuntimeError("** setup.py **: CMake build failed. Check logs above.")

        # Build with all available cores
        print("Installing indigox...")
        install_proc = subprocess.Popen(
            ["cmake", "--install", ".", "--config", "Release"],
            cwd=build_subdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        for line in install_proc.stdout:
            print(line, end="")
        install_proc.wait()
        if install_proc.returncode != 0:
            raise RuntimeError("** setup.py **: CMake install failed. Check logs above.")

        # Find the built library
        possible_libs = [
            # os.path.join(build_subdir, "libindigox.so"),
            # os.path.join(build_subdir, "Release", "libindigox.so"),
            # os.path.join(build_subdir, "src", "libindigox.so"),
            os.path.join(build_subdir, "src/python", "pyindigox.cpython-310-x86_64-linux-gnu.so"),
        ]
        built_lib = next((p for p in possible_libs if os.path.exists(p)), None)
        if not built_lib:
            print("** setup.py **: Could not find pyindigox.cpython-310-x86_64-linux-gnu.so. Build directory contents:")
            print("** setup.py **: build_subdir:", build_subdir)
            for root, dirs, files in os.walk(build_subdir):
                for f in files:
                    print(os.path.join(root, f))
            raise FileNotFoundError("** setup.py **: pyindigox.cpython-310-x86_64-linux-gnu.so not found after build.")

        # Copy library to Python package
        target_dir = os.path.join("topology", "indigox_lib")
        os.makedirs(target_dir, exist_ok=True)
        print(f"Copying {built_lib} -> {target_dir}")
        shutil.copy2(built_lib, target_dir)

        super().run()

# Try to import Cython — fall back to .c files if not available
try:
    from Cython.Build import cythonize
    use_cython = True
except ImportError:
    def cythonize(*args, **kwargs):
        pass
    use_cython = False

# Find .pyx files inside topology/ext_libc
extlib_dir = "topology/ext_libc"
pyx_files = glob.glob(os.path.join(extlib_dir, "*.pyx"))

# If Cython is not installed, use the corresponding .c files
if use_cython:
    sources = pyx_files
else:
    sources = [f.replace(".pyx", ".c") for f in pyx_files]

# Build one Extension object per .pyx/.c module
extensions = [
    Extension(
        name=f"topology.ext_libc.{os.path.splitext(os.path.basename(src))[0]}",
        sources=[src],
        include_dirs=[numpy.get_include()],
        extra_compile_args=["-O3", "-fopenmp"],
        extra_link_args=["-fopenmp"],
    )
    for src in sources
]

# Cythonize if available
if use_cython:
    ext_modules = cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
        },
        annotate=False,
)
else:
    ext_modules = extensions

setup(
    name="topology",
    version=version["__version__"],
    author="Javier Ramos",
    author_email="jrdcasa@gmail.com",
    description="Python toolkit to create topology information of molecules",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jrdcasa/topology",
    license="GPL-3.0-or-later",
    packages=find_packages()+["topology.indigox_lib"]+["topology.ext_libc"],
    include_package_data=True,
    install_requires=[
        "numpy>=1.21.2",
        "MDAnalysis>=2.9.0",
        "networkx>=2.2",
        "periodictable>=0.4.1",
        "psutil>=5.4.3",
        "rdkit>=2.2.1",
        "GitPython",
        "Cython"
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "topology_cmd = topology_cmd.topology_cmd:main_app",
            "check_topology_pdb_psf_cmd = topology_cmd.check_topology_pdb_psf_cmd:main_app",
            "label_pdb_withoutrenumber_cmd = topology_cmd.label_pdb_withoutrenumber_cmd:main_app",
            "typing_carb_cmd = topology_cmd.typing_carb_cmd:main_app",
            "multiplePDBs_to_XTC = topology_cmd.multiplePDBs_to_XTC:main_app",
            "multiplePDBs_to_XTC_parallel = topology_cmd.multiplePDBs_to_XTC_parallel:main_app",
            "create_index_dih_carb_cmd = topology_cmd.create_index_dih_carb_cmd:main_app",
            "typing_gases_to_VLE_cmd = topology_cmd.typing_gases_to_VLE_cmd:main_app",
        ]
    },
    ext_modules=ext_modules,
    cmdclass={
        "build_ext": BuildIndigox,
    },
)
