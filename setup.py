"""
setup.py: Used for packaging the Nexus Python wrapper.
"""
import os
import sys
from setuptools import setup, find_packages

# Defining the version of the library
MAJOR_VERSION_NUM = '1'
MINOR_VERSION_NUM = '0'
BUILD_INFO = '0'

def reportError(message):
    sys.stdout.write("ERROR: ")
    sys.stdout.write(message)
    sys.stdout.write("\nExiting\n")
    sys.exit(1)

def buildKeywordDictionary():
    # Finding the .pyd and .py files
    wrapper_py_path = os.path.join("out", "build", "wrapper", "Release","Nexus.py")
    wrapper_pyd_path = os.path.join("out", "build", "wrapper", "Release", "_Nexus.pyd")
    fftw_dll_path = os.path.join("out", "build", "wrapper", "Release", "libfftw3-3.dll")

    print(f"Wrapper Python Path: {wrapper_py_path}")
    print(f"Wrapper PYD Path: {wrapper_pyd_path}")
    print(f"Nexus FFTW DLL Path: {fftw_dll_path}")

    # Checking for the .py file
    if not os.path.exists(wrapper_py_path):
        reportError(f"Missing Nexus.py at {wrapper_py_path}")
    else:
        print(f"Found Nexus.py at {wrapper_py_path}")

    # Checking for the .pyd file
    if not os.path.exists(wrapper_pyd_path):
        reportError(f"Missing _Nexus.pyd at {wrapper_pyd_path}")
    else:
        print(f"Found _Nexus.pyd at {wrapper_pyd_path}")

    # Checking for FFTW3 DLL
    if not os.path.exists(fftw_dll_path):
        sys.stdout.write(f"WARNING: FFTW3 DLL not found at {fftw_dll_path}\n")
    else:
        print(f"Found libfftw3-3.dll at {fftw_dll_path}")

    # Defining package data and setup keywords
    setupKeywords = {
        "name": "Nexus",
        "version": f"{MAJOR_VERSION_NUM}.{MINOR_VERSION_NUM}.{BUILD_INFO}",
        "author": "Omid JM, IQB Group",
        "license": "BSD-like",
        "platforms": ["Linux", "Mac OS X", "Windows"],
        "description": "Python wrapper for the Nexus library",
        "long_description": """
        Nexus is a toolkit for molecular simulation, including Python wrappers
        for efficient execution and scripting.
        """,
        "packages": find_packages(),
        "package_dir": {"": "out/build/wrapper"},
        "package_data": {
            "": ["*.py", "*.pyd", "libfftw3-3.dll"],
        },
        "include_packag e_data": True,
        "data_files": [
            ("Nexus", [wrapper_py_path, wrapper_pyd_path,fftw_dll_path])
        ],
    }

    return setupKeywords

def main():
    if sys.version_info < (3, 6):
        reportError("Nexus requires Python 3.6 or better.")

    setupKeywords = buildKeywordDictionary()
    setup(**setupKeywords)

if __name__ == '__main__':
    main()
