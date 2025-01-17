import os
import sys
from setuptools import setup, find_packages

# Define the version of the library
MAJOR_VERSION_NUM = '1'
MINOR_VERSION_NUM = '0'
BUILD_INFO = '0'


def reportError(message):
    sys.stdout.write("ERROR: ")
    sys.stdout.write(message)
    sys.stdout.write("\nExiting\n")
    sys.exit(1)


def create_init_files():
    """
    Create `__init__.py` files in required directories.
    """
    # Paths to the directories where __init__.py needs to be created
    base_wrapper_dir = os.path.join("out", "build", "wrapper")
    release_dir = os.path.join(base_wrapper_dir, "Release")

    # __init__.py content for `out/build/wrapper`
    wrapper_init_content = '''"""
Nexus: Python wrapper for the Nexus library.
"""

from .Nexus import *  # Import symbols from Nexus.py
from .Release._Nexus import *  # Import symbols from _Nexus.pyd
'''

    # __init__.py content for `out/build/wrapper/Release`
    release_init_content = '''import os
import sys

# Get the absolute path of the current directory (Release folder)
current_dir = os.path.dirname(os.path.abspath(__file__))

# Ensure the current directory is in sys.path
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Import the .pyd module
try:
    from _Nexus import *  # Import all symbols from the .pyd file
except ImportError as e:
    raise ImportError(f"Could not load _Nexus.pyd. Ensure it is in the Release folder: {current_dir}. Error: {e}")
'''

    # Create __init__.py in `out/build/wrapper`
    os.makedirs(base_wrapper_dir, exist_ok=True)
    with open(os.path.join(base_wrapper_dir, "__init__.py"), "w") as f:
        f.write(wrapper_init_content)

    # Create __init__.py in `out/build/wrapper/Release`
    os.makedirs(release_dir, exist_ok=True)
    with open(os.path.join(release_dir, "__init__.py"), "w") as f:
        f.write(release_init_content)


def buildKeywordDictionary():
    """
    Build the dictionary of keywords for the setup process.
    """
    # Find the .pyd and .py files (already built by CMake)
    wrapper_py_path = os.path.join("out", "build", "wrapper", "Nexus.py")
    wrapper_pyd_path = os.path.join("out", "build", "wrapper", "Release", "_Nexus.pyd")

    if not os.path.exists(wrapper_py_path):
        reportError(f"Missing Nexus.py at {wrapper_py_path}")
    if not os.path.exists(wrapper_pyd_path):
        reportError(f"Missing _Nexus.pyd at {wrapper_pyd_path}")

    print(f"Found Nexus.py at {wrapper_py_path}")
    print(f"Found _Nexus.pyd at {wrapper_pyd_path}")

    # Define package data and setup keywords
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
        "packages": find_packages(where="out/build/wrapper"),
        "package_dir": {"": "out/build/wrapper"},
        "package_data": {
            "": ["*.py", "*.pyd"],  # Include .py and .pyd files
        },
        "include_package_data": True,
        "data_files": [
            ("", [wrapper_py_path, wrapper_pyd_path])  # Copy these files to the package root
        ],
    }

    return setupKeywords


def main():
    if sys.version_info < (3, 6):
        reportError("Nexus requires Python 3.6 or better.")

    # Create necessary __init__.py files
    create_init_files()

    # Build the setup keywords and run setup
    setupKeywords = buildKeywordDictionary()
    setup(**setupKeywords)


if __name__ == '__main__':
    main()
