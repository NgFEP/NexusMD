# Nexus

A demo version of the Python user interface for molecular dynamics simulation software package called Nexus
## Description

Nexus is the next-generation molecular dynamics simulation software package to predict protein-ligand binding affinities leading to drug discovery. Nexus is equipped with modern free energy methods, novel sampling algorithms, and parallel programming to support diverse scientific applications.

## Getting Started

### Dependencies

* Python 3.8 or above
* OS: Windows 10, 11
* Visual Studio 2022 for execution via the source code
* (Optional) Anaconda for execution via the published version

### Dependencies

* **Python**: Version 3.8 or above
* **Operating Systems**: 
  * Windows 10, 11
  * Linux (Tested on Ubuntu 20.04 or later)
* **Execution via the source code**:
  * For Windows: Visual Studio 2022
  * For Linux: GCC or Clang (latest version recommended)
* **SWIG**: Required for code generation.
* **Execution the published version**:
  *Anaconda*
* **Python Packages**:
  - `setuptools`, `wheel`, and `cython`

---



### Installation Instructions (Source code)

#### Windows

1. **Install Visual Studio 2022**:
   - Ensure the following components are installed:
     * Desktop development with C++
     * Python development tools
2. **Install SWIG**:
   - Download the SWIG `.zip` file for Windows from the [official website](http://www.swig.org/download.html).
   - Extract the downloaded file to a folder, e.g., `C:\swigwin`.

3. **Add SWIG to the System PATH**:
   - Locate the extracted SWIG folder, e.g., `C:\swigwin`.
   - Open **Environment Variables**:
     1. Press `Win + S` and type **Environment Variables**.
     2. Select **Edit the system environment variables**.
     3. In the **System Properties** window, click **Environment Variables**.
   - Edit the `Path` variable:
     1. In the **System Variables** section, find and select the variable named `Path`.
     2. Click **Edit**.
     3. Click **New** and enter the full path to the SWIG folder (e.g., `C:\swigwin`).
     4. Click **OK** to save changes.
   - Restart any open terminals or your system for the changes to take effect.
   - Verify the installation by running:
     ```bash
     swig -version
     ```
4. **Clone the Repository**:
   - Clone the repository to your local machine.
   - Open the main folder in Visual Studio 2022.

5. **Generate Build Files**:
   - Open the terminal in Visual Studio:
     * Go to `View` → `Terminal` → `Developer PowerShell`.
   - Navigate to the `out/build` directory:
     ```bash
     cd .\out\build\
     ```
   - Run the following command to generate the build files:
     ```bash
     cmake ../../.
     ```

6. **Build the Project**:
   - Open the generated `Nexus.sln` file located in the `out/build` folder.
   - Ensure the solution configuration is set to **Release** (not Debug).
   - In the Solution Explorer:
     * Right-click on `ALL_BUILD` and select `Build`.

7. **Switch Solution Explorer View**:
   - Change the Solution Explorer view to "Folder View":
     * Click the "Switch between solution and available views" button (top-left corner of the Solution Explorer).

8. **Install Required Python Packages**:
   - Open the terminal again and navigate to the main folder.
   - Install the following Python packages:
     ```bash
     pip install setuptools wheel cython
     ```

9. **Build and Install the Nexus Library**:
   - Run the following commands in the terminal:
     ```bash
     python setup.py build -cmingw32
     python setup.py bdist_wheel
     ```
   - Install the Nexus library:
     ```bash
     pip install dist/Nexus-1.0.0-py3-none-any.whl --force-reinstall
     ```

10. **Run the Test File**:
   - Open `NexusTest.py` inside the test folder containing .pdb and system, state .xml files.
   - Adjust the options as desired.
   - Run the test file using the following command:
     ```bash
	 cd .\test\
     python .\NexusTest.py
     ```

#### Linux

1. **Install Required Tools**:
   - Ensure the following packages are installed:
     ```bash
     sudo apt update
     sudo apt install build-essential cmake swig python3-dev python3-pip
     pip3 install --upgrade pip
     ```

2. **Clone the Repository**:
   - Clone the repository to your local machine:
     ```bash
     git clone <repository_url>
     ```
   - Navigate to the main folder:
     ```bash
     cd <repository_folder>
     ```

3. **Generate Build Files**:
   - Navigate to the build directory:
     ```bash
     cd out/build
     ```
   - Run the following command to generate the build files:
     ```bash
     cmake ../../.
     ```

4. **Build the Project**:
   - Compile the project in **Release** mode:
     ```bash
     cmake --build . --config Release
     ```

5. **Install Required Python Packages**:
   - Navigate back to the main folder:
     ```bash
     cd ../../
     ```
   - Install the required Python packages:
     ```bash
     pip3 install setuptools wheel cython
     ```

6. **Build and Install the Nexus Library**:
   - Build the Python package:
     ```bash
     python3 setup.py build
     python3 setup.py bdist_wheel
     ```
   - Install the Nexus library:
     ```bash
     pip install dist/Nexus-1.0.0-py3-none-any.whl --force-reinstall
     ```

7. **Run the Test File**:
   - Open `NexusTest.py` inside the test folder containing .pdb and system, state .xml files.
   - Adjust the options as desired.
   - Run the test file using the following command:
     ```bash
	 cd test
     python3 NexusTest.py
     ```

### Installation Instructions (Published Version)

#### Windows

1. **Install Anaconda**:
   - Download and install **Anaconda** from the [official website](https://www.anaconda.com/).

2. **Create a New Anaconda Environment**:
   - Open the Anaconda Prompt.
   - Create a new environment with Python 3.8 or later:
     ```bash
     conda create -n nexus_env python=3.8
     ```
   - Activate the environment:
     ```bash
     conda activate nexus_env
     ```

3. **Install Dependencies**:
   - While in the `nexus_env` environment, install the required Python packages:
     ```bash
     pip install setuptools wheel cython
     ```
   - Install Nexus version 1.0.0:
     ```bash
     pip install -i https://test.pypi.org/simple/ Nexus==1.0.0
     ```
4. **Prepare the Folder**:
   - Create a folder to contain all necessary files (available in the test folder (github repository)):
     - `Sample_Protein.pdb`: The PDB file for the system.
     - `system_Sample_Protein.xml`: The system XML configuration file.
     - `state_Sample_Protein.xml`: The state XML configuration file.
     - `NexusTest.py`: Test script for Nexus.
   - Copy these files into the folder.

5. **Adjust Options in NexusTest.py**:
   - Open `NexusTest.py` in any text editor or IDE.
   - Modify options such as file paths, configurations, or parameters as needed.

6. **Run NexusTest.py**:
   - Navigate to the folder in the terminal.
   - Run the test script:
     ```bash
     python NexusTest.py
     ```

---

#### Linux

1. **Install Required Tools**:
   - Ensure the following packages are installed:
     ```bash
     sudo apt update
     sudo apt install python3 python3-pip
     pip3 install --upgrade pip
     ```

2. **Install Nexus**:
   - Install Nexus version 1.0.0 using pip:
     ```bash
     pip3 install -i https://test.pypi.org/simple/ Nexus==1.0.0
     ```

3. **Prepare the Folder**:
   - Create a folder to contain all necessary files (available in the test folder (github repository)):
     - `input.pdb`: The PDB file for the system.
     - `system.xml`: The system XML configuration file.
     - `state.xml`: The state XML configuration file.
     - `NexusTest.py`: Test script for Nexus.
   - Copy these files into the folder.

4. **Adjust Options in NexusTest.py**:
   - Open `NexusTest.py` using a text editor (e.g., nano, vim, or VS Code):
     ```bash
     nano NexusTest.py
     ```
   - Modify options such as file paths, configurations, or parameters as needed.

5. **Run NexusTest.py**:
   - Navigate to the folder in the terminal:
     ```bash
     cd /path/to/your/folder
     ```
   - Execute the test script:
     ```bash
     python3 NexusTest.py
     ```

---

### Key Notes
- Ensure the `.pdb`, `.xml`, and `NexusTest.py` files are all in the same folder for convenience.
- Adjust any parameters in the `NexusTest.py` file to match your use case (e.g., file paths, simulation parameters, etc.).
- Output and results will depend on the configuration of the provided `.xml` and `.pdb` files.


## Expected output!

a simple print of:
- "..."

## Authors

Contributors names and contact info
- [Omid Jahanmahin](https://github.com/ozj1)
- [Taisung Lee](https://github.com/taisung)

## Version History

* 1.0.0
    * Initial Release

## License

## Acknowledgments

