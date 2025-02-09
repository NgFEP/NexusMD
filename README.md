# NexusMD

A demo version of the Python user interface for molecular dynamics simulation software package called NexusMD
## Description

NexusMD is the next-generation molecular dynamics simulation software package to predict protein-ligand binding affinities leading to drug discovery. NexusMD is equipped with modern free energy methods, novel sampling algorithms, and parallel programming to support diverse scientific applications.

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
   - Ensure the following components are installed [Link](https://visualstudio.microsoft.com/downloads/):
     * Desktop development with C++
     * Python development package 
2. **Install SWIG**:
   - Download the swigwin latest version `.zip` file for Windows from the [official website](http://www.swig.org/download.html).
   - Extract the downloaded file to a folder, e.g., `C:\Users\[YourUsername]\swigwin-4.3.0`.

3. **Add SWIG to the System PATH**:
   - Locate the extracted SWIG folder containing **swig.exe**, e.g., `C:\Users\[YourUsername]\swigwin-4.3.0`.
   - Open **Environment Variables**:
     1. Press `Win + S` and type **Environment Variables**.
     2. Select **Edit the system environment variables**.
     3. In the **System Properties** window and **Advanced** tab, click **Environment Variables**.
   - Edit the `Path` variable:
     1. In the **System Variables** section, find and select the variable named `Path`.
     2. Click **Edit**.
     3. Click **New** and enter the full path to the SWIG folder containing **swig.exe** (e.g., `C:\Users\[YourUsername]\swigwin-4.3.0`).
     4. Click **OK** to save changes.
   - Restart any open terminals or your system for the changes to take effect.
   - Verify the installation by running:
     ```
     Windows key + R (to open the Run dialog), then type "cmd" and press Enter to open command prompt, type
     swig -version
     ```

3. **Install CUDA Toolkit and Add to System PATH**:
   
   - Download **CUDA Toolkit** from [NVIDIA's official site](https://developer.nvidia.com/cuda-downloads).
   - Choose the appropriate version for your **GPU and operating system**.
   - Install it following the **on-screen instructions**.
   
   - **Add CUDA binaries to the system's PATH**:
     1. Press `Win + S`, type **Environment Variables**, and open it.
     2. In **System Properties**, click **Environment Variables**.
     3. Under **System variables**, find and edit the **Path** variable.
     4. Click **New**, then add the following paths (adjust for your CUDA version):
        ```
        C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin
        C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\libnvvm\bin
        C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\libnvvm\lib
        ```
     5. Click **OK** to save changes.

   - **Verify CUDA Installation**:
     ```bash
     nvcc --version
     ```
     ```bash
     nvidia-smi
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
   - Open the generated `NexusMD.sln` file located in the `out/build` folder.
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

9. **Build and Install the NexusMD Library**:
   - Run the following commands in the terminal:
     ```bash
     python setup.py build -cmingw32
     python setup.py bdist_wheel
     ```
   - Install the NexusMD library:
     ```bash
     pip install dist/NexusMD-1.0.0-py3-none-any.whl --force-reinstall
     ```

10. **Run the Test File**:
   - Open `NexusMDTest.py` inside the test folder containing .pdb and system, state .xml files.
   - Adjust the options as desired.
   - Run the test file using the following command:
     ```bash
	 cd .\test\
     python .\NexusMDTest.py
     ```

### Installation Instructions (Published Version)

#### Windows

1. **Install Anaconda**:
   - Download and install **Anaconda** from the [official website](https://docs.anaconda.com/anaconda/install/).

2. **Create a New Anaconda Environment**:
   - Press `Win + S` and type **Environment Variables**.
   - Open the Anaconda Prompt.
   - Create a new environment with Python 3.11.5:
     ```bash
     conda create -n nexusMD_env python=3.11.5
     ```
   - Activate the environment:
     ```bash
     conda activate nexusMD_env
     ```
3. **Install NexusMDMD library**:
   - run the following command to install NexusMDMD 1.0.0 MD software package:
     ```bash
     pip install -i https://test.pypi.org/simple/ NexusMD
     ```
4. **Prepare the Folder**:
   - Create a folder (preferebly on your Desktop for easy navigation on Anaconda Prompt) to contain all necessary files (available in the test folder (NexusMD github repository)):
     - `Sample_Protein.pdb`: The PDB file for the system.
     - `system_Sample_Protein.xml`: The system XML configuration file.
     - `state_Sample_Protein.xml`: The state XML configuration file.
     - `NexusMDTest.py`: Test script for NexusMD.
   - Copy these files into the folder.

5. **Adjust Options in NexusMDTest.py**:
   - Open `NexusMDTest.py` in any text editor or IDE. A straightforward option's [Notepad++](https://notepad-plus-plus.org/downloads/)
   - Uncomment any or all of the following 4 forces inside the `NexusMDTest.py` to enable their calculations during the MD simulation.
     ```bash
     #dispatcher.EnableHarmonicBondForce()
     #dispatcher.EnableHarmonicAngleForce()
     #dispatcher.EnablePeriodicTorsionForce()
     #dispatcher.EnableNonbondedForce()
     ```

6. **Run NexusMDTest.py**:
   - Navigate to the folder in the Anaconda Prompt make sure nexusMD_env environment is activated.
   - Run the test script:
     ```bash
     python NexusMDTest.py
     ```
7. **Expected output!**:
   - terminal
     ```bash
     	Enabling forces...
	Running simulation...
	Harmonic Bond Force is enabled
	Harmonic Angle Force is enabled
	Periodic Torsion Force is enabled
	Nonbonded Force is enabled
	1   Simulation Step #1 is completed
	2   Simulation Step #2 is completed
	3   Simulation Step #3 is completed
	4   Simulation Step #4 is completed
	5   Simulation Step #5 is completed
	6   Simulation Step #6 is completed
	7   Simulation Step #7 is completed
	8   Simulation Step #8 is completed
	9   Simulation Step #9 is completed
	10   Simulation Step #10 is completed
	Simulation completed successfully with selected forces enabled.
     ```
     - generated files: output_SimRun_Sample_Protein.pdb, PVFReportoutput_SimRun_Sample_Protein.pdb (position, velocity and force quantities), TotalEnergyoutput_SimRun_Sample_Protein.pdb (Total Kinetic, Potential, Total Energy (kJ/mol) quantities)
---

### Key Notes
- Ensure the `.pdb`, `.xml`, and `NexusMDTest.py` files are all in the same folder for convenience.
- Adjust any parameters in the `NexusMDTest.py` file to match your use case (e.g., simulation parameters, desired forces to enable, etc.).
- Output and results will depend on the configuration of the provided `.xml` and `.pdb` files.

## Authors

Contributors names and contact info
- [Omid Jahanmahin](https://github.com/ozj1)
- [Taisung Lee](https://github.com/taisung)

## Version History

* 1.0.0
    * Initial Release

## License

## Acknowledgments

