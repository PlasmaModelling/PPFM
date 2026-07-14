<p align="center">
  <img src="logo.png" alt="PPFM logo" width="200"/>
</p>

<h1 align="center"> Plasma Properties For Many </h1>

**PPFM** is a modular and scalable tool designed in C++ for computing the thermodynamic and transport properties of complex plasma mixtures. 

**PPFM** can also compute and print : 

 - Species Partition Functions from spectroscopic data
 
 - Equilibrium and 2T composition
 
 - Thermodynamic LTE and NLTE 2T properties
 
 - Transport cross sections up to the 4th order (l=1..4) from interaction potentials and other data-sources.
 
 - Interaction collision integrals up to the 4th order needed for the Chapman-Enskog expansion (l=1..4 , s = l..2lmax-l)
 
 - Transport LTE and NLTE 2T properties.

The building is done using **CMake**.

## Setup

Start by **cloning this repository** the way you prefer from **git** or **integrated IDE supports**.

## Prerequisites

Before building PPFM, **ensure the following tools are installed** on your system.

### 🐧 Ubuntu / Linux

Install required tools via terminal through the following command lines:

```bash
sudo apt update
sudo apt install g++ cmake ninja-build libcurl4-openssl-dev libomp-dev python3
```

to verify installations you can always type 
```bash
toolName --version
```

on your terminal.

### 🪟 Windows (much more painful)

1. Install **Microsoft Visual Studio Build Tools** from: https://visualstudio.microsoft.com/visual-cpp-build-tools/ .
1. a) Double-click on program installer and go for MSVC install,
1. b) **CHECK** Desktop development with C++
1. c) **CHECK** From the bar appearing to the right of the installer window : 
    - MSVC Latest Version available
    - Windows 10 SDK if you use Windows 10, 11 for Windows 11, both if you don't know. For Windows 7, Xp and previouses just close the PC and go walking outside questioning yourself about your life choices.
    - C++ CMake Tools for Windows (very important)
    - Reboot the PC.
2. Download and install Python 3.xx from https://www.python.org/downloads/windows/ . 
  - When installing **CHECK** **"Add Python to Path"** (very important)
3. (I think this is not mandatory but try if the program won't work without). Install Ninja, download a zip (recommended) from here: https://github.com/ninja-build/ninja/releases and place it in a folder where you're sure it won't be deleted (write down the folder path).
  - Add Ninja to path by pressing Win+R then digit SystemPropertiesAdvanced, go to "Ambient variables" and click on "System variables". Find a voice named path, click on it then on modify. Now add C:\YourPersonalPathTo\Ninja

## Running a calculation

To run a calculation a main file has to be written, and represent the user environment for the usage of PPFM.
In there, you can include the header files you need for your computations, initialize variables, classes, and use their methods to perform the calculations you need with (in principle) little knowledge of C++ syntax. 

Output classes are provided for writing of computed properties.

See the demo files in the demo folder, or documentation for further details. 

Next, you will need to configure and build the file you implemented as follows

## Build Instructions

Building and launching can be done through terminal command line in the main folder of the program. 

### 🐧 Ubuntu / Linux - command line

The default build will have the ```main.cpp``` file in the main folder of PPFM as its target.
You can do it by command line with 

```bash
cmake -B build -S .
cmake --build build
```

This will generate an executable named ```main.out``` in the ```/executables``` folder.
Then, to execute, you just have to run that file. 
You can do it in the executables folder 

```bash
./main.out
```

or from the main one

```bash
./executables/main.out
```

When building, you can designate the main file to be build by specifying a ```-DMAIN``` argument, the main file can be labeled as you prefer, as long as it contain an int main() {} section

```bash
cmake -B build -S . -DMAIN=Your/Path/To/MainFileName.cpp
cmake --build build
```

For reason of simplicity, the generated executable will be labeled as ```MainFileName.out``` and put in the ```/executables``` folder.

In case of printing of properties, output will be put on an ```/out``` folder

### ⌨️ Ubuntu / Linux - bash scripts

Some bash .sh files are provided to give PPFM a quick try. They can be executed from a terminal into PPFM source folder :   

```bash
./BuildAndRun.sh
```

will build the default main.cpp into PPFM source folder, generate and run the related executable in the executables folder.

You can also specify an argument 

```bash
./BuildAndRun.sh -a demo/mainPerformances.cpp
```

Another bash script in the demo folder 

```bash
./demo/DemoBuilding.sh
```

will build all the demos in the demo folder and generate the related executables in the executables folder

```bash
./executables/RunAll.sh
```

will run all the executables in the executables folder.

```bash
./docs/GenerateDocumentation.sh
```

will generate the doxygen documentation of PPFM in html and latex formats.

You should be able to run the .sh scripts also by right clicking on them and clicking on "Run as a program". 

### ⌨️ Ubuntu / Linux - IDE !

IDE that supports CMake can be used to build, execute and debug the code much easier! Development and testing has been carried out through **VisualStudio Code** with **CMake, C++ and Python extension** installed. 
VScode extensions usually don't install things on your PC, they just tell Code how to behave in some cases. So, remember to satisfy prerequisites for Ubuntu described before. 
Building, and, Release and Debug activities can be run through the cmake toolbar appearing in below the VScode GUI but remember only default option i.e. the main.cpp file in the main folder, will be compiled and executed. 

### 🪟 Windows - command line

1. After installing Visual Studio you'll be able to find **Developer Command Prompt for VS** searching for it in the Win menù.
2. Place yourself in the project main directory by calling:
  - ```bash cd C:\YourPersonalPathTo\PPFM ```
3. Configure (just the first time, minor modifications can usually simply re-build) with the command:
  - ```bash cmake -B build -S . ```
4. Build with the command:
  - ```bash cmake --build build ```
5. rooster.exe will be created in the mainfolder, just double-click on it.

### 🪟 Windows - IDE
The same configuration of VScode used for Linux came very handy in addressing the code compiling and execution. Unfortunally, debugging is still disabled. If you're used in working on debugging you should be able to implement your personal launch.json and properties.json and run debug sessions yourself.

If you familiar with Windows development and want to improve PPFM capabilities please contact me I'll be truly grateful.

Remember, your favourite IDE has to support CMake to run the program.

#### Notes on building : 

PPFM uses a compile-time selection of chemical species through a generated `std::variant` type (`AcceptedSpecies.h`).  
Each executable is associated with a specific set of species inferred from the selected main source file, this is the reason every main file needs its own configure and building to generate an executable.

## 📂 Directory Structure

For being in its very first release PPFM handle data paths by hard-coding them into strings and macros. Further Development will see better data handling. 
Please, do not rename any existing folder, do whatever inside the out/ folder for further data processing. PPFM also include routines to sistematically print data in desired folders. Check demo mainArH2.cpp for an example. 

```text
PPFM/
├── src/              # Source code (.h/.cpp) + 1 python .py script 
  ├── src/alglib/     # Internal Alglib static library
├── data/             # Input CSV data
  ├──Collision_Integrals 
  ├──Differential_Cross_Sections
  ├──Electronic_Configurations
  ├──Momentum_Transfer_Cross_Section
  ├──Partition_Functions
  ├──Phase_Shifts
├── demo/             # Some main.cpp files for demo running
├── docs/
├── build/            # Build directory (ignored by git)
├── executables/      # Executables will be placed here
├── CMakeLists.txt    # Project configuration
├── main.cpp          # The default environment for users, default "-DMAIN=" argument to build.
├── README.md
├── LICENSE.md
├── AUTHORS.md

```

## Author disclaimer

I'm proud of this project but it is still full of flaws due to flaws in my education and skills, so please read carefully the followings.

Basically, PPFM saw light to gather toghether ab-initio calculation routines for facing this painful problem of determining plasma transport properties. That is, to the author knowledge, it is the only tool up to date that actually include routines to compute deflection angles, transport cross sections and collision integrals required for the Chapman-Enskog approximation. 

If you're just interested in computed properties for your particular case, I'm sorry to disappoint you but you won't find other data apart from those used for demos. 
For now. 
BUT, PPFM offers the most scalable and powerful environment to compute them from data gathered from the scientific literature or many different data-sources, do not surrender here! Computing Transport Properties for plasmas is complex but its' never been this easy either!

AND, if you're interested in: computed properties, collaborating, PPFM development, becoming a developer, helping, sharing knowledge or you have any kind of questions related to this work, please contact me, I'll be very keen to help.

alberto.vagnoni3@unibo.it

albertovagnoni97@gmail.com

## 🤝 Contributing

Contributions are welcome!

If you'd like to contribute to PPFM:

1. Fork this repository
2. Clone your fork locally
3. Make your changes (on `master` or a separate branch)
4. Commit and push
5. Open a Pull Request to `master` on this repository

Please ensure your contributions are well-documented and tested. All changes will be reviewed and must be approved by the owner before being merged.

## Authors

- **Alberto Vagnoni** — University of Bologna, Italy  
- **Emanuele Ghedini** — University of Bologna, Italy  
- **Paolo Sanibondi** (past contributor, former Research Assistant University of Bologna, Italy)

## Scientific roots

This project has its roots in the 2008-2011 activities performed by some of the contributors in the Group for Industrial Applications of Plasmas, led by Prof. Vittorio Colombo, whose contributions to plasma modelling and transport theory have been an important source of inspiration for the present developments.

The authors gratefully acknowledge this scientific influence on the theoretical and methodological foundations of the present work.

## License

This project is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

<a href="https://creativecommons.org/licenses/by/4.0/">
PPFM © 2025 by Emanuele Ghedini, Alberto Vagnoni (University of Bologna, Italy)
</a> is licensed under 
<a href="https://creativecommons.org/licenses/by/4.0/">
Creative Commons Attribution 4.0 International
</a>
<img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="CC" style="max-width: 1em;max-height:1em;margin-left: .2em;">
<img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="BY" style="max-width: 1em;max-height:1em;margin-left: .2em;">