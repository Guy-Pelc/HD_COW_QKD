## Installation

To install the package, please download the entire repository and unzip to a preferred directory.

Before running the program, one must first install:
- Yalmip (https://yalmip.github.io/download/) optimization package for Matlab, and add it to the MATLAB path using
    ```
	addpath(genpath('YALMIP_DIRECTORY'))
    ```
- Mosek (https://docs.mosek.com/10.0/toolbox/install-interface.html) solver package for Matlab, and add it to the MATLAB path using
    ```
	addpath(genpath('MOSKE_DIRECTORY'))
    ```
## Usage
- Set the desired parameters in main.m and run to get the secure key rates, further details for usage are in the main.m file.

- The figures for the lower bounds in the paper can be reconstructed by running gen_figs.m