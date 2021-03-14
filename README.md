# FPGA_SC
This library contains the resources for Stochastic Computing (SC) implementation in FPGA.

# What is Stochastic Computing?
Stochastic computing (SC) is another computing domain in contrast to ubiquitous binary computing. Unlike binary computing, SC exploits the probability mathematics to perform calculation with single logic gate.
It was proposed in the 1960s when the combinational logic gates for binary computing was expensive at that time, then SC was abandoned when the binary computing becomes feasible and efficient with silicon chip technology. Nowadays, the bring up of SC is mostly due to the AI compute and edge couputing application such as image processing and CNN inferencing, in which SC could provide the alternatives for more efficient coumputation in those field of applications.

Do not confuse with quantum computing, SC just make use of probabilistic sampling to compute, but the bit information itself still follows the nature of digital binary bit. Of course there are pros and cons in the SC. It is only suited for some applications, especially for image processing and CNN.

Do check out my paper (Lee, Abdul Halim, DOI:10.7717/peerj-cs.309) for more detailed review regarding SC, especially in the CNN use case. Citations are welcomed.

# FPGA library in SC is far lacking and unoptimized!
I am personally doing research on implementing SC elements in FPGA, and not much open source can be found in this field of study. Most of the SC research are biased to Application Specific Integrated Circuit (ASIC) design. Thus, I made this library for personal collection in the process of opimizing SC im for FPGA. Meanwhile, I also look for improvement through the power of open source community.

# Development environment
Software:
- Xilinx Vivado HLx 2019
  - Additional libraries: RapidXML (please download and include it in the compiler).
- MATLAB 2020b (for simulation, import and export, and analysis)
  - Deep Learning Toolbox
  - Parallel Computing Toolbox
  - Add-On: Deep Learning Toolbox Importer for Tensorflow-Keras Models, MATLAB support for MinGW-w64 C/C++ Compiler (use add-on manager to get them)

Hardware:
- Xilinx FPGA development board (I am using Digilent Zybo (Xilinx Zynq Z7010 FPGA SoC), any compatible boards are welcomed to be verified)
- PC:
  - more core is better, stochastic computation simulation will use parallel workers!
  - 64GB RAM or above are recommended for large-scale SC simulation.

## Currently developed components
1) FPGA Hardware:
- Stochastic Number Generator (SNG)
  - Weighted Binary Generator (WBG)(4-bit and 8-bit)(original ASIC transcoded logic circuit)
    - LFSR + WBG frontend(M. Yang, B. Li, et.al., DOI:10.1109/ISVLSI.2018.00037) with permutated pair output(Salehi, DOI:10.1109/TVLSI.2019.2963678)
    - WBG backend
  - Weighted Binary Converter (WBC)(4-bit and 8-bit)(novel FPGA-optimized implementation)(Paper accepted)
    - LFSR + WBC with permutated pair output
    - MUX SNG
- Note: LFSR (Linear Feedback Shift Register), MUX (Multiplexer)

2) Software / simualtion:
- LFSR and WBG simulation
  - Simulate LFSR and WBG (circular shift, permute)
- MUX select random function
  - Simulate MUX select random programming to operate MUX SUC Adder.
- MUX SUC CNN multi-layer computation simulation (still working on it...)
  - compute CNN from 1 layer down to 5 layers, both in binary domain and stochastic domain, only linear, ReLU and Tanh activation function available.
  - only 3x3 kernel weight is supported, only simulate on one weight. 
- MUX SUC CNN weight scheduler (requires at least MATLAB v2020b, older than that cannot work.) (still working on it...)
  - Import HDF5 Keras CNN model, partition and compute the timing for MUX SUC for each weight and bias.
  - Analyze the MUX timing requirement and optimize timing to reduce required LFSR resources.
  - Export the timing and MUX I/O allocation in XML file, to be used in Xilinx Vivado HLS C++ compiler.

## Components under development
1) FPGA Hardware:
- MUX-based Shifted Unary Coded (SUC) Adder, 8-bit resolution
- Stanh (stochastic TanH) FSM
- End counter
- updating...

2) Software / simualtion:
- MUX SUC random function generator (C++)
- updating...

# Usage
The resources contains MATLAB folder and Xilinx Vivado folder. MATLAB codes are mainly used for simulation, data synthesis and analysis, and file import/export, whilst Xilinx C++ codes are mainly used for Xilinx Vivado High Level Synthesis (HLS) to generate HDL and RTL file from C++ code function.

## MATLAB coding
THe MATLAB folder contains collection of simulation and import/export fuctions.
- Random_function.m -> To simulate the programming of converting list of weight to random MUX select input for MUX SUC Adder. More info could be found in the code itself.
- SUC_simulation_nxpx.m -> To simulate the CNN in both binary and stochastic domains, up to 5 layers depth. 
- LFSR_simulation.m -> To simulate the generation of LFSR arrays, including circular shifting and permutation techniques in LFSR sharing scheme.
- rest with .mat format is the data for "SUC_simulation_nxpx" simulation.
- More info (code usage and functions) could be found in the code itself.

## C++ coding
- updating...
