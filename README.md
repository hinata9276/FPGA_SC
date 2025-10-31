# FPGA_SC
This library contains the resources for Stochastic Computing (SC) implementation in FPGA.

Note: Just uploaded the simulation code in MATLAB(New). You could simulate MNIST CNN classification in the SC domain. Still updating the main readme...

Before continue reading, this is a video of demonstrating a multiplexer (MUX) can be an AI neuron. If you find it interesting, you are in the right place.

[![Watch the video](https://img.youtube.com/vi/hdA5uVT7YiA/0.jpg)](https://www.youtube.com/watch?v=hdA5uVT7YiA)

# What is Stochastic Computing?
Stochastic computing (SC) is another computing domain in contrast to ubiquitous binary computing. Unlike binary computing, SC exploits probability mathematics to perform calculations with a single logic gate.
It was proposed in the 1960s, when combinational logic gates for binary computing were expensive; SC was abandoned as binary computing became feasible and efficient with silicon chip technology. Nowadays, the rise of SC is primarily due to AI and edge computing applications, such as image processing and CNN inference, where SC can provide alternatives for more efficient computation.

Do not confuse with quantum computing; SC uses probabilistic sampling to compute, but the bit information itself still follows the nature of a digital binary bit. Of course, there are give-and-takes in the SC. It is only suited for some applications, especially for image processing and CNN.

Do check out my following paper for a more detailed review and research regarding SC, especially in the CNN use case. Citations are welcome.

https://doi.org/10.1109/ACCESS.2025.3539986 (Toward Universal Multiplexer Multiply-Accumulate Architecture In Stochastic Computing, 2025)

https://doi.org/10.34133/research.0307 (Stochastic Computing Convolution Neural Network Architecture Reinvented For Highly Efficient Artificial Intelligence Workload on Field Programmable Gate Array, 2024)

https://doi.org/10.1007/978-981-16-8129-5_94 (Novel FPGA-Optimized Stochastic Number Generator for Stochastic Computing, 2022)

https://doi.org/10.7717/peerj-cs.309 (Stochastic computing in convolutional neural network implementation: a review, 2020)


# FPGA library in SC is far lacking and unoptimized!
I am personally doing research on implementing SC elements in FPGA, and not much open source can be found in this field of study. Most of the SC research is biased to Application Specific Integrated Circuit (ASIC) design. Thus, I made this library for personal collection in the process of opimizing SC for FPGA.

# Development environment
Software:
- Xilinx Vivado HLx 2019 / Vitis HLS 2020
  - Additional libraries: none
- MATLAB 2020b (for simulation, import and export, and analysis)
  - Deep Learning Toolbox (for binary CNN)
  - Parallel Computing Toolbox (for SC CNN theoretical simulation)
  - Add-On: Deep Learning Toolbox Importer for Tensorflow-Keras Models, MATLAB support for MinGW-w64 C/C++ Compiler (use add-on manager to get them)
  - Community add-on: multiwaitbar, RouletteWheelSelection, 

Hardware:
- Xilinx FPGA development board (I am using Digilent Zybo (Xilinx Zynq Z7010 FPGA SoC), and official Xilinx Kintex7 KC705 dev kit)
- PC:
  - More cores are better, stochastic computation simulation will use up all parallel workers for vector computing!
  - 64GB RAM or above are recommended for large-scale SC simulation.

## Currently developed components
1) FPGA Hardware:
- Stochastic Number Generator (SNG)
  - Weighted Binary Generator (WBG)(4-bit and 8-bit)(original ASIC transcoded logic circuit)
    - LFSR (circular shifting, permuted pair, reset timing injection, EDT capability)
      - WBG frontend partial sharing(M. Yang, B. Li, et.al., DOI:10.1109/ISVLSI.2018.00037) 
      - permutated pair output(Salehi, DOI:10.1109/TVLSI.2019.2963678)
    - WBG backend
  - Weighted Binary Converter (WBC)(4-bit and 8-bit)(novel FPGA-optimized implementation)(Lee, Abdul Halim, Ab Wahab, DOI: 10.1007/978-981-16-8129-5_94)
    - LFSR + WBC with permutated pair output
    - MUX SNG
- Input MUX SNG array with serial interface
- SC CNN Conv layer 1 (MAC-FG, SC MAC, PC, Novel BReLU)
- SC CNN Conv layer 2 (MAC-FG, SC MAC, PC, Novel BReLU)
- SC CNN Conv layer 3 (MAC-FG, SC MAC, PC, Novel BReLU)
- SC CNN Conv layer 3 counter with serial interface
- IP integration (implemented successfully @ 90MHz FPGA)
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

2) Software / simulation:
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
