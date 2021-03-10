# FPGA_SC
This library contains the resources for Stochastic Computing (SC) implementation in FPGA.

# What is Stochastic Computing?
Stochastic computing (SC) is another computing domain in contrast to ubiquitous binary computing. Unlike binary computing, SC exploits the probability mathematics to perform calculation with single logic gate.
It was proposed in the 1960s when the combinational logic gates for binary computing was expensive at that time, then SC was abandoned when the binary computing becomes feasible and efficient with silicon chip technology. Nowadays, the bring up is mostly due to the AI compute bottlenecking and edge couputing application such as image processing and signal processing, in which SC could provide the alternatives for more efficient coumputation in those field of applications.

It seems like quantum computing but it is not, it just make use of probabilistic sampling to compute, but the bit information itself still follows the nature of binary bit. Of course there are pros and cons in the SC. It is only suited for some applications, especially for image processing and CNN.

Do check out my paper (Lee, Abdul Halim, DOI:10.7717/peerj-cs.309) for more information. Citations are welcomed.

# FPGA library in SC is far lacking!
I am personally doing research on implmeneting SC elements in FPGA, and not much open source can be found in this field of study. Thus, I made this library for personal collection as well as looking for improvement through the power of open source.

# Development environment
- Xilinx Vivado HLx 2019
- Xilinx FPGA development board (I am using Zybo)
- MATLAB 2020b (for simulation, import and export, and analysis)
  - Deep Learning Toolbox
  - Parallel Computing Toolbox
  - Add-On: Deep Learning Toolbox Importer for Tensorflow-Keras Models, MATLAB support for MinGW-w64 C/C++ Compiler
- PC:
  - more core is better, stochastic computation simulation will use parallel workers!
  - 64GB RAM and above are recommended for simulation!

# Currently developed components
1) Hardware, 
Note: LFSR (Linear Feedback Shift Register), MUX (Multiplexer)

- Stochastic Number Generator (SNG)
  - Weighted Binary Generator (WBG)(4-bit and 8-bit)(original ASIC transcoded logic circuit)
    - LFSR + WBG frontend(M. Yang, B. Li, et.al., DOI:10.1109/ISVLSI.2018.00037) with permutated pair output(Salehi, DOI:10.1109/TVLSI.2019.2963678)
    - WBG backend
  - Weighted Binary Converter (WBC)(4-bit and 8-bit)(novel FPGA-optimized implementation)(Paper accepted)
    - LFSR + WBC with permutated pair output
    - MUX SNG

2) Software / simualtion
- LFSR and WBG simulation
  - Simulate LFSR and WBG (circular shift, peermute)
- MUX SUC CNN multi-layer computation simulation
  - compute CNN from 1 layer down to 5 layers, both in binary domain and stochastic domain, only linear, ReLU and Tanh activation function available.
  - only 3x3 kernel weight is supported, only simulate on one weight. 
- MUX SUC CNN weight scheduler (requires at least MATLAB version 2020b, older than that cannot work.)
  - Import HDL5 CNN model and compute the timing for MUX SUC
  - Analyze the MUX timing requirement and optimize timing to reduce required LFSR.
  - Export the timing and MUX I/O allocation in XML file, to be used in Xilinx Vivado HLS C++ compiler.

# Components under development
- Shifted Unary Coded (SUC) Adder
- MUX-based SUC Adder, 8-bit resolution
- 
