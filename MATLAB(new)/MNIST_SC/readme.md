# Steps to reproduce the MATLAB simulation

1) Extract the zip files:

   -> 'schedule_MNISTnet2_trimmedOpClk' stored the MNIST CNN information in the SC domain

   -> 'SCclassifyMNIST2_runlength_reluState_256' stored the previously simulated output result. You could unzip and load it into MATLAB to examine the resultant bitstream in the dimension of (class, runlength, batch). Summing the runlength dimension gives the total for each class, and the class with the maximum value is the predicted output for each batch. You could also examine the early convergence from the runlength; in most cases, 32 clocks already have the output converged, i.e. cut off computation at 32 clocks and still gives an accurate prediction! The entire CNN could be computed in just 32 clocks in the SC domain because the bitstream did not wait for accumulation, passing down to the next layer immediately, which is not possible in binary computing!

   -> skip the gz files as they will be automatically unzipped in the simulation script

3) Open the SC_MNIST_Sim.mlx main script. You may need to change the associated directories before running the live script. (only .mlx file for now, have to convert to .m file if need to open from VS Code or other apps)

4) All the SC functions are contained in the Appendix section in the live script.

## File description

- 'multiwaitbar.m' is for showing progress and estimate completion time (to lessen the programmer's anxiety)

- 'Sc Edt Mnist Cm-1.mp4' is a video of the results compiled into frames of a confusion matrix. For every clock cycle, it converges to the actual class.

- 'MNISTnet_small2.mat' is the MNIST CNN in binary, to be used by CPU as ground reference computation.

- 'processLabelMNIST.mat" is a function needed to process MNIST dataset.

## Expected results:

Binary CNN accuracy = 0.9836

Stochastic Computing CNN accuracy = 0.9826 (0.1% accuracy degradation with noisy LFSR!)

## How does it simulate MUX operation efficiently?

In MATLAB, a clock-for-clock simulation will be too slow and time-consuming. Instead, the vector computing technique is employed. It lists all selected bitstreams in rows, then takes the diagonal element of the matrix, which is functionally equivalent to a multiplexing operation across the entire runlength at once.

![image](https://raw.githubusercontent.com/hinata9276/FPGA_SC/refs/heads/main/MATLAB(new)/MNIST_SC/images/vectorComputing2.png)

## Decoding the SC CNN object

The binary CNN model has been custom-converted into an "MUXSchedule" object that can be read and executed in the SC domain. Note that the IO of '-1' and '0' correspond to bias and zero, respectively. The index will be offset by two in the simulation, and the IO kernel will be expanded with 'ones' and 'zeros' arrays to match the indexing. The MUXSize encodes the types of SC hardware to be implemented for a specific kernel.

![image](https://raw.githubusercontent.com/hinata9276/FPGA_SC/refs/heads/main/MATLAB(new)/MNIST_SC/images/AppendixF.jpg)

![image](https://raw.githubusercontent.com/hinata9276/FPGA_SC/refs/heads/main/MATLAB(new)/MNIST_SC/images/scheduleMap.jpg)

![image](https://raw.githubusercontent.com/hinata9276/FPGA_SC/refs/heads/main/MATLAB(new)/MNIST_SC/images/muxSizeEncoding.jpg)

There are several levels of the MACFG (MAC function generator) signals, S0 to S4. Each corresponds to the conditional probability of the MUX select input.

![image](https://raw.githubusercontent.com/hinata9276/FPGA_SC/refs/heads/main/MATLAB(new)/MNIST_SC/images/encodedWeights.jpg)

## How does it perform ReLU activation while you still have no idea of the exact magnitude in the SC domain?

Many SC research samples the bitstream back to binary data to do ReLU activation, but that would be too late since it has already incurred compute latency. A special ReLU function is needed to be executed exclusively in the SC domain, i.e., to perform activation while the data is still in the probability domain. So, how the heck do you know if the bitstream is of positive value while you still have no idea of the final magnitude?

The solution is quite elegant, i.e., exploiting biases in random distributions. If the summation is biased to negativity, then it is deemed to be zero. But how can this logic be translated into hardware architecture? [Please refer to my paper for detailed explanation.](https://doi.org/10.34133/research.0307)

![image](https://raw.githubusercontent.com/hinata9276/FPGA_SC/refs/heads/main/MATLAB(new)/MNIST_SC/images/BReLU.jpg)
