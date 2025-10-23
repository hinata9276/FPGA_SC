# Steps to reproduce the MATLAB simulation

1) Extract the zip files:

   -> 'multiwaitbar.m' is for showing progress (to lessen the programmer's anxiety)

   -> 'schedule_MNISTnet2_trimmedOpClk' stored the MNIST CNN information in the SC domain

   -> 'SCclassifyMNIST2_runlength_reluState_256' stored the previously simulated output result. You could unzip and load it into MATLAB to examine the resultant bitstream in the dimension of (class, runlength, batch). Summing the runlength dimension gives the total for each class, and the class with the maximum value is the predicted output for each batch. You could also examine the early convergence from the runlength; in most cases, 32 clocks already have the output converged, i.e. cut off computation at 32 clocks and still gives an accurate prediction! The entire CNN could be computed in just 32 clocks in the SC domain because the bitstream did not wait for accumulation, passing down to the next layer immediately, which is not possible in binary computing!

   -> skip the gz files as they will be automatically unzipped in the simulation script

3) Open the SC_MNIST_Sim.mlx main script. You may need to change the associated directories before running the live script.

4) All the SC functions are contained in the Appendix section in the live script.

5) 'Sc Edt Mnist Cm-1' is a video of the results compiled into frames of a confusion matrix. For every clock cycle, it converges to the actual class.

## Expected results:

Binary CNN accuracy = 0.9836

Stochastic Computing CNN accuracy = 0.9826 (0.1% accuracy degradation!)

## How does it simulate MUX operation efficiently?

In MATLAB, a clock-for-clock simulation will be too slow and time-consuming. Instead, the vector computing technique is employed. It lists all selected bitstreams in rows, then takes the diagonal element of the matrix, which is functionally equivalent to a multiplexing operation across the entire runlength at once.


