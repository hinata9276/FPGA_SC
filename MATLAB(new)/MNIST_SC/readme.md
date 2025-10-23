Steps to reproduce the MATLAB simulation

1) Extract the zip files:

   -> 'multiwaitbar.m' is for showing progress (to lessen the programmer's anxiety)

   -> 'schedule_MNISTnet2_trimmedOpClk' stored the MNIST CNN information in the SC domain

   -> 'SCclassifyMNIST2_runlength_reluState_256' stored the previously simulated output result. You could unzip and load it into MATLAB to examine the resultant bitstream in the dimension of (class, runlength, batch). Summing the runlength dimension gives the total for each class, and the class with the maximum value is the predicted output for each batch. You could also examine the early convergence from the runlength; in most cases, 32 clocks already have the output converged, i.e. cut off computation at 32 clocks and still gives an accurate prediction!

   -> skip the gz files as they will be automatically unzipped in the simulation script

3) Open the SC_MNIST_Sim.mlx. You may need to change the associated directories before running the live script.

4) All the SC functions are contained in the Appendix section in the live script. 

Expected results:

Binary CNN accuracy = 0.9836

Stochastic Computing CNN accuracy = 0.9826 (0.1% accuracy degradation!)
