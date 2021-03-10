%This code attempts to automate the scheduler events for any weight input
clear;
w=[48 44 -54; 60 25 -106; 60 -9 -75]; bias = -1;       %weight of n1p1 
%w=[-105 31 69; -108 17 44; -77 8 23] ;bias=-27;        %weight of n2p1
%w=[45 -76 -111; 69 16 -144; 95 113 -14]; bias = -7;    %weight of n2p2
%w=[-70 63 112; -121 -227 -107; 60 62 -53]; bias = 38;   %weight of n3p2 
%w=[-35 42 97; 31 15 -5; 4 -23 -108]; bias = -6;         %new weight 1
%w=[93 65 -52; 100 31 -160;-4 -167 -83]; bias = 14;      %new weight 2
%w=[-75 -217 -114; -187 -2 57; -64 43 125]; bias = 33;     %new weight 3
%w=[49 44 98; 32 60 -19; -80 -117 -87]; bias = -1;       %new weight 4
%w=[2 100 86; -119 14 49; -169 -43 64]; bias = -2;       %new wright 5
%w=[48 58 60; -33 65 50; -131 -117 4]; bias = -3;        %new weight 6

runlength = 510;
sample_size = 600;  %predefined matrix run will rewrite this parameter
LFSR8 = [1 0 0 0 0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loading predefined matrices for algorithm comparison
%dir = 'D:\PhD\design\matlab';
dir = 'C:\Users\Sunny\Desktop\PhD\design\matlab';
%load(strcat(dir,'\sample_matrix.mat'));    %600 samples 
%load(strcat(dir,'\sample_matrix_3k.mat')); %3000 samples
%load(strcat(dir,'\sample_matrix_5400.mat'));   %L1:L2 = 5400:600 samples
%load(strcat(dir,'\sample_matrix_8100.mat'));   %L1:L2:L3 = 8100:900:100
%load(strcat(dir,'\sample_matrix_48600.mat'));  %L1:L2:L3=48600:5400:600
load(strcat(dir,'\sample_matrix_437400.mat')); %L1:L2:L3:L4 = 437400:48600:5400:600
%load(strcat(dir,'\sample_matrix_1968k.mat'));   %L1:L2:L3:L4:L5 = 1968300:218700:24300:2700:300
%comment out above lines if requires random matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
for j = 1:sample_size  %for generating predefined matrix
   x=round(rand(3)*255);
   sample(:,:,j)=x;
end
dir = 'D:\PhD\design\matlab';
%save('C:\Users\Sunny\Desktop\PhD\design\matlab\sample_matrix.mat', 'sample');
save(strcat(dir,'\sample_matrix_437400.mat'), 'sample');
%}
for i = 1:runlength-1 %linear feedback shift register
    LFSR8(i+1,1) = xor(xor(LFSR8(i,8),LFSR8(i,6)),xor(LFSR8(i,5),LFSR8(i,4)));
    LFSR8(i+1,2) = LFSR8(i,1);    LFSR8(i+1,3) = LFSR8(i,2);    LFSR8(i+1,4) = LFSR8(i,3);
    LFSR8(i+1,5) = LFSR8(i,4);    LFSR8(i+1,6) = LFSR8(i,5);    LFSR8(i+1,7) = LFSR8(i,6);
    LFSR8(i+1,8) = LFSR8(i,7);
end
LFSR8_wbg = setup_LFSR(LFSR8, runlength);

%scheduling MUX selection
sche = schedule_mux(w,bias);    %return struct array of n and p
to_extract_n = 2*round(reshape(sum(sche.n(2,:,:)>0), [size(sche.n,3),1])/2);
to_extract_p = 2*round(reshape(sum(sche.p(2,:,:)>0), [size(sche.p,3),1])/2);

n1_f=0; n2_f=0; n3_f=0; p1_f=0; p2_f=0; p3_f=0;
if size(to_extract_n,1) >0
    mux_n1 = sche.n(:,1:to_extract_n(1),1);
    idx = [ceil(mux_n1(1,:)/3); mod(mux_n1(1,:),3)];
    idx(2,idx(2,:)==0)=3; n1_f=1;
    mux_n1 = [idx ; mux_n1(2,:)]
end
if size(to_extract_n,1) >1
    mux_n2 = sche.n(:,1:to_extract_n(2),2);
    idx = [ceil(mux_n2(1,:)/3); mod(mux_n2(1,:),3)];
    idx(2,idx(2,:)==0)=3; n2_f=1;
    mux_n2 = [idx ; mux_n2(2,:)]
end
if size(to_extract_n,1) >2
    mux_n3 = sche.n(:,1:to_extract_n(3),3);
    idx = [ceil(mux_n3(1,:)/3); mod(mux_n3(1,:),3)];
    idx(2,idx(2,:)==0)=3; n3_f=1;
    mux_n3 = [idx ; mux_n3(2,:)]
end
if size(to_extract_p,1) >0
    mux_p1 = sche.p(:,1:to_extract_p(1),1);
    idx = [ceil(mux_p1(1,:)/3); mod(mux_p1(1,:),3)];
    idx(2,idx(2,:)==0)=3; p1_f=1;
    mux_p1 = [idx ; mux_p1(2,:)]
end
if size(to_extract_p,1) >1
    mux_p2 = sche.p(:,1:to_extract_p(2),2);
    idx = [ceil(mux_p2(1,:)/3); mod(mux_p2(1,:),3)];
    idx(2,idx(2,:)==0)=3; p2_f=1;
    mux_p2 = [idx ; mux_p2(2,:)]
end
if size(to_extract_p,1) >2
    mux_p3 = sche.p(:,1:to_extract_p(3),3);
    idx = [ceil(mux_p3(1,:)/3); mod(mux_p3(1,:),3)];
    idx(2,idx(2,:)==0)=3; p3_f=1;
    mux_p3 = [idx ; mux_p3(2,:)]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start sampling data point
plot_response = 0; %to plat intermediate response of each layer
conv2D_L2 = 1;  %0 = no 2nd layer SC neuron, 1 = run 2nd layer SC neuron simulation,
conv2D_L3 = 1;  %0 = no 3rd layer SC neuron, 1 = run 3rd layer SC neuron simulation,
conv2D_L4 = 1;  %0 = no 4th layer SC neuron, 1 = run 4th layer SC neuron simulation,
conv2D_L5 = 0;  %0 = no 5th layer SC neuron, 1 = run 5th layer SC neuron simulation,
                %using the same weight as the main config for simplicity.
                %Timing mode will be forced to SC random programming mode.
mode = 1;               % choose 1 or 2
layer_reschedule = 0;   % choose 0 or 1
%mode:   1) SC random programming timing mode
%        2) MATLAB pure random probability distribution 
%layer_reschedule:  0) use the same MUX schedule for all layers
%                   1) use different MUX schedule for each layer
optimize_sel = 1;   % choose 0 or 1
%use optimized LFSR settings,   0) no optimization
%                               1) use optimization clock for MUX
FSM_bit = 2;    %define the capacity of the FSM memory of 2^n
fraction_cycle_max = 2; %define max fraction cycle, max out to make FSM a Stanh activation
%//////////////////////////////////////////////////////
%{
%Experimental, formula from Kim.et.al. Calculate for nearest multiple of 2 for FSM_add
s = 2/(2^FSM_bit);
n = 2;
r = 2*round(((2*(1-s)*(n-1))/(s*(1-1.835*(2*n)^-0.5552))+2*n)/2);
bipolar_FSM_add = r-2^FSM_bit; %experimental modification to TanH FSM of bipolar input for L2 and L3 Conv
                                %in addition of origial FSM max state
%}
r = 2^(FSM_bit+1);
bipolar_FSM_add = r-2^FSM_bit;
%//////////////////////////////////////////////////////                                
binact = 1;     %layer 1 activation function
binact_L2 = binact;  %layer 2 activation function
binact_L3 = binact; binact_L4 = binact; binact_L5 = binact; 
%Select binary activation function:   
%   1) no activation, linear (currently only can match the SC weight integral
%   with the binary counterpart at good accuracy, still testing solution for pipelining)
%   2) clipped ReLU activation (need to set "fraction_cycle_max")
%   3) tanh activation following Stanh(k,x) = tanh(kx/2) function, where k is
%   the number of FSM states ("fraction_cycle_max" will be auto maxed out)

%   Note
%   1) Linear and tanh activation to be compared in SC bipolar mode, the 
%   clipped ReLU activation is compared in SC unipolar mode.
%   2) Tanh activation is known to be inaccurate in low number of states, 
%   i.e. 1~3 bit FSM. 4 bit and above is more accurate.

%////////////////////////////////////////////////////////////////////////
%error checking commands to ensure codes run as expected
if conv2D_L5 == 1
    fprintf('Code will run through 5 layers convolution simulation\n');
    conv2D_L2 = 1; %force 2nd ~ 4th layer execution.
    conv2D_L3 = 1;
    conv2D_L4 = 1;
elseif conv2D_L4 == 1
    fprintf('Code will run through 4 layers convolution simulation\n');
    conv2D_L2 = 1; %force 2nd ~ 3th layer execution.
    conv2D_L3 = 1;
elseif conv2D_L3 == 1
    fprintf('Code will run through 3 layers convolution simulation\n');
    conv2D_L2 = 1; %force 2nd layer execution.
elseif conv2D_L2 == 1 
    fprintf('Code will run through 2 layers convolution simulation\n');
else
    fprintf('Code will run single layer convolution simulation\n');
end
if mode == 1
    fprintf('SC random programming timing mode for MUX sel.\n');
elseif mode == 2
    fprintf('MATLAB pure random permutation mode for MUX sel.\n');
else
    mode = 1;
    fprintf('Default to SC random programming timing mode for MUX sel.\n');
end
if exist('sample','var')
    fprintf('Using predefined matrices, ');
    sample_size = size(sample,3);
else
    fprintf('Using random generated matrices, ');
    sample_rand = zeros(3,3,sample_size);
end
if optimize_sel == 0 && mode == 1
    fprintf('MUX scheduling unoptimized, ');
elseif optimize_sel == 1 && mode == 1
    fprintf('MUX scheduling is optimized, ');
end
if layer_reschedule == 0
    fprintf('same MUX schedule for all layers.\n');
else
    fprintf('MUX will be rescheduled for each layer.\n');
end
if binact == 3
    fraction_cycle_max = runlength;
    fprintf('Stanh FSM max state = %d, modded bipolar FSM state for L2 and L3 = %d \n', 2^FSM_bit-1, r);
elseif binact == 2
    fprintf('ReLU FSM max state = %d, FSM-1 in every %d SC cycle \n',2^FSM_bit-1,fraction_cycle_max);
end
%disable parfor warning if any
warning_id = 'MATLAB:mir_warning_maybe_uninitialized_temporary';
warning('off',warning_id)
%end of script setup
%//////////////////////////////////////////////////////////////////////////
%                   START OF LAYER 1 CONVOLUTION OPERATION
%//////////////////////////////////////////////////////////////////////////
binary_sample = zeros(3,sample_size);        %array to store binary output
fprintf('Conv2D L1 runlength = %d, samples = %d \n', runlength, sample_size);
for j = 1:sample_size %binary compute here
    if exist('sample','var')
        x=sample(:,:,j);  
    else
        x=round(rand(3)*255); 
    end
    % start of neuron iteration in binary domain
    %w is in SC scaled binary and has to be multiplied with x as fraction of 255
    wf = w/255; %convert weight to fraction of 255
    y_binary = wf.*x;
    y_binary = sum(y_binary(:))+bias;
    %w_binary is in scaled binary and has to bey converted to decimals
    y_binary = y_binary/255;
    y_weight = y_binary;
    
    if binact == 1
        %clipping output
        if y_binary < -1 
            y_binary = -1;
        elseif y_binary > 1
            y_binary = 1;
        end
    elseif binact == 2 %clipping binary output to SC range limit
        if y_binary < 0 %ReLU
            y_binary = 0;
        elseif y_binary > 1 %clipped ReLU
            y_binary = 1;
        end
    elseif binact == 3
        %it is unknown why the SC did not exactly follow the Stanh function
        %it could be the unipolar nature of the initial input, or it could
        %be the small number of parallel counting input.
        %somehow the behaviour has to be proven mathematically.
        if FSM_bit < 5
            y_binary = tanh((2^FSM_bit)*y_binary); 
        else
            y_binary = tanh((2^FSM_bit)*y_binary/2);
        end
    end
    %shift and scale binary to compare SC in bipolar mode, clipped ReLU
    %will remain the same as the original output
    if binact == 1 || binact == 3    
        y_binary_s = (y_binary+1)/2;
    else
        y_binary_s = y_binary;
    end
    
    binary_sample(1,j)= y_binary;   %activated output
    binary_sample(2,j)= y_weight;   %original convolution output
    binary_sample(3,j)= y_binary_s; %scaled output to match stochastic output if any
    % end of neuron iteration in binary domain
end
fprintf('Binary computing ended, proceeding to SC.\n');
%////////////////////////////////////////////////////////////////////////////////////////////
% start of neuron iteration in stochastic domain    
%scheduler strategy 5: SC random programming mode
stochastic_sample = zeros(2,sample_size);    %store SUC adder (mux based) sample data points
SUC_record1 = cast(zeros(1,runlength,sample_size),'int16');    %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
SUC_record2 = cast(zeros(1,runlength,sample_size),'int16'); 
SUC_record3 = cast(zeros(1,runlength,sample_size),'int16'); 
SUC_record = cast(zeros(3,runlength,sample_size),'int16');
if n1_f
    sel_n1 = gen_sel(mux_n1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n1');
end
if n2_f
    sel_n2 = gen_sel(mux_n2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n2');
end
if n3_f
    sel_n3 = gen_sel(mux_n3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n3');
end
if p1_f
    sel_p1 = gen_sel(mux_p1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p1');
end
if p2_f
    sel_p2 = gen_sel(mux_p2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p2');
end
if p3_f
    sel_p3 = gen_sel(mux_p3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p3');
end
fprintf('Emulating stochastic computing of Conv_L1');
%testing on parallel processing to speed up emulation, but the actual
%response is not guaranteed to be faster than single thread, depending on
%the number of workers available in the local machine.
%Could be changed to single-threaded for() loop at any time.
tic
%for j = 1:sample_size %stochastic compute here
parfor j = 1:sample_size %stochastic compute here
    x=sample(:,:,j);  %in parfor loop, exist() cannot be used
    %{
    if exist('sample','var')
        x=sample(:,:,j);  
    else
        x=round(rand(3)*255); 
    end
    %}
    if n1_f
        mux_in_n1 = gen_feed_bin(x, mux_n1, sel_n1, runlength);
    end
    if n2_f
        mux_in_n2 = gen_feed_bin(x, mux_n2, sel_n2, runlength);
    end
    if n3_f
        mux_in_n3 = gen_feed_bin(x, mux_n3, sel_n3, runlength);
    end
    if p1_f
        mux_in_p1 = gen_feed_bin(x, mux_p1, sel_p1, runlength);
    end
    if p2_f
        mux_in_p2 = gen_feed_bin(x, mux_p2, sel_p2, runlength);
    end
    if p3_f
        mux_in_p3 = gen_feed_bin(x, mux_p3, sel_p3, runlength);
    end
    %///////////////////////////////////until here, will continue
    %Schedule selector for the first loop to save simulation runtime.
    %This scheduler coding might get lengthy, using while loop for code
    %folding.
    %end of scheduler strategy 5
    %///////////////////////////////////////////////////////////////////////////////////////////
    %start of SC cycle
    %settings for the SC output
    SUC_sum = 0;                %reset SUC sum
    fraction_cycle = 0;         %reset fraction cycle
    FSM_max = 2^FSM_bit-1;      %define FSM max state, must be 2^n-1
    FSM = round(FSM_max/2);     %reset FSM state
    FSM_half = FSM-1;
    for i = 1:runlength
        temp1 = []; temp2 = []; temp3 = [];
        %start computing     
        toAND = LFSR8_wbg(i);
        if n1_f
            SUCn_1 = bitand(mux_in_n1(i),toAND) > 0;
        else
            SUCn_1 = 0;
        end
        if n2_f
            SUCn_2 = bitand(mux_in_n2(i),toAND) > 0;
        else
            SUCn_2 = 0;
        end
        if n3_f
            SUCn_3 = bitand(mux_in_n3(i),toAND) > 0;
        else
            SUCn_3 = 0;
        end
        if p1_f
            SUCp_1 = bitand(mux_in_p1(i),toAND) > 0;
        else
            SUCp_1 = 0;
        end
        if p2_f
            SUCp_2 = bitand(mux_in_p2(i),toAND) > 0;
        else
            SUCp_2 = 0;
        end
        if p3_f
            SUCp_3 = bitand(mux_in_p3(i),toAND) > 0;
        else
            SUCp_3 = 0;   
        end
        w3 = bitand(128,toAND) > 0;      %1/2
        SUC_weight = SUCp_1 + SUCp_2 + SUCp_3 -...
                    SUCn_1 - SUCn_2 - SUCn_3;
        SUC_sum = SUC_sum + SUC_weight;     %this is the accurate sum or integral of the weight over time
        %test sign FSM as filtering to lower error margin
        if binact > 1
            fraction_cycle = fraction_cycle + 1;
            FSM = FSM + SUC_weight; %
            if fraction_cycle >= fraction_cycle_max
                fraction_cycle = 0;
                FSM = FSM - 1;
            end
            if FSM<0        %set lower upper limit
                FSM=0;
            elseif FSM>FSM_max
                FSM=FSM_max;
            end
            SUC_out = FSM>FSM_half;
        else
            if SUC_weight <0
                SUC_out = 0;
            elseif SUC_weight == 0
                SUC_out = w3;
            else
                SUC_out = 1;
            end
        end
        SUC_record1(1,i,j) = SUC_out;
        SUC_record2(1,i,j) = SUC_weight;
        SUC_record3(1,i,j) = SUC_sum;
    end
    if mod(j,sample_size/10)==0
        fprintf('.');
    end
end
fprintf('\n');
toc
%finishing SC calculation, aggregating sentitive variables out of parfor loop
SUC_record(1,:,:) = SUC_record1(1,:,:); SUC_record1=[]; 
SUC_record(2,:,:) = SUC_record2(1,:,:); SUC_record2=[];
SUC_record(3,:,:) = SUC_record3(1,:,:); SUC_record3=[];
for j = 1:sample_size
    stochastic_sample(1,j)= sum(SUC_record(1,:,j))/runlength;%2*sum(stochastic_val(1,:) == 1)/runlength - 1;
    stochastic_sample(2,j)= binary_sample(2,j);
end
%//////////////////////////////////////////////////////////////////////////
%                   END OF LAYER 1 CONVOLUTION OPERATION
%//////////////////////////////////////////////////////////////////////////
%proceed to layer 2 convolution if any, fixed to mode 5 (SC random
%progamming) timing for SC domain since it shows most promising result.
%Conv2D_L2 is intended to test the pipelining capabiliity of the MUX SUC
%assume using same hardware circuit as conv2D_L1

if conv2D_L2 == 1 %binary computing of Conv2D_L2
    sample_L2 = [];
    if mod(sample_size,9)==0
        sample_size = floor(sample_size/9);
    else
        sample_size = floor(sample_size/9)-1;       %determine new sample size in accordance to the previous sample size
    end
    binary_sample_L2 = zeros(2,sample_size);        %array to store binary output
    fprintf('Conv2D L2 runlength = %d, samples = %d \n', runlength, sample_size);
    %variables borrowed from layer 1 convolution
    %wf,sel_n, sel_p, SUC_record(1,runlength,sample)
    
    %create input data for L2 input using the activated output
    for i = 1:sample_size
       x2 = [binary_sample(1,9*(i-1)+1) binary_sample(1,9*(i-1)+2) binary_sample(1,9*(i-1)+3);...
             binary_sample(1,9*(i-1)+4) binary_sample(1,9*(i-1)+5) binary_sample(1,9*(i-1)+6);...
             binary_sample(1,9*(i-1)+7) binary_sample(1,9*(i-1)+8) binary_sample(1,9*(i-1)+9)];
       sample_L2(:,:,i)=x2;  
    end
    
    %start conv2D_L2 binary computing
    for j = 1:sample_size
        x = sample_L2(:,:,j);  
        % start of neuron iteration in binary domain
        wf = w/255; %convert weight to fraction of 255
        %2nd layer x are all in decimal point, thus the bias has to be in
        %decimal. Conv2D output no longer require /255
        y2_binary = wf.*x;
        y2_binary = sum(y2_binary(:))+bias/255;
        y2_weight = y2_binary;
        
        if binact_L2 == 1
        %clipping output
        if y2_binary < -1 
            y2_binary = -1;
        elseif y2_binary > 1
            y2_binary = 1;
        end
    elseif binact_L2 == 2 %clipping binary output to SC range limit
        if y2_binary < 0 %ReLU
            y2_binary = 0;
        elseif y2_binary > 1 %clipped ReLU
            y2_binary = 1;
        end
        elseif binact_L2 == 3
            %this Conv2D follows the exact Stanh formula
            y2_binary = tanh((2^FSM_bit)*y2_binary/2);
        end
        %shift and scale binary to compare SC in bipolar mode
        if binact_L2 == 1 || binact_L2 == 3    
            y2_binary_s = (y2_binary+1)/2;
        else
            y2_binary_s = y2_binary;
        end
        binary_sample_L2(1,j)= y2_binary;
        binary_sample_L2(2,j)= y2_weight;
        binary_sample_L2(3,j)= y2_binary_s;
    end
    fprintf('Binary computing for Conv_L2 ended, proceeding to SC.\n');
    %Binary part coding ended here
    % end of neuron iteration in binary domain
end
if conv2D_L2 == 1 %stochastic computing of Conv2D_L2
    %////////////////////////////////////////////////////////////////////////////////////////////   
    % start of stochastic domain
    % different from layer 1 convolution, now have to reuse the
    % generated stochastic stream for the computation, but first need
    % to reindex the SUC_record in accourdance to the MUX input
    % sequence
    
    SUC_record = permute(SUC_record, [1,3,2]);  %flip dimension order
    stochastic_sample_L2 = zeros(2,sample_size);    %store SUC adder (mux based) sample data points
    SUC_record_L2 = cast(zeros(3,runlength,sample_size),'int16'); %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record1 = cast(zeros(1,runlength,sample_size),'int16');    %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record2 = cast(zeros(1,runlength,sample_size),'int16'); 
    SUC_record3 = cast(zeros(1,runlength,sample_size),'int16'); 
    SC_x = cast(zeros(3,3,runlength,sample_size),'logical');
    SC_zero = cast(zeros(1,1,runlength),'logical');
    %experimental modification for tanh function
    if binact_L2 == 3
        FSM_max = 2^(FSM_bit)-1+bipolar_FSM_add;  %define FSM max state, must be 2^n-1
    else
        %experimental modification for clipped ReLU function
        FSM_max = 2^(FSM_bit)-1;
    end 
    %define SC zero according to bipolar/unipolar nature
    for i = 1:runlength
        if binact == 2  %unipolar
            SC_zero(i) = 0;
        else            %bipolar
            SC_zero(i) = bitand(128,LFSR8_wbg(i)) > 0; %1/2=0 for bipolar mode
            %SC_zero(i) = 0;
        end
    end
    for j = 1:sample_size %cast the SC output as if it is of matrix input
        SC_x(:,:,:,j) = [SUC_record(1,9*(j-1)+1,:) SUC_record(1,9*(j-1)+2,:) SUC_record(1,9*(j-1)+3,:);...
                        SUC_record(1,9*(j-1)+4,:) SUC_record(1,9*(j-1)+5,:) SUC_record(1,9*(j-1)+6,:);...
                        SUC_record(1,9*(j-1)+7,:) SUC_record(1,9*(j-1)+8,:) SUC_record(1,9*(j-1)+9,:)];
    end
    if layer_reschedule
        fprintf('Rescheduling MUX select for Conv_L2.\n');
        if optimize_sel > 0
            optimize_sel = optimize_sel+1;
        end
        if n1_f
            sel_n1 = gen_sel(mux_n1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n1');
        end
        if n2_f
            sel_n2 = gen_sel(mux_n2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n2');
        end
        if n3_f
            sel_n3 = gen_sel(mux_n3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n3');
        end
        if p1_f
            sel_p1 = gen_sel(mux_p1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p1');
        end
        if p2_f
            sel_p2 = gen_sel(mux_p2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p2');
        end
        if p3_f
            sel_p3 = gen_sel(mux_p3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p3');
        end
    end
    fprintf('Emulating stochastic computing of Conv_L2');
    tic
    %for j = 1:sample_size
    parfor j = 1:sample_size
        %settings for the SC output
        SUC_sum = 0;            %reset SUC sum
        fraction_cycle = 0;     %reset fraction cycle
        FSM = round(FSM_max/2); %reset FSM state
        FSM_half = FSM-1;
        if n1_f
            mux_in_n1 = gen_feed_SC(SC_x(:,:,:,j), mux_n1, sel_n1, SC_zero, runlength);
        end
        if n2_f
            mux_in_n2 = gen_feed_SC(SC_x(:,:,:,j), mux_n2, sel_n2, SC_zero, runlength);
        end
        if n3_f
            mux_in_n3 = gen_feed_SC(SC_x(:,:,:,j), mux_n3, sel_n3, SC_zero, runlength);
        end
        if p1_f
            mux_in_p1 = gen_feed_SC(SC_x(:,:,:,j), mux_p1, sel_p1, SC_zero, runlength);
        end
        if p2_f
            mux_in_p2 = gen_feed_SC(SC_x(:,:,:,j), mux_p2, sel_p2, SC_zero, runlength);
        end
        if p3_f
            mux_in_p3 = gen_feed_SC(SC_x(:,:,:,j), mux_p3, sel_p3, SC_zero, runlength);
        end
        for i = 1:runlength
            %It is worth noting that the SC now is in bipolar mode, except for
            %ReLU activation mode. There is a need to change the SC function
            %in accordance to the previous activation mode     
            if n1_f
                SUCn_1 = mux_in_n1(i);
            else
                SUCn_1 = 0;
            end
            if n2_f
                SUCn_2 = mux_in_n2(i);
            else
                SUCn_2 = 0;
            end
            if n3_f
                SUCn_3 = mux_in_n3(i);
            else
                SUCn_3 = 0;
            end
            if p1_f
                SUCp_1 = mux_in_p1(i);
            else
                SUCp_1 = 0;
            end
            if p2_f
                SUCp_2 = mux_in_p2(i);
            else
                SUCp_2 = 0;
            end
            if p3_f
                SUCp_3 = mux_in_p3(i);
            else
                SUCp_3 = 0;
            end
            w3 = bitand(128,LFSR8_wbg(i)) > 0;
            %has problem assigning bipolar addition
            if binact == 2
                SUC_weight = SUCp_1 + SUCp_2 + SUCp_3 -...
                            SUCn_1 - SUCn_2 - SUCn_3;
            else
                %this part is currently unknown, no study had been found
                %that uses SUC MUX to process bipolar stochastic stream.
                %This code is designed to fit wide range of SUC weights by
                %disabling specific value for non-existing streams.
                SUC_weight = p1_f*(2*SUCp_1-1) + p2_f*(2*SUCp_2-1) + p3_f*(2*SUCp_3-1) - ...
                            n1_f*(2*SUCn_1-1) - n2_f*(2*SUCn_2-1) - n3_f*(2*SUCn_3-1);
            end
            SUC_sum = SUC_sum + SUC_weight;
            %test sign FSM as filtering to lower error margin
            if binact_L2 > 1  
                fraction_cycle = fraction_cycle + 1;
                FSM = FSM + SUC_weight; %for usual SUC
                %FSM = FSM + (2*SUC_weight-1); %for modded SUC after MUX subtractor
                if fraction_cycle >= fraction_cycle_max
                    fraction_cycle = 0;
                    FSM = FSM - 1;
                end
                if FSM<0        %set lower upper limit
                    FSM=0;
                elseif FSM>FSM_max
                    FSM=FSM_max;
                end
                SUC_out = FSM>FSM_half;
            else
                if SUC_weight <0
                    SUC_out = 0;
                elseif SUC_weight == 0
                    SUC_out = w3;
                else
                    SUC_out = 1;
                end
            end
            SUC_record1(1,i,j) = SUC_out;
            SUC_record2(1,i,j) = SUC_weight;
            SUC_record3(1,i,j) = SUC_sum;
        end
        if mod(j,sample_size/10)<1
            fprintf('.');
        end
    end
    fprintf('\n');
    toc
    %finishing SC calculation, aggregating sentitive variables out of parfor loop
    SUC_record_L2(1,:,:) = SUC_record1(1,:,:); SUC_record1=[]; 
    SUC_record_L2(2,:,:) = SUC_record2(1,:,:); SUC_record2=[];
    SUC_record_L2(3,:,:) = SUC_record3(1,:,:); SUC_record3=[];
    for j = 1:sample_size
        stochastic_sample_L2(1,j)= sum(SUC_record_L2(1,:,j))/runlength;%2*sum(stochastic_val(1,:) == 1)/runlength - 1;
        stochastic_sample_L2(2,j)= binary_sample_L2(2,j);
    end
end
%//////////////////////////////////////////////////////////////////////////
%                   END OF LAYER 2 CONVOLUTION OPERATION
%//////////////////////////////////////////////////////////////////////////
if conv2D_L3 == 1 %binary computing of Conv2D_L3
    sample_L3 = [];
    if mod(sample_size,9)==0
        sample_size = floor(sample_size/9);
    else
        sample_size = floor(sample_size/9)-1;       %determine new sample size in accordance to the previous sample size
    end
    binary_sample_L3 = zeros(2,sample_size);        %array to store binary output
    fprintf('Conv2D L3 runlength = %d, samples = %d \n', runlength, sample_size);
    %variables borrowed from layer 1 convolution
    %wf,sel_n, sel_p, SUC_record(1,runlength,sample)
    
    %create input data for L2 input using the activated output
    for i = 1:sample_size
       x3 = [binary_sample_L2(1,9*(i-1)+1) binary_sample_L2(1,9*(i-1)+2) binary_sample_L2(1,9*(i-1)+3);...
             binary_sample_L2(1,9*(i-1)+4) binary_sample_L2(1,9*(i-1)+5) binary_sample_L2(1,9*(i-1)+6);...
             binary_sample_L2(1,9*(i-1)+7) binary_sample_L2(1,9*(i-1)+8) binary_sample_L2(1,9*(i-1)+9)];
       sample_L3(:,:,i)=x3;  
    end
    
    %start conv2D_L2 binary computing
    for j = 1:sample_size
        x = sample_L3(:,:,j);  
        % start of neuron iteration in binary domain
        wf = w/255; %convert weight to fraction of 255
        %2nd layer x are all in decimal point, thus the bias has to be in
        %decimal. Conv2D output no longer require /255
        y3_binary = wf.*x;
        y3_binary = sum(y3_binary(:))+bias/255;
        y3_weight = y3_binary;
        
        if binact_L3 == 1
        %clipping output
        if y3_binary < -1 
            y3_binary = -1;
        elseif y3_binary > 1
            y3_binary = 1;
        end
    elseif binact_L3 == 2 %clipping binary output to SC range limit
        if y3_binary < 0 %ReLU
            y3_binary = 0;
        elseif y3_binary > 1 %clipped ReLU
            y3_binary = 1;
        end
        elseif binact_L3 == 3
            %this Conv2D follows the exact Stanh formula
            y3_binary = tanh((2^FSM_bit)*y3_binary/2);
        end
        %shift and scale binary to compare SC in bipolar mode
        if binact_L3 == 1 || binact_L3 == 3    
            y3_binary_s = (y3_binary+1)/2;
        else
            y3_binary_s = y3_binary;
        end
        binary_sample_L3(1,j)= y3_binary;
        binary_sample_L3(2,j)= y3_weight;
        binary_sample_L3(3,j)= y3_binary_s;
    end
    fprintf('Binary computing for Conv_L3 ended, proceeding to SC.\n');
    %Binary part coding ended here
    % end of neuron iteration in binary domain
end
if conv2D_L3 == 1 %stochastic computing of Conv2D_L2
    %////////////////////////////////////////////////////////////////////////////////////////////   
    % start of stochastic domain
    % different from layer 1 convolution, now have to reuse the
    % generated stochastic stream for the computation, but first need
    % to reindex the SUC_record in accourdance to the MUX input
    % sequence
    
    SUC_record_L2 = permute(SUC_record_L2, [1,3,2]);  %flip dimension order
    stochastic_sample_L3 = zeros(2,sample_size);    %store SUC adder (mux based) sample data points
    SUC_record_L3 = cast(zeros(3,runlength,sample_size),'int16'); %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record1 = cast(zeros(1,runlength,sample_size),'int16');    %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record2 = cast(zeros(1,runlength,sample_size),'int16'); 
    SUC_record3 = cast(zeros(1,runlength,sample_size),'int16'); 
    SC_x = cast(zeros(3,3,runlength,sample_size),'logical');
    SC_zero = cast(zeros(1,1,runlength),'logical');
    %experimental modification for tanh function
    if binact_L3 == 3
        FSM_max = 2^(FSM_bit)-1+bipolar_FSM_add;  %define FSM max state, must be 2^n-1
    else
        %experimental modification for clipped ReLU function
        FSM_max = 2^(FSM_bit)-1;
    end 
    %define SC zero according to bipolar/unipolar nature
    for i = 1:runlength
        if binact_L2 == 2  %unipolar
            SC_zero(i) = 0;
        else            %bipolar
            SC_zero(i) = bitand(128,LFSR8_wbg(i)) > 0; %1/2=0 for bipolar mode
            %SC_zero(i) = 0;
        end
    end
    for j = 1:sample_size %cast the SC output as if it is of matrix input
        SC_x(:,:,:,j) = [SUC_record_L2(1,9*(j-1)+1,:) SUC_record_L2(1,9*(j-1)+2,:) SUC_record_L2(1,9*(j-1)+3,:);...
                SUC_record_L2(1,9*(j-1)+4,:) SUC_record_L2(1,9*(j-1)+5,:) SUC_record_L2(1,9*(j-1)+6,:);...
                SUC_record_L2(1,9*(j-1)+7,:) SUC_record_L2(1,9*(j-1)+8,:) SUC_record_L2(1,9*(j-1)+9,:)];
    end
    if layer_reschedule
        fprintf('Rescheduling MUX select for Conv_L3.\n');
        if optimize_sel > 0
            optimize_sel = optimize_sel+1;
        end
        if n1_f
            sel_n1 = gen_sel(mux_n1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n1');
        end
        if n2_f
            sel_n2 = gen_sel(mux_n2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n2');
        end
        if n3_f
            sel_n3 = gen_sel(mux_n3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n3');
        end
        if p1_f
            sel_p1 = gen_sel(mux_p1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p1');
        end
        if p2_f
            sel_p2 = gen_sel(mux_p2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p2');
        end
        if p3_f
            sel_p3 = gen_sel(mux_p3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p3');
        end
    end
    fprintf('Emulating stochastic computing of Conv_L3');
    tic
    %for j = 1:sample_size
    parfor j = 1:sample_size
        %settings for the SC output
        SUC_sum = 0;            %reset SUC sum
        fraction_cycle = 0;     %reset fraction cycle
        FSM = round(FSM_max/2); %reset FSM state
        FSM_half = FSM-1;
        
        if n1_f
            mux_in_n1 = gen_feed_SC(SC_x(:,:,:,j), mux_n1, sel_n1, SC_zero, runlength);
        end
        if n2_f
            mux_in_n2 = gen_feed_SC(SC_x(:,:,:,j), mux_n2, sel_n2, SC_zero, runlength);
        end
        if n3_f
            mux_in_n3 = gen_feed_SC(SC_x(:,:,:,j), mux_n3, sel_n3, SC_zero, runlength);
        end
        if p1_f
            mux_in_p1 = gen_feed_SC(SC_x(:,:,:,j), mux_p1, sel_p1, SC_zero, runlength);
        end
        if p2_f
            mux_in_p2 = gen_feed_SC(SC_x(:,:,:,j), mux_p2, sel_p2, SC_zero, runlength);
        end
        if p3_f
            mux_in_p3 = gen_feed_SC(SC_x(:,:,:,j), mux_p3, sel_p3, SC_zero, runlength);
        end
        for i = 1:runlength
            %It is worth noting that the SC now is in bipolar mode, except for
            %ReLU activation mode. There is a need to change the SC function
            %in accordance to the previous activation mode     
            if n1_f
                SUCn_1 = mux_in_n1(i);
            else
                SUCn_1 = 0;
            end
            if n2_f
                SUCn_2 = mux_in_n2(i);
            else
                SUCn_2 = 0;
            end
            if n3_f
                SUCn_3 = mux_in_n3(i);
            else
                SUCn_3 = 0;
            end
            if p1_f
                SUCp_1 = mux_in_p1(i);
            else
                SUCp_1 = 0;
            end
            if p2_f
                SUCp_2 = mux_in_p2(i);
            else
                SUCp_2 = 0;
            end
            if p3_f
                SUCp_3 = mux_in_p3(i);
            else
                SUCp_3 = 0;
            end
            w3 = bitand(128,LFSR8_wbg(i)) > 0;
            %has problem assigning bipolar addition
            if binact_L2 == 2
                SUC_weight = SUCp_1 + SUCp_2 + SUCp_3 -...
                            SUCn_1 - SUCn_2 - SUCn_3;
            else
                %this part is currently unknown, no study had been found
                %that uses SUC MUX to process bipolar stochastic stream.
                %This code is designed to fit wide range of SUC weights by
                %disabling specific value for non-existing streams.
                SUC_weight = p1_f*(2*SUCp_1-1) + p2_f*(2*SUCp_2-1) + p3_f*(2*SUCp_3-1) - ...
                            n1_f*(2*SUCn_1-1) - n2_f*(2*SUCn_2-1) - n3_f*(2*SUCn_3-1);
            end
            SUC_sum = SUC_sum + SUC_weight;
            %test sign FSM as filtering to lower error margin
            if binact_L3 > 1  
                fraction_cycle = fraction_cycle + 1;
                FSM = FSM + SUC_weight; %for usual SUC
                %FSM = FSM + (2*SUC_weight-1); %for modded SUC after MUX subtractor
                if fraction_cycle >= fraction_cycle_max
                    fraction_cycle = 0;
                    FSM = FSM - 1;
                end
                if FSM<0        %set lower upper limit
                    FSM=0;
                elseif FSM>FSM_max
                    FSM=FSM_max;
                end
                SUC_out = FSM>FSM_half;
            else
                if SUC_weight <0
                    SUC_out = 0;
                elseif SUC_weight == 0
                    SUC_out = w3;
                else
                    SUC_out = 1;
                end
            end
            SUC_record1(1,i,j) = SUC_out;
            SUC_record2(1,i,j) = SUC_weight;
            SUC_record3(1,i,j) = SUC_sum;
        end
        if mod(j,sample_size/10)<1
            fprintf('.');
        end
    end
    fprintf('\n');
    toc
    %finishing SC calculation, aggregating sentitive variables out of parfor loop
    SUC_record_L3(1,:,:) = SUC_record1(1,:,:); SUC_record1=[]; 
    SUC_record_L3(2,:,:) = SUC_record2(1,:,:); SUC_record2=[];
    SUC_record_L3(3,:,:) = SUC_record3(1,:,:); SUC_record3=[];
    for j = 1:sample_size
        stochastic_sample_L3(1,j)= sum(SUC_record_L3(1,:,j))/runlength;%2*sum(stochastic_val(1,:) == 1)/runlength - 1;
        stochastic_sample_L3(2,j)= binary_sample_L3(2,j);
    end
end
%//////////////////////////////////////////////////////////////////////////
%                   END OF LAYER 3 CONVOLUTION OPERATION
%//////////////////////////////////////////////////////////////////////////
if conv2D_L4 == 1 %binary computing of Conv2D_L4
    sample_L4 = [];
    if mod(sample_size,9)==0
        sample_size = floor(sample_size/9);
    else
        sample_size = floor(sample_size/9)-1;       %determine new sample size in accordance to the previous sample size
    end
    binary_sample_L4 = zeros(2,sample_size);        %array to store binary output
    fprintf('Conv2D L4 runlength = %d, samples = %d \n', runlength, sample_size);
    %variables borrowed from layer 1 convolution
    %wf,sel_n, sel_p, SUC_record(1,runlength,sample)
    
    %create input data for L2 input using the activated output
    for i = 1:sample_size
       x4 = [binary_sample_L3(1,9*(i-1)+1) binary_sample_L3(1,9*(i-1)+2) binary_sample_L3(1,9*(i-1)+3);...
             binary_sample_L3(1,9*(i-1)+4) binary_sample_L3(1,9*(i-1)+5) binary_sample_L3(1,9*(i-1)+6);...
             binary_sample_L3(1,9*(i-1)+7) binary_sample_L3(1,9*(i-1)+8) binary_sample_L3(1,9*(i-1)+9)];
       sample_L4(:,:,i)=x4;  
    end
    
    %start conv2D_L2 binary computing
    for j = 1:sample_size
        x = sample_L4(:,:,j);  
        % start of neuron iteration in binary domain
        wf = w/255; %convert weight to fraction of 255
        %2nd layer x are all in decimal point, thus the bias has to be in
        %decimal. Conv2D output no longer require /255
        y4_binary = wf.*x;
        y4_binary = sum(y4_binary(:))+bias/255;
        y4_weight = y4_binary;
        
        if binact_L4 == 1
        %clipping output
        if y4_binary < -1 
            y4_binary = -1;
        elseif y4_binary > 1
            y4_binary = 1;
        end
    elseif binact_L4 == 2 %clipping binary output to SC range limit
        if y4_binary < 0 %ReLU
            y4_binary = 0;
        elseif y4_binary > 1 %clipped ReLU
            y4_binary = 1;
        end
        elseif binact_L4 == 3
            %this Conv2D follows the exact Stanh formula
            y4_binary = tanh((2^FSM_bit)*y4_binary/2);
        end
        %shift and scale binary to compare SC in bipolar mode
        if binact_L4 == 1 || binact_L4 == 3    
            y4_binary_s = (y4_binary+1)/2;
        else
            y4_binary_s = y4_binary;
        end
        binary_sample_L4(1,j)= y4_binary;
        binary_sample_L4(2,j)= y4_weight;
        binary_sample_L4(3,j)= y4_binary_s;
    end
    fprintf('Binary computing for Conv_L4 ended, proceeding to SC.\n');
    %Binary part coding ended here
    % end of neuron iteration in binary domain
end
if conv2D_L4 == 1 %stochastic computing of Conv2D_L4
    %////////////////////////////////////////////////////////////////////////////////////////////   
    % start of stochastic domain
    % different from layer 1 convolution, now have to reuse the
    % generated stochastic stream for the computation, but first need
    % to reindex the SUC_record in accourdance to the MUX input
    % sequence
    SUC_record_L3 = permute(SUC_record_L3, [1,3,2]);  %flip dimension order
    stochastic_sample_L4 = zeros(2,sample_size);    %store SUC adder (mux based) sample data points
    SUC_record_L4 = cast(zeros(3,runlength,sample_size),'int16'); %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record1 = cast(zeros(1,runlength,sample_size),'int16');    %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record2 = cast(zeros(1,runlength,sample_size),'int16'); 
    SUC_record3 = cast(zeros(1,runlength,sample_size),'int16'); 
    SC_x = cast(zeros(3,3,runlength,sample_size),'logical');
    SC_zero = cast(zeros(1,1,runlength),'logical');
    %experimental modification for tanh function
    if binact_L4 == 3
        FSM_max = 2^(FSM_bit)-1+bipolar_FSM_add;  %define FSM max state, must be 2^n-1
    else
        %experimental modification for clipped ReLU function
        FSM_max = 2^(FSM_bit)-1;
    end 
    %define SC zero according to bipolar/unipolar nature
    for i = 1:runlength
        if binact_L3 == 2  %unipolar
            SC_zero(i) = 0;
        else            %bipolar
            SC_zero(i) = bitand(128,LFSR8_wbg(i)) > 0; %1/2=0 for bipolar mode
            %SC_zero(i) = 0;
        end
    end
    for j = 1:sample_size
        SC_x(:,:,:,j) = [SUC_record_L3(1,9*(j-1)+1,:) SUC_record_L3(1,9*(j-1)+2,:) SUC_record_L3(1,9*(j-1)+3,:);...
                SUC_record_L3(1,9*(j-1)+4,:) SUC_record_L3(1,9*(j-1)+5,:) SUC_record_L3(1,9*(j-1)+6,:);...
                SUC_record_L3(1,9*(j-1)+7,:) SUC_record_L3(1,9*(j-1)+8,:) SUC_record_L3(1,9*(j-1)+9,:)];
        %cast the SC output as if it is of matrix input
    end
    if layer_reschedule
        fprintf('Rescheduling MUX select for Conv_L4.\n');
        if optimize_sel > 0
            optimize_sel = optimize_sel+1;
        end
        if n1_f
            sel_n1 = gen_sel(mux_n1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n1');
        end
        if n2_f
            sel_n2 = gen_sel(mux_n2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n2');
        end
        if n3_f
            sel_n3 = gen_sel(mux_n3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n3');
        end
        if p1_f
            sel_p1 = gen_sel(mux_p1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p1');
        end
        if p2_f
            sel_p2 = gen_sel(mux_p2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p2');
        end
        if p3_f
            sel_p3 = gen_sel(mux_p3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p3');
        end
    end
    fprintf('Emulating stochastic computing of Conv_L4');
    tic
    %for j = 1:sample_size
    parfor j = 1:sample_size
        %settings for the SC output
        SUC_sum = 0;            %reset SUC sum
        fraction_cycle = 0;     %reset fraction cycle
        FSM = round(FSM_max/2); %reset FSM state
        FSM_half = FSM-1;
        if n1_f
            mux_in_n1 = gen_feed_SC(SC_x(:,:,:,j), mux_n1, sel_n1, SC_zero, runlength);
        end
        if n2_f
            mux_in_n2 = gen_feed_SC(SC_x(:,:,:,j), mux_n2, sel_n2, SC_zero, runlength);
        end
        if n3_f
            mux_in_n3 = gen_feed_SC(SC_x(:,:,:,j), mux_n3, sel_n3, SC_zero, runlength);
        end
        if p1_f
            mux_in_p1 = gen_feed_SC(SC_x(:,:,:,j), mux_p1, sel_p1, SC_zero, runlength);
        end
        if p2_f
            mux_in_p2 = gen_feed_SC(SC_x(:,:,:,j), mux_p2, sel_p2, SC_zero, runlength);
        end
        if p3_f
            mux_in_p3 = gen_feed_SC(SC_x(:,:,:,j), mux_p3, sel_p3, SC_zero, runlength);
        end
        for i = 1:runlength
            %It is worth noting that the SC now is in bipolar mode, except for
            %ReLU activation mode. There is a need to change the SC function
            %in accordance to the previous activation mode     
            if n1_f
                SUCn_1 = mux_in_n1(i);
            else
                SUCn_1 = 0;
            end
            if n2_f
                SUCn_2 = mux_in_n2(i);
            else
                SUCn_2 = 0;
            end
            if n3_f
                SUCn_3 = mux_in_n3(i);
            else
                SUCn_3 = 0;
            end
            if p1_f
                SUCp_1 = mux_in_p1(i);
            else
                SUCp_1 = 0;
            end
            if p2_f
                SUCp_2 = mux_in_p2(i);
            else
                SUCp_2 = 0;
            end
            if p3_f
                SUCp_3 = mux_in_p3(i);
            else
                SUCp_3 = 0;
            end
            w3 = bitand(128,LFSR8_wbg(i)) > 0;
            %has problem assigning bipolar addition
            if binact_L3 == 2
                SUC_weight = SUCp_1 + SUCp_2 + SUCp_3 -...
                            SUCn_1 - SUCn_2 - SUCn_3;
            else
                %this part is currently unknown, no study had been found
                %that uses SUC MUX to process bipolar stochastic stream.
                %This code is designed to fit wide range of SUC weights by
                %disabling specific value for non-existing streams.
                SUC_weight = p1_f*(2*SUCp_1-1) + p2_f*(2*SUCp_2-1) + p3_f*(2*SUCp_3-1) - ...
                            n1_f*(2*SUCn_1-1) - n2_f*(2*SUCn_2-1) - n3_f*(2*SUCn_3-1);
            end
            SUC_sum = SUC_sum + SUC_weight;
            %test sign FSM as filtering to lower error margin
            if binact_L4 > 1  
                fraction_cycle = fraction_cycle + 1;
                FSM = FSM + SUC_weight; %for usual SUC
                %FSM = FSM + (2*SUC_weight-1); %for modded SUC after MUX subtractor
                if fraction_cycle >= fraction_cycle_max
                    fraction_cycle = 0;
                    FSM = FSM - 1;
                end
                if FSM<0        %set lower upper limit
                    FSM=0;
                elseif FSM>FSM_max
                    FSM=FSM_max;
                end
                SUC_out = FSM>FSM_half;
            else
                if SUC_weight <0
                    SUC_out = 0;
                elseif SUC_weight == 0
                    SUC_out = w3;
                else
                    SUC_out = 1;
                end
            end
            SUC_record1(1,i,j) = SUC_out;
            SUC_record2(1,i,j) = SUC_weight;
            SUC_record3(1,i,j) = SUC_sum;
        end
        if mod(j,sample_size/10)<1
            fprintf('.');
        end
    end
    fprintf('\n');
    toc
    %finishing SC calculation, aggregating sentitive variables out of parfor loop
    SUC_record_L4(1,:,:) = SUC_record1(1,:,:); SUC_record1=[]; 
    SUC_record_L4(2,:,:) = SUC_record2(1,:,:); SUC_record2=[];
    SUC_record_L4(3,:,:) = SUC_record3(1,:,:); SUC_record3=[];
    for j = 1:sample_size
        stochastic_sample_L4(1,j)= sum(SUC_record_L4(1,:,j))/runlength;%2*sum(stochastic_val(1,:) == 1)/runlength - 1;
        stochastic_sample_L4(2,j)= binary_sample_L4(2,j);
    end
end

%//////////////////////////////////////////////////////////////////////////
%                   END OF LAYER 4 CONVOLUTION OPERATION
%//////////////////////////////////////////////////////////////////////////
if conv2D_L5 == 1 %binary computing of Conv2D_L5
    sample_L5 = [];
    if mod(sample_size,9)==0
        sample_size = floor(sample_size/9);
    else
        sample_size = floor(sample_size/9)-1;       %determine new sample size in accordance to the previous sample size
    end
    binary_sample_L5 = zeros(2,sample_size);        %array to store binary output
    fprintf('Conv2D L5 runlength = %d, samples = %d \n', runlength, sample_size);
    %variables borrowed from layer 1 convolution
    %wf,sel_n, sel_p, SUC_record(1,runlength,sample)
    
    %create input data for L2 input using the activated output
    for i = 1:sample_size
       x5 = [binary_sample_L4(1,9*(i-1)+1) binary_sample_L4(1,9*(i-1)+2) binary_sample_L4(1,9*(i-1)+3);...
             binary_sample_L4(1,9*(i-1)+4) binary_sample_L4(1,9*(i-1)+5) binary_sample_L4(1,9*(i-1)+6);...
             binary_sample_L4(1,9*(i-1)+7) binary_sample_L4(1,9*(i-1)+8) binary_sample_L4(1,9*(i-1)+9)];
       sample_L5(:,:,i)=x5;  
    end
    
    %start conv2D_L2 binary computing
    for j = 1:sample_size
        x = sample_L5(:,:,j);  
        % start of neuron iteration in binary domain
        wf = w/255; %convert weight to fraction of 255
        %2nd layer x are all in decimal point, thus the bias has to be in
        %decimal. Conv2D output no longer require /255
        y5_binary = wf.*x;
        y5_binary = sum(y5_binary(:))+bias/255;
        y5_weight = y5_binary;
        
        if binact_L5 == 1
        %clipping output
        if y5_binary < -1 
            y5_binary = -1;
        elseif y5_binary > 1
            y5_binary = 1;
        end
    elseif binact_L5 == 2 %clipping binary output to SC range limit
        if y5_binary < 0 %ReLU
            y5_binary = 0;
        elseif y5_binary > 1 %clipped ReLU
            y5_binary = 1;
        end
        elseif binact_L5 == 3
            %this Conv2D follows the exact Stanh formula
            y5_binary = tanh((2^FSM_bit)*y5_binary/2);
        end
        %shift and scale binary to compare SC in bipolar mode
        if binact_L5 == 1 || binact_L5 == 3    
            y5_binary_s = (y5_binary+1)/2;
        else
            y5_binary_s = y5_binary;
        end
        binary_sample_L5(1,j)= y5_binary;
        binary_sample_L5(2,j)= y5_weight;
        binary_sample_L5(3,j)= y5_binary_s;
    end
    fprintf('Binary computing for Conv_L5 ended, proceeding to SC.\n');
    %Binary part coding ended here
    % end of neuron iteration in binary domain
end
if conv2D_L5 == 1 %stochastic computing of Conv2D_L5
    %////////////////////////////////////////////////////////////////////////////////////////////   
    % start of stochastic domain
    % different from layer 1 convolution, now have to reuse the
    % generated stochastic stream for the computation, but first need
    % to reindex the SUC_record in accourdance to the MUX input
    % sequence
    
    SUC_record_L4 = permute(SUC_record_L4, [1,3,2]);  %flip dimension order
    stochastic_sample_L5 = zeros(2,sample_size);    %store SUC adder (mux based) sample data points
    SUC_record_L5 = cast(zeros(3,runlength,sample_size),'int16'); %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record1 = cast(zeros(1,runlength,sample_size),'int16');    %store SUC adder output (1=suc out, 2=suc weight, 3=suc integral
    SUC_record2 = cast(zeros(1,runlength,sample_size),'int16'); 
    SUC_record3 = cast(zeros(1,runlength,sample_size),'int16'); 
    SC_x = cast(zeros(3,3,runlength,sample_size),'logical');
    SC_zero = cast(zeros(1,1,runlength),'logical');
    %experimental modification for tanh function
    if binact_L5 == 3
        FSM_max = 2^(FSM_bit)-1+bipolar_FSM_add;  %define FSM max state, must be 2^n-1
    else
        %experimental modification for clipped ReLU function
        FSM_max = 2^(FSM_bit)-1;
    end 
    %define SC zero according to bipolar/unipolar nature
    for i = 1:runlength
        if binact_L4 == 2  %unipolar
            SC_zero(i) = 0;
        else            %bipolar
            SC_zero(i) = bitand(128,LFSR8_wbg(i)) > 0; %1/2=0 for bipolar mode
            %SC_zero(i) = 0;
        end
    end
    for j = 1:sample_size
        SC_x(:,:,:,j) = [SUC_record_L4(1,9*(j-1)+1,:) SUC_record_L4(1,9*(j-1)+2,:) SUC_record_L4(1,9*(j-1)+3,:);...
                SUC_record_L4(1,9*(j-1)+4,:) SUC_record_L4(1,9*(j-1)+5,:) SUC_record_L4(1,9*(j-1)+6,:);...
                SUC_record_L4(1,9*(j-1)+7,:) SUC_record_L4(1,9*(j-1)+8,:) SUC_record_L4(1,9*(j-1)+9,:)];
        %cast the SC output as if it is of matrix input
    end
    if layer_reschedule
        fprintf('Rescheduling MUX select for Conv_L5.\n');
        if optimize_sel > 0
            optimize_sel = optimize_sel+1;
        end
        if n1_f
            sel_n1 = gen_sel(mux_n1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n1');
        end
        if n2_f
            sel_n2 = gen_sel(mux_n2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n2');
        end
        if n3_f
            sel_n3 = gen_sel(mux_n3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_n3');
        end
        if p1_f
            sel_p1 = gen_sel(mux_p1(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p1');
        end
        if p2_f
            sel_p2 = gen_sel(mux_p2(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p2');
        end
        if p3_f
            sel_p3 = gen_sel(mux_p3(3,:), LFSR8, runlength, mode, optimize_sel, 'MUX_p3');
        end
    end
    fprintf('Emulating stochastic computing of Conv_L5');
    tic
    for j = 1:sample_size
    %parfor j = 1:sample_size
        %settings for the SC output
        SUC_sum = 0;            %reset SUC sum
        fraction_cycle = 0;     %reset fraction cycle
        FSM = round(FSM_max/2); %reset FSM state
        FSM_half = FSM-1;      
        if n1_f
            mux_in_n1 = gen_feed_SC(SC_x(:,:,:,j), mux_n1, sel_n1, SC_zero, runlength);
        end
        if n2_f
            mux_in_n2 = gen_feed_SC(SC_x(:,:,:,j), mux_n2, sel_n2, SC_zero, runlength);
        end
        if n3_f
            mux_in_n3 = gen_feed_SC(SC_x(:,:,:,j), mux_n3, sel_n3, SC_zero, runlength);
        end
        if p1_f
            mux_in_p1 = gen_feed_SC(SC_x(:,:,:,j), mux_p1, sel_p1, SC_zero, runlength);
        end
        if p2_f
            mux_in_p2 = gen_feed_SC(SC_x(:,:,:,j), mux_p2, sel_p2, SC_zero, runlength);
        end
        if p3_f
            mux_in_p3 = gen_feed_SC(SC_x(:,:,:,j), mux_p3, sel_p3, SC_zero, runlength);
        end
        for i = 1:runlength
            %It is worth noting that the SC now is in bipolar mode, except for
            %ReLU activation mode. There is a need to change the SC function
            %in accordance to the previous activation mode     
            if n1_f
                SUCn_1 = mux_in_n1(i);
            else
                SUCn_1 = 0;
            end
            if n2_f
                SUCn_2 = mux_in_n2(i);
            else
                SUCn_2 = 0;
            end
            if n3_f
                SUCn_3 = mux_in_n3(i);
            else
                SUCn_3 = 0;
            end
            if p1_f
                SUCp_1 = mux_in_p1(i);
            else
                SUCp_1 = 0;
            end
            if p2_f
                SUCp_2 = mux_in_p2(i);
            else
                SUCp_2 = 0;
            end
            if p3_f
                SUCp_3 = mux_in_p3(i);
            else
                SUCp_3 = 0;
            end
            w3 = bitand(128,LFSR8_wbg(i)) > 0;
            %has problem assigning bipolar addition
            if binact_L4 == 2
                SUC_weight = SUCp_1 + SUCp_2 + SUCp_3 -...
                            SUCn_1 - SUCn_2 - SUCn_3;
            else
                %this part is currently unknown, no study had been found
                %that uses SUC MUX to process bipolar stochastic stream.
                %This code is designed to fit wide range of SUC weights by
                %disabling specific value for non-existing streams.
                SUC_weight = p1_f*(2*SUCp_1-1) + p2_f*(2*SUCp_2-1) + p3_f*(2*SUCp_3-1) - ...
                            n1_f*(2*SUCn_1-1) - n2_f*(2*SUCn_2-1) - n3_f*(2*SUCn_3-1);
            end
            SUC_sum = SUC_sum + SUC_weight;
            %test sign FSM as filtering to lower error margin
            if binact_L5 > 1  
                fraction_cycle = fraction_cycle + 1;
                FSM = FSM + SUC_weight; %for usual SUC
                %FSM = FSM + (2*SUC_weight-1); %for modded SUC after MUX subtractor
                if fraction_cycle >= fraction_cycle_max
                    fraction_cycle = 0;
                    FSM = FSM - 1;
                end
                if FSM<0        %set lower upper limit
                    FSM=0;
                elseif FSM>FSM_max
                    FSM=FSM_max;
                end
                SUC_out = FSM>FSM_half;
            else
                if SUC_weight <0
                    SUC_out = 0;
                elseif SUC_weight == 0
                    SUC_out = w3;
                else
                    SUC_out = 1;
                end
            end
            SUC_record1(1,i,j) = SUC_out;
            SUC_record2(1,i,j) = SUC_weight;
            SUC_record3(1,i,j) = SUC_sum;
        end
        if mod(j,sample_size/10)<1
            fprintf('.');
        end
    end
    fprintf('\n');
    toc
    %finishing SC calculation, aggregating sentitive variables out of parfor loop
    SUC_record_L5(1,:,:) = SUC_record1(1,:,:); SUC_record1=[]; 
    SUC_record_L5(2,:,:) = SUC_record2(1,:,:); SUC_record2=[];
    SUC_record_L5(3,:,:) = SUC_record3(1,:,:); SUC_record3=[];
    for j = 1:sample_size
        stochastic_sample_L5(1,j)= sum(SUC_record_L5(1,:,j))/runlength;%2*sum(stochastic_val(1,:) == 1)/runlength - 1;
        stochastic_sample_L5(2,j)= binary_sample_L5(2,j);
    end
end
%//////////////////////////////////////////////////////////////////////////
%                   END OF LAYER 5 CONVOLUTION OPERATION
%//////////////////////////////////////////////////////////////////////////
figure(1)
if conv2D_L5 == 1
    tiledlayout(5,2)
elseif conv2D_L4 == 1
    tiledlayout(4,2)
elseif conv2D_L3 == 1
    tiledlayout(3,2)
elseif conv2D_L2 == 1
    tiledlayout(2,2)
else
    tiledlayout(1,2)
end
nexttile
scatter(stochastic_sample(2,:),stochastic_sample(1,:),'.')  %plot SC output
hold on
scatter(binary_sample(2,:),binary_sample(3,:),'x') %plot scaled binary output
hold off
title('Conv2D L1 SUC MUX vs binary scatter')
nexttile
histogram(stochastic_sample(1,:)-binary_sample(3,:)) %plot error distribution of subject vs binary output
title('Conv2D L1 SUC MUX error distribution')

if conv2D_L2 == 1
    nexttile
    scatter(stochastic_sample_L2(2,:),stochastic_sample_L2(1,:),'.')  %plot SC output
    hold on
    scatter(binary_sample_L2(2,:),binary_sample_L2(3,:),'x') %plot binary output
    hold off
    title('Conv2D L2 SUC MUX vs binary scatter')
    nexttile
    histogram(stochastic_sample_L2(1,:)-binary_sample_L2(3,:)) %plot error distribution of subject vs binary output
    title('Conv2D L2 SUC MUX error distribution')
end
if conv2D_L3 == 1
    nexttile
    scatter(stochastic_sample_L3(2,:),stochastic_sample_L3(1,:),'.')  %plot SC output
    hold on
    scatter(binary_sample_L3(2,:),binary_sample_L3(3,:),'x') %plot binary output
    hold off
    title('Conv2D L3 SUC MUX vs binary scatter')
    nexttile
    histogram(stochastic_sample_L3(1,:)-binary_sample_L3(3,:)) %plot error distribution of subject vs binary output
    title('Conv2D L3 SUC MUX error distribution')
end
if conv2D_L4 == 1
    nexttile
    scatter(stochastic_sample_L4(2,:),stochastic_sample_L4(1,:),'.')  %plot SC output
    hold on
    scatter(binary_sample_L4(2,:),binary_sample_L4(3,:),'x') %plot binary output
    hold off
    title('Conv2D L4 SUC MUX vs binary scatter')
    nexttile
    histogram(stochastic_sample_L4(1,:)-binary_sample_L4(3,:)) %plot error distribution of subject vs binary output
    title('Conv2D L4 SUC MUX error distribution')
end
if conv2D_L5 == 1
    nexttile
    scatter(stochastic_sample_L5(2,:),stochastic_sample_L5(1,:),'.')  %plot SC output
    hold on
    scatter(binary_sample_L5(2,:),binary_sample_L5(3,:),'x') %plot binary output
    hold off
    title('Conv2D L5 SUC MUX vs binary scatter')
    nexttile
    histogram(stochastic_sample_L5(1,:)-binary_sample_L5(3,:)) %plot error distribution of subject vs binary output
    title('Conv2D L5 SUC MUX error distribution')
end
%plot the intermediate response of SUC adder of first layer
if plot_response == 1
    plot_SUC_response(SUC_record,'L1',10);
    if conv2D_L2 == 1
        plot_SUC_response(SUC_record_L2,'L2',11);
    end
    if conv2D_L3 == 1
        plot_SUC_response(SUC_record_L3,'L3',12);
    end
    if conv2D_L4 == 1
        plot_SUC_response(SUC_record_L4,'L4',13);
    end
    if conv2D_L5 == 1
        plot_SUC_response(SUC_record_L5,'L5',14);
    end
end

fprintf('L1 error standard deviation = %f, L1 mean error = %f \n', ...
    std(stochastic_sample(1,:)-binary_sample(3,:)), mean(stochastic_sample(1,:)-binary_sample(3,:)));
if conv2D_L2 == 1
    fprintf('L2 error standard deviation = %f, L2 mean error = %f \n', ...
    std(stochastic_sample_L2(1,:)-binary_sample_L2(3,:)), mean(stochastic_sample_L2(1,:)-binary_sample_L2(3,:)));
end
if conv2D_L3 == 1
    fprintf('L3 error standard deviation = %f, L3 mean error = %f \n', ...
    std(stochastic_sample_L3(1,:)-binary_sample_L3(3,:)), mean(stochastic_sample_L3(1,:)-binary_sample_L3(3,:)));
end
if conv2D_L4 == 1
    fprintf('L4 error standard deviation = %f, L4 mean error = %f \n', ...
    std(stochastic_sample_L4(1,:)-binary_sample_L4(3,:)), mean(stochastic_sample_L4(1,:)-binary_sample_L4(3,:)));
end
if conv2D_L5 == 1
    fprintf('L5 error standard deviation = %f, L5 mean error = %f \n', ...
    std(stochastic_sample_L5(1,:)-binary_sample_L5(3,:)), mean(stochastic_sample_L5(1,:)-binary_sample_L5(3,:)));
end
%function definition starts here
%//////////////////////////////////////////////////////////////////////////
function output = setup_LFSR(LFSR8, runs)
    LFSR8_val = zeros(1,runs, 'uint8');     LFSR8_wbg = zeros(1,runs, 'uint8');
    for i = 1:runs %setup LFSR registry and WBG of LFSR8
        LFSR8_val(i)=LFSR8(i,1)*128+LFSR8(i,2)*64+LFSR8(i,3)*32+LFSR8(i,4)*16+LFSR8(i,5)*8+LFSR8(i,6)*4+LFSR8(i,7)*2+LFSR8(i,8)*1;
        if LFSR8_val(i) <128
            if LFSR8_val(i) <64
                if LFSR8_val(i) <32
                    if LFSR8_val(i) <16
                        if LFSR8_val(i) <8
                            if LFSR8_val(i) <4
                                if LFSR8_val(i) <2
                                    if LFSR8_val(i) <1
                                        LFSR8_wbg(i) = cast(0,'uint8');
                                    else
                                        LFSR8_wbg(i) = cast(1,'uint8');
                                    end
                                else
                                    LFSR8_wbg(i) = cast(2,'uint8');
                                end
                            else
                                LFSR8_wbg(i) = cast(4,'uint8');
                            end
                        else
                            LFSR8_wbg(i) = cast(8,'uint8');
                        end
                    else
                        LFSR8_wbg(i) = cast(16,'uint8');
                    end
                else
                    LFSR8_wbg(i) = cast(32,'uint8');
                end
            else
                LFSR8_wbg(i) = cast(64,'uint8');
            end
        else
            LFSR8_wbg(i) = cast(128,'uint8');
        end
    end
    output = LFSR8_wbg;
end
function output = rand_sel_21(weight, wbg1, runs)
    sel = zeros(1,runs);
    pB0 = weight(2);
    for i = 1:runs
        b0_sc = bitand(pB0,wbg1(i)) > 0; 
        if b0_sc == 0
                sel(i) = 0;
        else
                sel(i) = 1;
        end
    end
    output = sel;
end
function output = rand_sel_41(weight, wbg1, wbg2, runs)
    sel = zeros(1,runs);
    pB1 = weight(3)+ weight(4);
    pB0_0 = round(weight(2)/(weight(1)+weight(2))*255);
    pB0_1 = round(weight(4)/(weight(3)+weight(4))*255);
    for i = 1:runs
        b0_0_sc = bitand(pB0_0,wbg1(i)) > 0; 
        b0_1_sc = bitand(pB0_1,wbg1(i)) > 0; 
        b1_sc = bitand(pB1,wbg2(i)) > 0;
        
        if b1_sc == 0
            if b0_0_sc == 0
                sel(i) = 0;
            else
                sel(i) = 1;
            end
        else
            if b0_1_sc == 0
                sel(i) = 2;
            else
                sel(i) = 3;
            end
        end
    end
    output = sel;
end
function output = rand_sel_61(weight, wbg1, wbg2, wbg3, runs)
    sel = zeros(1,runs);
    pB2 = weight(5)+ weight(6);
    pB1_0 = round((weight(3)+weight(4))/(weight(1)+weight(2)+weight(3)+weight(4))*255);
    pB1_1 = 0;
    pB0_0 = round(weight(2)/(weight(1)+weight(2))*255);
    pB0_1 = round(weight(4)/(weight(3)+weight(4))*255);
    pB0_2 = round(weight(6)/(weight(5)+weight(6))*255);
    for i = 1:runs
        b0_0_sc = bitand(pB0_0,wbg1(i)) > 0; 
        b0_1_sc = bitand(pB0_1,wbg1(i)) > 0; 
        b0_2_sc = bitand(pB0_2,wbg1(i)) > 0; 
        b1_0_sc = bitand(pB1_0,wbg2(i)) > 0; 
        b1_1_sc = bitand(pB1_1,wbg2(i)) > 0; 
        b2_sc = bitand(pB2,wbg3(i)) > 0;
        if b2_sc == 0
            if b1_0_sc == 0
                if b0_0_sc == 0
                    sel(i) = 0;
                else
                    sel(i) = 1;
                end
            else
                if b0_1_sc == 0
                    sel(i) = 2;
                else
                    sel(i) = 3;
                end
            end
        else
            if b1_1_sc == 0
                if b0_2_sc == 0
                    sel(i) = 4;
                else
                    sel(i) = 5;
                end
            end
        end
    end
    output = sel;
end
function output = rand_sel_81(weight, wbg1, wbg2, wbg3, runs)
    sel = zeros(1,runs);
    pB2 = weight(5)+ weight(6)+weight(7)+ weight(8);
    pB1_0 = round((weight(3)+weight(4))/(weight(1)+weight(2)+weight(3)+weight(4))*255);
    pB1_1 = round((weight(7)+weight(8))/(weight(5)+weight(6)+weight(7)+weight(8))*255);
    pB0_0 = round(weight(2)/(weight(1)+weight(2))*255);
    pB0_1 = round(weight(4)/(weight(3)+weight(4))*255);
    pB0_2 = round(weight(6)/(weight(5)+weight(6))*255);
    pB0_3 = round(weight(8)/(weight(7)+weight(8))*255);
    for i = 1:runs
        b0_0_sc = bitand(pB0_0,wbg1(i)) > 0; 
        b0_1_sc = bitand(pB0_1,wbg1(i)) > 0; 
        b0_2_sc = bitand(pB0_2,wbg1(i)) > 0; 
        b0_3_sc = bitand(pB0_3,wbg1(i)) > 0;
        b1_0_sc = bitand(pB1_0,wbg2(i)) > 0; 
        b1_1_sc = bitand(pB1_1,wbg2(i)) > 0; 
        b2_sc = bitand(pB2,wbg3(i)) > 0;
        if b2_sc == 0
            if b1_0_sc == 0
                if b0_0_sc == 0
                    sel(i) = 0;
                else
                    sel(i) = 1;
                end
            else
                if b0_1_sc == 0
                    sel(i) = 2;
                else
                    sel(i) = 3;
                end
            end
        else
            if b1_1_sc == 0
                if b0_2_sc == 0
                    sel(i) = 4;
                else
                    sel(i) = 5;
                end
            else
                if b0_3_sc == 0
                    sel(i) = 6;
                else
                    sel(i) = 7;
                end
            end
        end
    end
    output = sel;
end
function output = schedule_mux(weight, bias)
    neg = []; pos = [];
    [ws,idx] = sort(reshape(weight.',[1,size(weight(:))]));
    list = [idx;ws];
    for i = 1:size(list,2)
        if list(2,i)<0
            neg = [neg,list(:,i)];
        else
            pos = [pos,list(:,i)];
        end
    end
    neg = abs(neg);
    pos = flip(pos,2);
    %add bias information
    if bias<0
        neg = [neg,[-1;abs(bias)]];
    else
        pos = [pos,[-1;abs(bias)]];
    end
    %bias input will be indexed as '-1'
    %zero input will be indexed as '0'
    %start scheduling
    sche_n = zeros(2,8,ceil(sum(neg(2,:)/255)));    %preallocate slot
    sche_p = zeros(2,8,ceil(sum(pos(2,:)/255)));    %preallocate slot
    %negative scheduling starts here
    MUXc = 1; flag = 0; clk = 255; next = 0; swap = [];
    for i = 1:size(neg,2)   %allocate all weight/bias into the array
        clk = clk - neg(2,i);   
        if clk>0           %if still got clock left 
            sche_n(:,i-flag,MUXc) = neg(:,i);
        elseif clk == 0     %if exactly no clock left
            sche_n(:,i-flag,MUXc) = neg(:,i);
            MUXc = MUXc+1; clk = 255; flag = i;
        else                %if running out of timing slot
            clk = clk + neg(2,i);           %restore last clock track 
            %check for next weight whether or not got smaller wheight that 
            %can be squeezed to the timing requirement
            for j = i+1:size(neg,2)         %run requirement check
                if j>size(neg,2)            %ignore last element check
                    break
                end
                if neg(2,j)>clk             %if next weight still bigger than required timing
                    continue
                else
                    swap = neg(:,i:j);      %record track length
                    swap(:,end) = [];       %remove last one
                    neg(:,i) = neg(:,j);    %replace current with found one
                    neg(:,i+1:j) = swap;    %reconstruct the new pos                    
                    next = 1;               %flag to next schedule
                    clk = clk - neg(2,i);
                    sche_n(:,i-flag,MUXc) = neg(:,i);
                    break
                end
            end
            if next == 1    %continue to next schedule
                next = 0;   %reset the flag
                continue
            else            %add zero weight and move to next MUX
                sche_n(:,i-flag,MUXc) = [0;clk];%fill in last slot with zero weight
                MUXc = MUXc+1; clk = 255-neg(2,i); flag = i-1;
                sche_n(:,1,MUXc) = neg(:,i);    %fill in next MUX with current weight
            end
        end
        if i == size(neg,2) && clk>0        %if timing slot is not used up after scheduling
            sche_n(:,i-flag+1,MUXc) = [0;clk];   %fill in last slot with zero weight
        end
    end
    %positive scheduling starts here
    MUXc = 1; flag = 0; clk = 255; next = 0; swap = [];
    for i = 1:size(pos,2)   %allocate all weight/bias into the array
        clk = clk - pos(2,i);   
        if clk>0           %if still got clock left 
            sche_p(:,i-flag,MUXc) = pos(:,i);
        elseif clk == 0     %if exactly no clock left
            sche_p(:,i-flag,MUXc) = pos(:,i);
            MUXc = MUXc+1; clk = 255; flag = i;
        else                %if running out of timing slot
            clk = clk + pos(2,i);           %restore last clock track 
            %check for next weight whether or not got smaller wheight that 
            %can be squeezed to the timing requirement
            for j = i+1:size(pos,2)         %run requirement check
                if j>size(pos,2)            %ignore last element check
                    break
                end
                if pos(2,j)>clk             %if next weight still bigger than required timing
                    continue
                else
                    swap = pos(:,i:j);      %record track length
                    swap(:,end) = [];       %remove last one
                    pos(:,i) = pos(:,j);    %replace current with found one
                    pos(:,i+1:j) = swap;    %reconstruct the new pos                    
                    next = 1;               %flag to next move
                    clk = clk - pos(2,i);
                    sche_p(:,i-flag,MUXc) = pos(:,i);
                    break
                end
            end
            if next == 1
                next = 0;   %continue the loop
                continue
            else            %add zero weight and move to next MUX
                sche_p(:,i-flag,MUXc) = [0;clk];%fill in last slot with zero weight
                MUXc = MUXc+1; clk = 255-pos(2,i); flag = i-1;
                sche_p(:,1,MUXc) = pos(:,i);    %fill in next MUX with current weight
            end
        end
        if i == size(pos,2) && clk>0        %if timing slot is not used up after scheduling
            sche_p(:,i-flag+1,MUXc) = [0;clk];   %fill in last slot with zero weight
        end
    end
    output.n = sche_n;
    output.p = sche_p;
end
function output = optimize_clk(weight, LFSR8)   %the runlength must be set to 255
    runs = 255; choice = [];
    sel_count = zeros(1,size(weight,2));
    list_mse = zeros(runs-1,1);
    LFSR8 = LFSR8(1:runs,:);
    LFSR8a = flip(LFSR8,2);
    LFSR8_wbga = setup_LFSR(LFSR8a, runs);
    if size(sel_count,2) > 2    %if MUX is more than 2 input
        for cycle = 1:runs-1
            LFSR8b = circshift(LFSR8,cycle); %(x,cycle)
            LFSR8c = flip(LFSR8b,2);
            LFSR8_wbgb = setup_LFSR(LFSR8b, runs);
            LFSR8_wbgc = setup_LFSR(LFSR8c, runs);          
            if size(sel_count,2) == 4
                sel = rand_sel_41(weight, LFSR8_wbga, LFSR8_wbgb, runs);
            elseif size(sel_count,2) == 6
                sel = rand_sel_61(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runs);
            elseif size(sel_count,2) == 8
                sel = rand_sel_81(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runs);
            end
            for i = 1:size(sel_count,2)
                sel_count(i) = sum(sel(:) == i-1);
            end
            list_mse(cycle) =  immse(sel_count, weight);
        end  
        [error, sortIdx] = sort(list_mse);
        for i = 1:50
            if error(i) == error(1)
                choice = [choice, sortIdx(i)];
            end
        end
    else
        for i = 1:50
            choice = [choice, 0];
        end
    end
    output = choice;
end
function output = gen_sel(weight, LFSR8, runs, mode, optimize, char)
    word = "Finding the best 2nd LFSR circular shift clk for ";
    if mode == 1
        if optimize
            fprintf(strcat(word,char,"... "));
            best_shift_clk = optimize_clk(weight, LFSR8);
            if optimize>size(best_shift_clk,2)  %limit selection option
                optimize = size(best_shift_clk,2);
            end
            min_clk = best_shift_clk(optimize);
            fprintf('%d clk\n',min_clk);
            LFSR8a = flip(LFSR8);
            LFSR8b = circshift(LFSR8,min_clk); %(x,cycle)
            LFSR8c = flip(LFSR8b,2);
            LFSR8_wbga = setup_LFSR(LFSR8a, runs);
            LFSR8_wbgb = setup_LFSR(LFSR8b, runs);
            LFSR8_wbgc = setup_LFSR(LFSR8c, runs);
        else
            LFSR8a = flip(LFSR8);
            LFSR8b = circshift(LFSR8,10); %(x,cycle)
            LFSR8c = flip(LFSR8b,2);
            LFSR8_wbga = setup_LFSR(LFSR8a, runs);
            LFSR8_wbgb = setup_LFSR(LFSR8b, runs);
            LFSR8_wbgc = setup_LFSR(LFSR8c, runs);
        end
        mux_size = size(weight,2);
        if mux_size==2
            sel = rand_sel_21(weight, LFSR8_wbga, runs);
        elseif mux_size==4
            sel = rand_sel_41(weight, LFSR8_wbga, LFSR8_wbgb, runs);
        elseif mux_size==6
            sel = rand_sel_61(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runs);
        elseif mux_size==8
            sel = rand_sel_81(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runs);
        end
    elseif mode == 2    %if select MATLAB random permutation for testing
        sel = []; idx=1;
        for i = 1:255
            if weight(idx)-1 >=0
                sel = [sel, idx-1];
                weight(idx) = weight(idx)-1;
            else
                idx = idx+1;
                sel = [sel, idx-1];
                weight(idx) = weight(idx)-1;
            end
        end
        sel = sel(randperm(length(sel)));  %randomize muxn
        if runs>255
            selx = sel;
            for i = 2:round(runs/255)
                sel = [sel,selx];
            end
        end
    end
    output = sel;
end
function output = gen_feed_bin(to_feed, schedule, sel, runs)
    mux_feed = []; mux_in = [];
    for i = 1:size(schedule,2)
        if schedule(1,i) > 0          %if the 1st row index >0, it is the weight
            mux_feed = [mux_feed, cast(to_feed(schedule(1,i),schedule(2,i)),'uint8')];
        elseif schedule(2,i) == 2     %if the 2nd row index is 2, it is bias
            mux_feed = [mux_feed, cast(255,'uint8')];
        elseif schedule(2,i) == 3     %if the 2nd row index is 3, it is zero
            mux_feed = [mux_feed, cast(0,'uint8')];
        end
    end
    for i = 1:runs
        mux_in = [mux_in, mux_feed(sel(i)+1)];
    end
    output = mux_in;
end
function output = gen_feed_SC(to_feed, schedule, sel, SC_zero, runs)
    mux_feed = []; mux_in = [];
    for i = 1:size(schedule,2)
        if schedule(1,i) > 0        %if the 1st row index >0, it is the weight
            mux_feed = [mux_feed, to_feed(schedule(1,i),schedule(2,i),:)];
        elseif schedule(2,i) == 2     %if the 2nd row index is 2, it is bias
            mux_feed = [mux_feed, ones(1,1,runs)];
        elseif schedule(2,i) == 3     %if the 2nd row index is 3, it is zero
            mux_feed = [mux_feed, SC_zero];
        end
    end
    for i = 1:runs
        mux_in = [mux_in, mux_feed(1,sel(i)+1,i)];
    end
    output = mux_in;
end
function plot_SUC_response(record,layer,fig_id)
    figure(fig_id)
    tiledlayout(3,6)
    nexttile
    plot(record(1,:,1))   %plot SC real output
    title(strcat('SUC-',layer,' output 1'))
    nexttile
    plot(record(1,:,2))   %plot SC real output
    title(strcat('SUC-',layer,' output 2'))
    nexttile
    plot(record(1,:,3))   %plot SC real output
    title(strcat('SUC-',layer,' output 3'))
    nexttile
    plot(record(1,:,4))   %plot SC real output
    title(strcat('SUC-',layer,' output 4'))
    nexttile
    plot(record(1,:,5))   %plot SC real output
    title(strcat('SUC-',layer,' output 5'))
    nexttile
    plot(record(1,:,6))   %plot SC real output
    title(strcat('SUC-',layer,' output 6'))
    
    nexttile
    plot(record(2,:,1))   %plot SC real output
    title(strcat('SUC-',layer,' inter-weight 1'))
    nexttile
    plot(record(2,:,2))   %plot SC real output
    title(strcat('SUC-',layer,' inter-weight 2'))
    nexttile
    plot(record(2,:,3))   %plot SC real output
    title(strcat('SUC-',layer,' inter-weight 3'))
    nexttile
    plot(record(2,:,4))   %plot SC real output
    title(strcat('SUC-',layer,' inter-weight 4'))
    nexttile
    plot(record(2,:,5))   %plot SC real output
    title(strcat('SUC-',layer,' inter-weight 5'))
    nexttile
    plot(record(2,:,6))   %plot SC real output
    title(strcat('SUC-',layer,' inter-weight 6'))
    
    nexttile
    plot(record(3,:,1))   %plot SC real output
    title(strcat('SUC-',layer,' integral 1'))
    nexttile
    plot(record(3,:,2))   %plot SC real output
    title(strcat('SUC-',layer,' integral 2'))
    nexttile
    plot(record(3,:,3))   %plot SC real output
    title(strcat('SUC-',layer,' integral 3'))
    nexttile
    plot(record(3,:,4))   %plot SC real output
    title(strcat('SUC-',layer,' integral 4'))
    nexttile
    plot(record(3,:,5))   %plot SC real output
    title(strcat('SUC-',layer,' integral 5'))
    nexttile
    plot(record(3,:,6))   %plot SC real output
    title(strcat('SUC-',layer,' integral 6'))
    
    %plot the intermediate response of SUC adder of second layer
end