%MUX selector programming optimization code
%Even if the scheduler could deliver zero error MUX selection frequency, the
%output of SUC MUX convolution is not guaranteed to be of higher accuracy.
%The reason is still under investigation, in addition to the failure of the
%Conv2D L2 SC pipelining operation.
clear;
runlength = 255;
LFSR8 = [1 0 0 0 0 0 0 0];
sel = [];

for i = 1:runlength-1 %linear feedback shift register
    LFSR8(i+1,1) = xor(xor(LFSR8(i,8),LFSR8(i,6)),xor(LFSR8(i,5),LFSR8(i,4)));
    LFSR8(i+1,2) = LFSR8(i,1);    LFSR8(i+1,3) = LFSR8(i,2);    LFSR8(i+1,4) = LFSR8(i,3);
    LFSR8(i+1,5) = LFSR8(i,4);    LFSR8(i+1,6) = LFSR8(i,5);    LFSR8(i+1,7) = LFSR8(i,6);
    LFSR8(i+1,8) = LFSR8(i,7);
end
LFSR8a = flip(LFSR8,2);
LFSR8b = circshift(LFSR8,124); %(x,cycle)
LFSR8c = flip(LFSR8b,2);
% get intermediate WBG value (WBG part1)
test = 3; %select test code
%test = 0, 4-to-1 MUX manual primitive scheduler
%       1, 4-to-1 MUX auto scheduler
%       2, 8-to-1 MUX auto scheduler
%       3, 6-to-1 MUX auto scheduler
error_rec = 1;  %to execute error recording to find the best LFSR config
                %This function is not available in test=0
if test == 0
    fprintf('Running 4-to-1 MUX manual code \n');
elseif test == 1
    fprintf('Running auto 4-to-1 MUX random scheduler code \n');
elseif test == 2
    fprintf('Running auto 8-to-1 MUX random scheduler code \n');
elseif test == 3
    fprintf('Running auto 6-to-1 MUX random scheduler code \n');
end

LFSR8_wbg = setup_LFSR(LFSR8, runlength);
LFSR8_wbga = setup_LFSR(LFSR8a, runlength);
LFSR8_wbgb = setup_LFSR(LFSR8b, runlength);
LFSR8_wbgc = setup_LFSR(LFSR8c, runlength);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%code for manual scheduling starts here

%////////////////////////////////////////////////////////////////////////////////////////
%testing primitive 4-to-1 mux code
%////////////////////////////////////////////////////////////////////////////////////////
if test == 0
    x1 = cast(116,'uint8');
    x2 = cast(11,'uint8');
    b1 = cast(47,'uint8');
    
    for i = 1:runlength
        x1_sc = bitand(x1,LFSR8_wbg(i)) > 0; 
        x2_sc = bitand(x2,LFSR8_wbg(i)) > 0; 
        b1_sc = bitand(b1,LFSR8_wbga(i)) > 0;
        
        if b1_sc == 0
            if x1_sc == 0
                sel = [sel,0];
            else
                sel = [sel,1];
            end
        else
            if x2_sc == 0
                sel = [sel,2];
            else
                sel = [sel,3];
            end
        end
    end
    fprintf('Number of element: \n');
    fprintf('Sel 0 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 0), weight(1), sum(sel(:) == 0) - weight(1));
    fprintf('Sel 1 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 1), weight(2), sum(sel(:) == 1) - weight(2));
    fprintf('Sel 2 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 2), weight(3), sum(sel(:) == 2) - weight(3));
    fprintf('Sel 3 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 3), weight(4), sum(sel(:) == 3) - weight(4));
    plot(sel)
    [h,p] = runstest(sel,median(sel));
    fprintf('Randomness test: h=%d, p=%d \n',h,p);
end

%////////////////////////////////////////////////////////////////////////////////////////
%testing automated 4-to-1 MUX select weight generation
%////////////////////////////////////////////////////////////////////////////////////////
if test == 1
    %weight = [85 85 69 16]; %in sequence of sel = 0,1,2,3
    %weight = [108 77 43 27];
    %weight = [108 77 27 43]; %testing to allocate zero weight to the end to ease scheduling coding
    %weight = [158 76 14 7];
    %weight = [76 14 7 158]; %testing to allocate zero weight to the end to ease scheduling coding
    %weight = [113 95 45 2];
    %weight = [170 69 16 0];
    %weight = [69 16 170 0];   %testing to allocate zero weight to the end to ease scheduling coding
    %weight = [121 107 27 0];
    %weight = [132 70 53 0];
    %weight = [70 53 132 0]; %testing to allocate zero weight to the end to ease scheduling coding
    %weight = [112 63 62 18];
    %weight = [157 60 38 0];
    %weight = [60 38 157 0];   %testing to allocate zero weight to the end to ease scheduling coding
    %weight = [105 108 27 15];   %no specific order
    weight = [108 105 42 0];
    %weight=weight(randperm(length(weight)));
    %The sequence of weight is insignificant in this scheduler method
    %the scheduler will assign the proper selection weight regardless of the
    %sequence or order of the weights
    if error_rec == 1
        fprintf('Finding the best LFSR circular shift number... \n');
        sel_count = zeros(1,4);
        list_mse = zeros(runlength-1,1);
        for cycle = 1:254
            LFSR8b = circshift(LFSR8,cycle); %(x,cycle)
            LFSR8_wbgb = setup_LFSR(LFSR8b, runlength);
            sel = rand_sel_41(weight, LFSR8_wbga, LFSR8_wbgb, runlength);
            [mse, sel_count] = sel_mse(sel,weight);
            list_mse(cycle) =  mse;
            %{
            for i = 1:4
                sel_count(i) = sum(sel(:) == i-1);
            end
            list_mse(cycle) =  immse(sel_count, weight);
            %}
        end
        [error, sortIdx] = sort(list_mse);
        minimum_error = error(1); best_shift_cycle = sortIdx(1);
        fprintf('Minimum error of %d achieved for circular shift of %d \n',minimum_error,best_shift_cycle);
        other_choice = [];
        for i = 2:60
            if error(i) == error(1)
                other_choice = [other_choice, sortIdx(i)];
            end
        end
        fprintf('Other shift cycles available = ');
        fprintf('%d,', other_choice);
        fprintf('\nOptimizing MUX sel output... \n');
        LFSR8b = circshift(LFSR8,best_shift_cycle); %(x,cycle)
        LFSR8_wbgb = setup_LFSR(LFSR8b, runlength);
    end
    sel = rand_sel_41(weight, LFSR8_wbga, LFSR8_wbgb, runlength);
    fprintf('Number of element: \n');
    fprintf('Sel 0 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 0), weight(1), sum(sel(:) == 0) - weight(1));
    fprintf('Sel 1 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 1), weight(2), sum(sel(:) == 1) - weight(2));
    fprintf('Sel 2 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 2), weight(3), sum(sel(:) == 2) - weight(3));
    fprintf('Sel 3 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 3), weight(4), sum(sel(:) == 3) - weight(4));
    plot(sel)
    [h,p] = runstest(sel,median(sel));
    fprintf('Randomness test: h=%d, p=%d \n',h,p);
    
end

%////////////////////////////////////////////////////////////////////////////////////////
%testing automated 8-to-1 MUX select weight generation
%////////////////////////////////////////////////////////////////////////////////////////
if test == 2
    %x(1,3)->69     %x(1,2)->31
    %zero  ->63     %x(1,3)->69
    %x(2,3)->44     %x(2,2)->17
    %x(1,2)->31     %x(2,3)->44
    %x(3,3)->23     %x(3,2)->8
    %x(2,2)->17     %x(3,3)->23
    %x(3,2)->8      %zero  ->63
    %weight = [69 63 44 31 23 17 8 0]; %in sequence of sel
    weight = [69 44 31 23 17 8 63 0]; %testing to allocate zero weight to the end to ease scheduling coding
    %weight = [31 69 17 44 8 23 63 0]; %no specific order
    %weight=weight(randperm(length(weight)));
    %The sequence of weight is insignificant in this scheduler method
    %the scheduler will assign the proper selection weight regardless of the
    %sequence or order of the weights
    if error_rec == 1
        fprintf('Finding the best LFSR circular shift number... \n');
        sel_count = zeros(1,8);
        list_mse = zeros(runlength-1,1);
        for cycle = 1:254
            LFSR8b = circshift(LFSR8,cycle); %(x,cycle)
            LFSR8c = flip(LFSR8b,2);
            LFSR8_wbgb = setup_LFSR(LFSR8b, runlength);
            LFSR8_wbgc = setup_LFSR(LFSR8c, runlength);
            sel = rand_sel_81(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runlength);
            [mse, sel_count] = sel_mse(sel,weight);
            list_mse(cycle) =  mse;
            %{
            for i = 1:8
                sel_count(i) = sum(sel(:) == i-1);
            end
            list_mse(cycle) =  immse(sel_count, weight);
            %}
        end
        [error, sortIdx] = sort(list_mse);
        minimum_error = error(1); best_shift_cycle = sortIdx(1);
        fprintf('Minimum error of %d achieved for circular shift of %d \n',minimum_error,best_shift_cycle);
        other_choice = [];
        for i = 2:50
            if error(i) == error(1)
                other_choice = [other_choice, sortIdx(i)];
            end
        end
        fprintf('Other shift cycles available = ');
        fprintf('%d,', other_choice);
        fprintf('\nOptimizing MUX sel output... \n');
        LFSR8b = circshift(LFSR8,best_shift_cycle); %(x,cycle)
        LFSR8c = flip(LFSR8b,2);
        LFSR8_wbgb = setup_LFSR(LFSR8b, runlength);
        LFSR8_wbgc = setup_LFSR(LFSR8c, runlength);
    else
        fprintf('Using default shift cycle number to check for MUX sel output... \n');
    end
    sel = rand_sel_81(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runlength);
    fprintf('Number of element: \n');
    fprintf('Sel 0 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 0), weight(1), sum(sel(:) == 0) - weight(1));
    fprintf('Sel 1 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 1), weight(2), sum(sel(:) == 1) - weight(2));
    fprintf('Sel 2 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 2), weight(3), sum(sel(:) == 2) - weight(3));
    fprintf('Sel 3 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 3), weight(4), sum(sel(:) == 3) - weight(4));
    fprintf('Sel 4 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 4), weight(5), sum(sel(:) == 4) - weight(5));
    fprintf('Sel 5 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 5), weight(6), sum(sel(:) == 5) - weight(6));
    fprintf('Sel 6 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 6), weight(7), sum(sel(:) == 6) - weight(7));
    fprintf('Sel 7 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 7), weight(8), sum(sel(:) == 7) - weight(8));
    plot(sel)
    [h,p] = runstest(sel,median(sel));
    fprintf('Randomness test: h=%d, p=%d \n',h,p);
    
end

%////////////////////////////////////////////////////////////////////////////////////////
%testing automated 6-to-1 MUX select weight generation
%////////////////////////////////////////////////////////////////////////////////////////
if test == 3
    %x(1,3)->69
    %zero  ->63
    %x(2,3)->44     
    %x(1,2)->31
    %x(3,3)->23
    %x(2,2)->17 
    %x(3,2)->8
    %newer finding shows that no specific order is required for the
    %scheduling. Arranging weight in series and try to fit as much weight as possible.
    %but it has difficulty in optimizing MUX timing allocation even though the
    %scheduling code could be a bit simpler.
    %positive sum   %negative sum
    %x(1,1)->48     x(1,3)->54
    %x(1,2)->44     x(2,3)->106
    %x(2,1)->60     x(3,2)->9
    %x(2,2)->25     x(3,3)->75
    %x(3,1)->60     bias  ->1
    %zero  ->18     zero  ->10
    %weight = [106 75 54 9 1 10]; %in sequence of sel
    %weight = [54 106 9 75 1 10]; %no specific order
    weight = [60 60 48 44 25 18]; %in sequence of sel
    %weight = [65 60 58 50 4 18]; %in sequence of sel
    %weight = [48 44 60 25 60 18]; %no specific order
    %weight = [67 61 58 52 7 10]; %in sequence of sel
    %weight = [94 50 39 31 17 24];
    %weight=weight(randperm(length(weight)));
    %The sequence of weight is insignificant in this scheduler method
    %the scheduler will assign the proper selection weight regardless of the
    %sequence or order of the weights
    if error_rec == 1
        fprintf('Finding the best LFSR circular shift numbering... \n');
        sel_count = zeros(1,6);
        list_mse = zeros(runlength-1,1);
        for cycle = 1:254
            LFSR8b = circshift(LFSR8,cycle); %(x,cycle)
            LFSR8c = flip(LFSR8b,2);
            LFSR8_wbgb = setup_LFSR(LFSR8b, runlength);
            LFSR8_wbgc = setup_LFSR(LFSR8c, runlength);
            sel = rand_sel_61(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runlength);
            [mse, sel_count] = sel_mse(sel,weight);
            list_mse(cycle) =  mse;
            %{
            for i = 1:6
                sel_count(i) = sum(sel(:) == i-1);
            end
            list_mse(cycle) = immse(sel_count, weight);
            %}
        end
        [error, sortIdx] = sort(list_mse);
        minimum_error = error(1); best_shift_cycle = sortIdx(1);
        fprintf('Minimum error of %d achieved for circular shift of %d \n',minimum_error,best_shift_cycle);
        other_choice = [];
        for i = 2:50
            if error(i) == error(1)
                other_choice = [other_choice, sortIdx(i)];
            end
        end
        fprintf('Other shift cycles available = ');
        fprintf('%d,', other_choice);
        fprintf('\nOptimizing MUX sel output... \n');
        %LFSR8b = circshift(LFSR8,best_shift_cycle); %(x,cycle)
        LFSR8b = circshift(LFSR8,best_shift_cycle); %(x,cycle)
        LFSR8c = flip(LFSR8b,2);
        LFSR8_wbgb = setup_LFSR(LFSR8b, runlength);
        LFSR8_wbgc = setup_LFSR(LFSR8c, runlength);
    else
        fprintf('Using default shift cycle number to check for MUX sel output... \n');
    end
    sel = rand_sel_61(weight, LFSR8_wbga, LFSR8_wbgb, LFSR8_wbgc, runlength);
    [mse, sel_count] = sel_mse(sel,weight);
    fprintf('Number of element: \n');
    fprintf('Sel 0 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 0), weight(1), sum(sel(:) == 0) - weight(1));
    fprintf('Sel 1 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 1), weight(2), sum(sel(:) == 1) - weight(2));
    fprintf('Sel 2 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 2), weight(3), sum(sel(:) == 2) - weight(3));
    fprintf('Sel 3 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 3), weight(4), sum(sel(:) == 3) - weight(4));
    fprintf('Sel 4 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 4), weight(5), sum(sel(:) == 4) - weight(5));
    fprintf('Sel 5 = %d, \t actual weight = %d, \t error = %d \n',sum(sel(:) == 5), weight(6), sum(sel(:) == 5) - weight(6));
    fprintf('Mean Square Error of select = %d \n',mse);
    plot(sel)
    [h,p] = runstest(sel.' ,median(sel));
    fprintf('Randomness test: h=%d, p=%d \n',h,p);
    disp(LFSR8(2,:))
    disp(LFSR8b(2,:))
    disp(LFSR8c(2,:))
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
function [output, sumsel] = sel_mse(sel,weight)
    if size(sel,2)>255
        sel = sel(1:255);
    end
    sel_size = size(weight,2);
    sel_count = zeros(1,sel_size);
    for i = 1:sel_size
        sel_count(i) = sum(sel(:) == i-1);
    end
    output =  immse(sel_count, weight);
    sumsel = sel_count;
end