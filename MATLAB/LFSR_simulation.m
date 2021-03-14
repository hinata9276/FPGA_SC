% Initialize all possible variables for performance run
a = 2; b = 3;       %test input for LFSR4 example
c = 12; d = 10;    %test input for LFSR8 example
initial_learn_rate1 = 1; final_learn_rate1 = 1; learn_step1 = 255;
e = 25; f =22;     %test input for LFSR6 example
initial_learn_rate2 = 1; final_learn_rate2 = 1; learn_step2 = 63;
fir_coef = 0.9;
runlength = 999;
filtered1 =  zeros(1,runlength+1);
filtered2 =  zeros(1,runlength+1);
filtered3 =  zeros(1,runlength+1);
LFSR4 = [1 0 0 0];      
LFSR8 = [1 0 0 0 0 0 0 0];
LFSR6 = [1 0 0 0 0 0];  
%for LFSR-4 study
FSM1 = 64/2; max1=FSM1*2-1 ;min1 = 0; half1=FSM1-1;      %subtraction based
FSM2 = 64/2; max2=FSM2*2-1 ;min2 = 0; half2=FSM2-1;      %XOR based
%for LFSR-8 study
FSM3 =64/2; max3=FSM3*2-1 ;min3 = 0; half3=FSM3-1;      %subtraction based
FSM4 =64/2; max4=FSM4*2-1 ;min4 = 0; half4=FSM4-1;      %XOR based
%for LFSR-6 study
FSM5 = 32/2; max5=FSM5*2-1 ;min5 = 0; half5=FSM5-1;      %subtraction based
FSM6 = 32/2; max6=FSM6*2-1 ;min6 = 0; half6=FSM6-1;      %XOR based
zero_p5 = 1;        %for subtraction
FSM1_rec =  zeros(1,runlength); logic1 = zeros(1,runlength); 
FSM2_rec =  zeros(1,runlength); logic2 = zeros(1,runlength);
FSM3_rec =  zeros(1,runlength); logic3 = zeros(255,255,runlength);
FSM4_rec =  zeros(1,runlength); logic4 = zeros(255,255,runlength);
%FSM3_rec =  zeros(1,runlength); logic3 = zeros(1,runlength);
%FSM4_rec =  zeros(1,runlength); logic4 = zeros(1,runlength);
FSM5_rec =  zeros(1,runlength); logic5 = zeros(63,63,runlength);
FSM6_rec =  zeros(1,runlength); logic6 = zeros(63,63,runlength);
%FSM5_rec =  zeros(1,runlength); logic5 = zeros(1,runlength);
%FSM6_rec =  zeros(1,runlength); logic6 = zeros(1,runlength);
conv_time3 = zeros(255,255);
conv_time4 = zeros(255,255);

LFSR4_val = zeros(1,runlength, 'uint8');     LFSR4_wbg = zeros(1,runlength, 'uint8');
LFSR4_perval = zeros(1,runlength, 'uint8');  LFSR4_perwbg = zeros(1,runlength, 'uint8'); 
LFSR4_permodval = zeros(1,runlength, 'uint8');LFSR4_permodwbg = zeros(1,runlength, 'uint8');

LFSR8_val = zeros(1,runlength, 'uint8');     LFSR8_wbg = zeros(1,runlength, 'uint8');
LFSR8_perval = zeros(1,runlength, 'uint8');  LFSR8_perwbg = zeros(1,runlength, 'uint8');
LFSR8_permodval = zeros(1,runlength, 'uint8');LFSR8_permodwbg = zeros(1,runlength, 'uint8');

LFSR6_val = zeros(1,runlength, 'uint8');     LFSR6_wbg = zeros(1,runlength, 'uint8');
LFSR6_perval = zeros(1,runlength, 'uint8');  LFSR6_perwbg = zeros(1,runlength, 'uint8');

cout1 = zeros(1,runlength, 'int8');          cout11 = zeros(1,runlength, 'int8');
SCout1 = zeros(1,runlength);                 SCout11 = zeros(1,runlength);     trueval1 = zeros(1,runlength);

cout2 = zeros(1,runlength, 'int8');          cout22 = zeros(1,runlength, 'int8');
SCout2 = zeros(1,runlength);                 SCout22 = zeros(1,runlength);       trueval2 = zeros(1,runlength);

cout3 = zeros(1,runlength, 'int8');          
SCout3 = zeros(1,runlength);                 
trueval3 = zeros(1,runlength);

% Linear feedback sequence for a runtime
for i = 1:runlength  %linear feedback shift register
    LFSR4(i+1,1) = xor(LFSR4(i,4),LFSR4(i,3));
    LFSR4(i+1,2) = LFSR4(i,1);    LFSR4(i+1,3) = LFSR4(i,2);    LFSR4(i+1,4) = LFSR4(i,3);
end
for i = 1:runlength %linear feedback shift register
    LFSR8(i+1,1) = xor(xor(LFSR8(i,8),LFSR8(i,6)),xor(LFSR8(i,5),LFSR8(i,4)));
    LFSR8(i+1,2) = LFSR8(i,1);    LFSR8(i+1,3) = LFSR8(i,2);    LFSR8(i+1,4) = LFSR8(i,3);
    LFSR8(i+1,5) = LFSR8(i,4);    LFSR8(i+1,6) = LFSR8(i,5);    LFSR8(i+1,7) = LFSR8(i,6);
    LFSR8(i+1,8) = LFSR8(i,7);
end
for i = 1:runlength  %linear feedback shift register
    LFSR6(i+1,1) = xor(LFSR6(i,6),LFSR6(i,5));
    LFSR6(i+1,2) = LFSR6(i,1);    LFSR6(i+1,3) = LFSR6(i,2);    LFSR6(i+1,4) = LFSR6(i,3);
    LFSR6(i+1,5) = LFSR6(i,4);    LFSR6(i+1,6) = LFSR6(i,5);  
end

% wire permutation
LFSR4_per = flip(LFSR4, 2); %permutation
LFSR4_permod = LFSR4_per;
LFSR8_per = LFSR8;
%LFSR8_per = flip(LFSR8, 2); %permutation
%LFSR8_per = circshift(LFSR8, 11); %circular shift
%LFSR8_permod = LFSR8_per;
%LFSR6_per = flip(LFSR6, 2); %permutation
LFSR6_per = LFSR6;
%LFSR6_per(:,3) = 1 - LFSR6_per(:,3);
%LFSR6_per = circshift(LFSR6, 20);
for i = 1:runlength %modified LFSR permutation for testing
    LFSR4_permod(i,1) = 1-LFSR4_per(i,1);
    LFSR8_permod(i,1) = 1-LFSR8_per(i,1);
end

% recover LFSR value for calculation purpose
for i = 1:runlength %get value of LFSR anf permutated LFSR
    LFSR4_val(i)=LFSR4(i,1)*8+LFSR4(i,2)*4+LFSR4(i,3)*2+LFSR4(i,4)*1;
    LFSR4_perval(i)=LFSR4_per(i,1)*8+LFSR4_per(i,2)*4+LFSR4_per(i,3)*2+LFSR4_per(i,4)*1;
    LFSR4_permodval(i)=LFSR4_permod(i,1)*8+LFSR4_permod(i,2)*4+LFSR4_permod(i,3)*2+LFSR4_permod(i,4)*1;
end
for i = 1:runlength  %get value
    LFSR8_val(i)=LFSR8(i,1)*128+LFSR8(i,2)*64+LFSR8(i,3)*32+LFSR8(i,4)*16+LFSR8(i,5)*8+LFSR8(i,6)*4+LFSR8(i,7)*2+LFSR8(i,8)*1;
    LFSR8_perval(i)=LFSR8_per(i,1)*128+LFSR8_per(i,2)*64+LFSR8_per(i,3)*32+LFSR8_per(i,4)*16+LFSR8_per(i,5)*8+LFSR8_per(i,6)*4+LFSR8_per(i,7)*2+LFSR8_per(i,8)*1;
    LFSR8_permodval(i)=LFSR8_permod(i,1)*128+LFSR8_permod(i,2)*64+LFSR8_permod(i,3)*32+LFSR8_permod(i,4)*16+LFSR8_permod(i,5)*8+LFSR8_permod(i,6)*4+LFSR8_permod(i,7)*2+LFSR8_permod(i,8)*1;
end
for i = 1:runlength %get value of LFSR anf permutated LFSR
    LFSR6_val(i)=LFSR6(i,1)*32+LFSR6(i,2)*16+LFSR6(i,3)*8+LFSR6(i,4)*4+LFSR6(i,5)*2+LFSR6(i,6)*1;
    LFSR6_perval(i)=LFSR6_per(i,1)*32+LFSR6_per(i,2)*16+LFSR6_per(i,3)*8+LFSR6_per(i,4)*4+LFSR6(i,5)*2+LFSR6(i,6)*1;
end

% get intermediate WBG value (WBG part1)
for i = 1:runlength  %WBG of LFSR4
    if LFSR4_val(i) <8
        if LFSR4_val(i) <4
            if LFSR4_val(i) <2
                if LFSR4_val(i) <1
                    LFSR4_wbg(i) = 0;
                else
                    LFSR4_wbg(i) = 1;
                end
            else
                LFSR4_wbg(i) = 2;
            end
        else
            LFSR4_wbg(i) = 4;
        end
    else
        LFSR4_wbg(i) = 8;
    end
end
for i = 1:runlength  %WBG of LFSR4_per
    if LFSR4_perval(i) <8
        if LFSR4_perval(i) <4
            if LFSR4_perval(i) <2
                if LFSR4_perval(i) <1
                    LFSR4_perwbg(i) = 0;
                else
                    LFSR4_perwbg(i) = 1;
                end
            else
                LFSR4_perwbg(i) = 2;
            end
        else
            LFSR4_perwbg(i) = 4;
        end
    else
        LFSR4_perwbg(i) = 8;
    end
end
for i = 1:runlength  %WBG of LFSR4_permod
    if LFSR4_permodval(i) <8
        if LFSR4_permodval(i) <4
            if LFSR4_permodval(i) <2
                if LFSR4_permodval(i) <1
                    LFSR4_permodwbg(i) = 0;
                else
                    LFSR4_permodwbg(i) = 1;
                end
            else
                LFSR4_permodwbg(i) = 2;
            end
        else
            LFSR4_permodwbg(i) = 4;
        end
    else
        LFSR4_permodwbg(i) = 8;
    end
end
for i = 1:runlength %WBG of LFSR8
    if LFSR8_val(i) <128
        if LFSR8_val(i) <64
            if LFSR8_val(i) <32
                if LFSR8_val(i) <16
                    if LFSR8_val(i) <8
                        if LFSR8_val(i) <4
                            if LFSR8_val(i) <2
                                if LFSR8_val(i) <1
                                    LFSR8_wbg(i) = 0;
                                else
                                    LFSR8_wbg(i) = 1;
                                end
                            else
                                LFSR8_wbg(i) = 2;
                            end
                        else
                            LFSR8_wbg(i) = 4;
                        end
                    else
                        LFSR8_wbg(i) = 8;
                    end
                else
                    LFSR8_wbg(i) = 16;
                end
            else
                LFSR8_wbg(i) = 32;
            end
        else
            LFSR8_wbg(i) = 64;
        end
    else
        LFSR8_wbg(i) = 128;
    end
end
for i = 1:runlength  %WBG of LFSR8_per
    if LFSR8_perval(i) <128
        if LFSR8_perval(i) <64
            if LFSR8_perval(i) <32
                if LFSR8_perval(i) <16
                    if LFSR8_perval(i) <8
                        if LFSR8_perval(i) <4
                            if LFSR8_perval(i) <2
                                if LFSR8_perval(i) <1
                                    LFSR8_perwbg(i) = 0;
                                else
                                    LFSR8_perwbg(i) = 1;
                                end
                            else
                                LFSR8_perwbg(i) = 2;
                            end
                        else
                            LFSR8_perwbg(i) = 4;
                        end
                    else
                        LFSR8_perwbg(i) = 8;
                    end
                else
                    LFSR8_perwbg(i) = 16;
                end
            else
                LFSR8_perwbg(i) = 32;
            end
        else
            LFSR8_perwbg(i) = 64;
        end
    else
        LFSR8_perwbg(i) = 128;
    end
end
for i = 1:runlength  %WBG of LFSR8_permod
    if LFSR8_permodval(i) <128
        if LFSR8_permodval(i) <64
            if LFSR8_permodval(i) <32
                if LFSR8_permodval(i) <16
                    if LFSR8_permodval(i) <8
                        if LFSR8_permodval(i) <4
                            if LFSR8_permodval(i) <2
                                if LFSR8_permodval(i) <1
                                    LFSR8_permodwbg(i) = 0;
                                else
                                    LFSR8_permodwbg(i) = 1;
                                end
                            else
                                LFSR8_permodwbg(i) = 2;
                            end
                        else
                            LFSR8_permodwbg(i) = 4;
                        end
                    else
                        LFSR8_permodwbg(i) = 8;
                    end
                else
                    LFSR8_permodwbg(i) = 16;
                end
            else
                LFSR8_permodwbg(i) = 32;
            end
        else
            LFSR8_permodwbg(i) = 64;
        end
    else
        LFSR8_permodwbg(i) = 128;
    end
end
for i = 1:runlength %WBG of LFSR6
    if LFSR6_val(i) <32
        if LFSR6_val(i) <16
            if LFSR6_val(i) <8
                if LFSR6_val(i) <4
                    if LFSR6_val(i) <2
                        if LFSR6_val(i) <1
                            LFSR6_wbg(i) = 0;
                        else
                            LFSR6_wbg(i) = 1;
                        end
                    else
                        LFSR6_wbg(i) = 2;
                    end
                else
                    LFSR6_wbg(i) = 4;
                end
            else
                LFSR6_wbg(i) = 8;
            end
        else
            LFSR6_wbg(i) = 16;
        end
    else
        LFSR6_wbg(i) = 32;
    end    
end
for i = 1:runlength %WBG of LFSR6_per
    if LFSR6_perval(i) <32
        if LFSR6_perval(i) <16
            if LFSR6_perval(i) <8
                if LFSR6_perval(i) <4
                    if LFSR6_perval(i) <2
                        if LFSR6_perval(i) <1
                            LFSR6_perwbg(i) = 0;
                        else
                            LFSR6_perwbg(i) = 1;
                        end
                    else
                        LFSR6_perwbg(i) = 2;
                    end
                else
                    LFSR6_perwbg(i) = 4;
                end
            else
                LFSR6_perwbg(i) = 8;
            end
        else
            LFSR6_perwbg(i) = 16;
        end
    else
        LFSR6_perwbg(i) = 32;
    end    
end

% time to get SC started!
%{
for i = 1:runlength  %SC on LFSR4
    sn1 = bitand(a,LFSR4_wbg(i)) > 0;
    sn2 = bitand(b,LFSR4_perwbg(i)) > 0;
    sn3 = bitand(b,LFSR4_permodwbg(i)) > 0;
    cout1(i) = bitand(sn1,sn2);
    cout11(i) = bitand(sn1,sn3);
    SCout1(i) = sum(cout1(:) == 1)/i;
    SCout11(i) = sum(cout11(:) == 1)/i;
    trueval1(i) = a/16 * b/16;
    filtered1(i+1) = filtered1(i)*fir_coef + SCout1(i)*(1-fir_coef);
    zero_p5 = 1 - zero_p5;              % subtraction based comparator
    if(zero_p5)
        if(sn1)
            FSM1 = FSM1 + 1;
        else
            FSM1 = FSM1 - 1;
        end
    else
        if(~sn2)
            FSM1 = FSM1 + 1;
        else
            FSM1 = FSM1 - 1;
        end
    end
    if FSM1> max1
        FSM1 = max1;
    end
    if FSM1< min1
        FSM1 = min1;
    end
    FSM1_rec(i) = FSM1;
    logic1(i) = FSM1 > half1;
    %{
    if i <128                           %experimental learning rate
        learn = 2;
    else
        learn = 1;
    end
    %}
    learn = 1;
    if(bitxor(sn1,sn2))                 %XOR based comparator
       if(sn1)
           FSM2 = FSM2 + learn;
       else
           FSM2 = FSM2 - learn;
       end
    end
    if FSM2> max2
        FSM2 = max2;
    end
    if FSM2< min2
        FSM2 = min2;
    end
    FSM2_rec(i) = FSM2;
    logic2(i) = FSM2 > half2;
end
%}
FSM3_buf = FSM3;
FSM4_buf = FSM4;

for p = 1:255
    row = uint8(p);
    for q = 1:255  
        col = uint8(q);
        FSM3 = FSM3_buf;
        FSM4 = FSM4_buf;

        for i = 1:runlength  %SC on LFSR8
            sn1 = bitand(row,LFSR8_wbg(i)) > 0;
            sn2 = bitand(col,LFSR8_perwbg(i)) > 0;
            %sn1 = bitand(c,LFSR8_wbg(i)) > 0;
            %sn2 = bitand(d,LFSR8_perwbg(i)) > 0;
            %cout2(i) = bitand(sn1,sn2);
            %SCout2(i) = sum(cout2(:) == 1)/i;
            %trueval2(i) = c/256 * d/256;
            
            zero_p5 = 1 - zero_p5;              % subtraction based comparator
            if(zero_p5)
                if(sn1)
                    FSM3 = FSM3 + 1;
                else
                    FSM3 = FSM3 - 1;
                end
            else
                if(~sn2)
                    FSM3 = FSM3 + 1;
                else
                    FSM3 = FSM3 - 1;
                end
            end
            if FSM3> max3
                FSM3 = max3;
            end
            if FSM3< min3
                FSM3 = min3;
            end
            FSM3_rec(i) = FSM3;
            %logic3(i) = FSM3 > half3;
            logic3(row,col,i) = FSM3 > half3;
            
            
            if i <=learn_step1                           %experimental learning rate
                learn = initial_learn_rate1;
            else
                learn = final_learn_rate1;
            end

            if(bitxor(sn1,sn2))                 %XOR based comparator
                if(sn1)
                    FSM4 = FSM4 + learn;
                else
                    FSM4 = FSM4 - learn;
                end
            end
            if FSM4> max4
                FSM4 = max4;
            end
            if FSM4< min4
                FSM4 = min4;
            end
            FSM4_rec(i) = FSM4;
            %logic4(i) = FSM4 > half4;
            logic4(row,col,i) = FSM4 > half4;
        end
    end
end

%{
str1 = 'SUB 64 state correlated';
v = VideoWriter(strcat('D:\a_',str1,'_surf.avi'));
open(v);
for k = 1:runlength
    surf(logic3(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);
v = VideoWriter(strcat('D:\a_',str1,'_flat.avi'));
open(v);
for k = 1:runlength
    imagesc(logic3(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);


str1 = 'XOR 64 state circshift EDT510';
v = VideoWriter(strcat('D:\b_',str1,'_surf.avi'));
open(v);
for k = 1:runlength
    surf(logic4(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);
v = VideoWriter(strcat('D:\b_',str1,'_flat.avi'));
open(v);
for k = 1:runlength
    imagesc(logic4(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);  
%}

for row = 1:255
    for col = 1:255
        for i = 1:runlength-1
            if(logic3(row,col,runlength-i)~=logic3(row,col,runlength-i+1))    %if previous logic is different 
                conv_time3(row,col)= runlength-i;%record current clock cycle
                break
            end
        end
    end
end
for row = 1:255
    for col = 1:255
        for i = 1:runlength-1
            if(logic4(row,col,runlength-i)~=logic4(row,col,runlength-i+1))    %if previous logic is different 
                conv_time4(row,col)= runlength-i;%record current clock cycle
                break
            end
        end
    end
end
fprintf('Mean convergence time for XOR FSM relative to subtractor FSM = %f \n', mean(conv_time4,'all')-mean(conv_time3,'all'));
fprintf('- number means shorter convergence time than subtractor FSM, otherwise means worse convergence');
relative = conv_time3-conv_time4;
relative(relative>127)=0;
relative(relative<-127)=0;
figure(1)
image(relative+127)            % +/-127 data points removed
colorbar
relative = conv_time3-conv_time4; 
figure(2)
image((relative/1000*128+128)) %all data points scaled
colorbar
%{
FSM5_buf = FSM5;
FSM6_buf = FSM6;
for p = 1:63                               
    row = uint8(p);
    for q = 1:63  
        col = uint8(q);
        FSM5 = FSM5_buf;
        FSM6 = FSM6_buf;

        for i = 1:runlength  %SC on LFSR6
            sn1 = bitand(row,LFSR6_wbg(i)) > 0;
            sn2 = bitand(col,LFSR6_perwbg(i)) > 0;
            %sn1 = bitand(c,LFSR6_wbg(i)) > 0;
            %sn2 = bitand(d,LFSR6_perwbg(i)) > 0;
            %cout3(i) = bitand(sn1,sn2);
            %SCout3(i) = sum(cout3(:) == 1)/i;
            %trueval3(i) = e/64 * f/64;
            
            zero_p5 = 1 - zero_p5;              % subtraction based comparator
            if(zero_p5)
                if(sn1)
                    FSM5 = FSM5 + 1;
                else
                    FSM5 = FSM5 - 1;
                end
            else
                if(~sn2)
                    FSM5 = FSM5 + 1;
                else
                    FSM5 = FSM5 - 1;
                end
            end
            if FSM5> max5
                FSM5 = max5;
            end
            if FSM5< min5
                FSM5 = min5;
            end
            FSM5_rec(i) = FSM5;
            %logic5(i) = FSM3 > half3;
            logic5(row,col,i) = FSM5 > half5;
            
            
            if i <=learn_step2                          %experimental learning rate
                learn = initial_learn_rate2;
            else
                learn = final_learn_rate2;
            end

            if(bitxor(sn1,sn2))                 %XOR based comparator
                if(sn1)
                    FSM6 = FSM6 + learn;
                else
                    FSM6 = FSM6 - learn;
                end
            end
            if FSM6> max6
                FSM6 = max6;
            end
            if FSM6< min6
                FSM6 = min6;
            end
            FSM6_rec(i) = FSM6;
            %logic6(i) = FSM4 > half4;
            logic6(row,col,i) = FSM6 > half6;
        end
    end
end
str1 = 'SUB 32 state 6bit correlated';
v = VideoWriter(strcat('D:\a_',str1,'_surf.avi'));
open(v);
for k = 1:runlength
    surf(logic5(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);
v = VideoWriter(strcat('D:\a_',str1,'_flat.avi'));
open(v);
for k = 1:runlength
    imagesc(logic5(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);


str1 = 'XOR 32 state 6bit correlated';
v = VideoWriter(strcat('D:\a_',str1,'_surf.avi'));
open(v);
for k = 1:runlength
    surf(logic6(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);
v = VideoWriter(strcat('D:\a_',str1,'_flat.avi'));
open(v);
for k = 1:runlength
    imagesc(logic6(:,:,k))
    str2 = strcat('SC cycle = ', num2str(k));    title({str1,str2})
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    frame = getframe(ax,rect); writeVideo(v,frame);
end
close(v);  
%}
% SC is done!
% plotting graph

%tiledlayout(2,2)

%{
nexttile
%plot(FSM1_rec)
plot(FSM5_rec)
title('Subtract FSM state 6-bit')
%legend({'y = original permutation','y = modded','y = true value'},'Location','northeast')

% Bottom plot
nexttile
%plot(FSM2_rec)
plot(FSM6_rec)
title('XOR FSM state 6-bit')
%}
%{
% Bottom plot
nexttile
%plot(FSM2_rec)
plot(FSM3_rec)
title('Subtract FSM state 8-bit')

% Bottom plot
nexttile
%plot(FSM2_rec)
plot(FSM4_rec)
title('XOR FSM state 8-bit')
%}
%/////////////////////////////////////////////
%{
nexttile
%plot(logic1)
plot(logic5)
axis([1 runlength -1 2])
title('Subtract FSM logic 6-bit')

nexttile
%plot(logic2)
plot(logic6)
axis([1 runlength -1 2])
title('XOR FSM logic 6-bit')
%}
%{
nexttile
%plot(logic1)
plot(logic3)
axis([1 runlength -1 2])
title('Subtract FSM logic 8-bit')

nexttile
%plot(logic2)
plot(logic4)
axis([1 runlength -1 2])
title('XOR FSM logic 8-bit')
%}