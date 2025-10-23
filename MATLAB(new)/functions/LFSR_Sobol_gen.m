runlength = 255;
LFSR8 = [true false false false false false false false];%initialize LFSR
for i = 1:runlength-1 %linear feedback shift register
    LFSR8(i+1,1) = xor(xor(LFSR8(i,8),LFSR8(i,6)),xor(LFSR8(i,5),LFSR8(i,4)));
    LFSR8(i+1,2) = LFSR8(i,1);    LFSR8(i+1,3) = LFSR8(i,2);    LFSR8(i+1,4) = LFSR8(i,3);
    LFSR8(i+1,5) = LFSR8(i,4);    LFSR8(i+1,6) = LFSR8(i,5);    LFSR8(i+1,7) = LFSR8(i,6);
    LFSR8(i+1,8) = LFSR8(i,7);
end
SOBOL8 = [false false false false false false false false];
for i = 1:runlength-1 %linear feedback shift register
    SOBOL8(i+1,8) = mod(i,2)>=1; SOBOL8(i+1,7) = mod(i,4)/2>=1;
    SOBOL8(i+1,6) = mod(i,8)/4>=1;SOBOL8(i+1,5) = mod(i,16)/8>=1;
    SOBOL8(i+1,4) = mod(i,32)/16>=1;SOBOL8(i+1,3) = mod(i,64)/32>=1;
    SOBOL8(i+1,2) = mod(i,128)/64>=1;SOBOL8(i+1,1) = mod(i,256)/128>=1;
end
SOBOL8 = flip(~SOBOL8,2);

WBC_LFSR = LFSR2WBC(LFSR8);
WBC_SOBOL = LFSR2WBC(SOBOL8);
subplot(2,1,1)
plot(WBC_LFSR)
xlim([0 255])
ylim([0 7])
subplot(2,1,2)
plot(WBC_SOBOL)
xlim([0 255])
ylim([0 7])

%%
function output = LFSR2WBC(x) % output will be int8 array of original size
    output = int8.empty(0,length(x));   %initialize array for speed
    for i = 1:length(x)                 %weighted binary converter to feed to MUX-based SNG
        output(i) = find(x(i,:),1)-1;   %return the position of the first non-zero element
    end
end