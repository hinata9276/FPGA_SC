%this function accept binary sum of parallel counted streams and perform advanced ReLU 
%the output is the logical bitstream with the size of the input binary stream,
%and the input binary stream is the result from parallel counters
function output = SCbinaryrelu2(nStream,pStream,stateSize,EDT) 
    maxState = stateSize-1; 
    trigger = stateSize/2-1;
    state = trigger;
    runlength = size(pStream,2);
    residue = 0;
    output = zeros([1 runlength],'logical');
    if EDT > 0 && EDT <=255
        for i = 1:EDT
            posSum = pStream(i)+state;
            negSum = nStream(i)+residue;
            state = posSum - negSum;
            %if bitget(state,8)>0 %if negative overflow occured
            if state<0 %remove low-level operation overhead
                state = 0;
            elseif state>maxState %bitget(state,MSBBit+1)>0
                state = maxState;
            end
            residue = state>trigger; %bitget(state,MSBBit);
            output(i) = residue; %bitget(state,MSBBit)>0;
        end
    else
        error('Error. Expected 0<EDT<=255, cannot be less or more. 255 = EDT at the end of runlength')
    end
end