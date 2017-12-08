function array = shift_array(y_data,shift)
%This function will shift a given data array, y_data, by a given number of
%indices, shift. It will preserve shape by setting appending an appropriate
%amount of zeros.
    [~,index] = max(y_data);
    if shift == 0 || floor(shift) ~= shift
        error('Shift must an integer and non-zero.')
    end
    array = zeros(length(y_data),1);
    
    if shift < 0
%         if index < shift
%             error('Shift amount too great, you will lose the peak of the feature')
%         end
        array(1:end + shift) = y_data(1-shift:end);
    elseif shift > 0
%         if index+shift > length(y_data)
%             error('Shift amount too great, you will lose the peak of the feature')
%         end
        array(1+shift:end) = y_data(1:end-shift);
    end
end
    
