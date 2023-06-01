%Extracts the sign and turn it into a text format with spaces
%used for output displays
function out = disp_sign(in)
    if sign(imag(in)) == -1
        out = ' - ';
    else
        out = ' + ';
    end
end