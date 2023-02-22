function out = disp_sign(in)
    if sign(imag(in)) == -1
        out = ' - ';
    else
        out = ' + ';
    end
end