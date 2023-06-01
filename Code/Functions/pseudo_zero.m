function out = pseudo_zero(in)
    out = in;
    if in == 0
        out = 1e-15;
    end
end