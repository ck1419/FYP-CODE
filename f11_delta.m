%Single-Phase Single-Arm Equation Derivative Matrix
function out = f11_delta(x, Xarm, R, Vgrid_RE, Vgrid_IM)
    A1 = [x(3) x(4) x(1) x(2) x(6) x(5)];
    A2 = [1 0 0 Xarm 0 0];
    A3 = [0 1 -Xarm 0 0 0];
    A4 = [0 0 0 0 1 R];
    A5 = [0 0 Vgrid_RE Vgrid_IM 0 0];
    A6 = [0 0 Vgrid_IM -Vgrid_RE 0 0];
    out = [A1; A2; A3; A4; A5; A6];
end


