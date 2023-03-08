function out = f11_delta(x, Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid)
    syms revac imvac reiac imiac vdc idc
    state_variables = [revac imvac reiac imiac vdc idc];

    eqn(1) = - (revac*reiac) - (imvac*imiac) + (vdc*idc) - Pcon;
    eqn(2) = vdc + (idc*R) - Vhvdc;
    eqn(3) = revac - (Xarm*imiac) - Vgrid_RE;
    eqn(4) = imvac + (Xarm*reiac) - Vgrid_IM;
    eqn(5) = (Vgrid_RE*reiac) - (Vgrid_IM*imiac) - Pgrid;
    eqn(6) = (Vgrid_IM*reiac) + (Vgrid_RE*imiac) - Qgrid;

    temp = jacobian(eqn, state_variables);
    out = subs(temp, state_variables, transpose(x));
    out = double(out);
end

% %Single-Phase Single-Arm Equation Derivative Matrix
% function out = f11_delta(x, Xarm, R, Vgrid_RE, Vgrid_IM)
%     A1 = [x(3) x(4) x(1) x(2) x(6) x(5)];
%     A2 = [1 0 0 Xarm 0 0];
%     A3 = [0 1 -Xarm 0 0 0];
%     A4 = [0 0 0 0 -1 R];
%     A5 = [0 0 Vgrid_RE Vgrid_IM 0 0];
%     A6 = [0 0 Vgrid_IM -Vgrid_RE 0 0];
%     out = [A1; A2; A3; A4; A5; A6];
% end
% 

% %Single-Phase Single-Arm Equation Derivative Matrix
% function out = f11_delta(x, Xarm, R, Vgrid_RE, Vgrid_IM)
%     A1 = [x(3) x(4) x(1) x(2) x(6) x(5)];
%     A2 = [1 0 0 Xarm 0 0];
%     A3 = [0 1 -Xarm 0 0 0];
%     A4 = [0 0 0 0 1 R];
%     A5 = [0 0 Vgrid_RE Vgrid_IM 0 0];
%     A6 = [0 0 Vgrid_IM -Vgrid_RE 0 0];
%     out = [A1; A2; A3; A4; A5; A6];
% end

