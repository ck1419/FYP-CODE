function out = f11(x, Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid)
    revac = x(1);
    imvac = x(2);
    reiac = x(3);
    imiac = x(4);
    vdc = x(5);
    idc = x(6);

    out = zeros(6,1);

    out(1) = - (revac*reiac) - (imvac*imiac) + (vdc*idc) - Pcon;
    out(2) = vdc + (idc*R) - Vhvdc;
    out(3) = revac - (Xarm*imiac) - Vgrid_RE;
    out(4) = imvac + (Xarm*reiac) - Vgrid_IM;
    out(5) = (Vgrid_RE*reiac) - (Vgrid_IM*imiac) - Pgrid;
    out(6) = (Vgrid_IM*reiac) + (Vgrid_RE*imiac) - Qgrid;
end

% %Single-Phase Single-Arm Equation Matrix
% function out = f11(x, Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid)
%     A1 = (x(5)*x(6)) + (x(1)*x(3)) + (x(2)*x(4)) - Pcon;
%     A2 = x(1) + (x(4)*Xarm) - Vgrid_RE;
%     A3 = x(2) - (x(3)*Xarm) - Vgrid_IM;
%     A4 = -x(5) + (x(6)*R) - Vhvdc;
%     A5 = (x(3)*Vgrid_RE) + (x(4)*Vgrid_IM) - Pgrid;
%     A6 = (x(3)*Vgrid_IM) - (x(4)*Vgrid_RE) - Qgrid;
%     out = [A1; A2; A3; A4; A5; A6];
% end


% %Single-Phase Single-Arm Equation Matrix
% function out = f11(x, Pcon, Xarm, R, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid)
%     A1 = (x(5)*x(6)) + (x(1)*x(3)) + (x(2)*x(4)) - Pcon;
%     A2 = x(1) + (x(4)*Xarm) - Vgrid_RE;
%     A3 = x(2) - (x(3)*Xarm) - Vgrid_IM;
%     A4 = x(5) + (x(6)*R) - Vhvdc;
%     A5 = (x(3)*Vgrid_RE) + (x(4)*Vgrid_IM) - Pgrid;
%     A6 = (x(3)*Vgrid_IM) - (x(4)*Vgrid_RE) - Qgrid;
%     out = [A1; A2; A3; A4; A5; A6];
% end
