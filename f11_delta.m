function out = f11_delta(x, Xarm, R, Rarm, Vgrid_RE, Vgrid_IM)

    revac = x(1);
    imvac = x(2);
    reiac = x(3);
    imiac = x(4);
    vdc = x(5);
    idc = x(6);

    out(1,:) = [reiac, imiac,    revac,     imvac, idc, vdc];
    out(2,:) = [    0,     0,        0,         0,  -1,   R];
    out(3,:) = [    1,     0,    -Rarm,      Xarm,   0,   0];
    out(4,:) = [    0,     1,    -Xarm,     -Rarm,   0,   0];
    out(5,:) = [    0,     0, Vgrid_RE,  Vgrid_IM,   0,   0];
    out(6,:) = [    0,     0, Vgrid_IM, -Vgrid_RE,   0,   0];

%     syms revac imvac reiac imiac vdc idc
%     syms Xarm R Rarm Vgrid_RE Vgrid_IM
%     state_variables = [revac imvac reiac imiac vdc idc];
% 
%     out(1) = (revac*reiac) + (imvac*imiac) + (vdc*idc) - Pcon;
%     out(2) = - vdc + (idc*R) - Vhvdc;
%     out(3) = revac + (Xarm*imiac) - (Rarm*reiac) - Vgrid_RE;
%     out(4) = imvac - (Xarm*reiac) - (Rarm*imiac) - Vgrid_IM;
%     out(5) = (Vgrid_RE*reiac) + (Vgrid_IM*imiac) - Pgrid;
%     out(6) = (Vgrid_IM*reiac) - (Vgrid_RE*imiac) - Qgrid;
% 
%     jacobian(out, state_variables);

end
