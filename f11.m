function out = f11(x, Pcon, Xarm, R, Rarm, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid)
    revac = x(1);
    imvac = x(2);
    reiac = x(3);
    imiac = x(4);
    vdc = x(5);
    idc = x(6);

    out = zeros(6,1);

    out(1) = - (revac*reiac) - (imvac*imiac) + (vdc*idc) - Pcon;
    out(2) = vdc + (idc*R) - Vhvdc;
    out(3) = revac - (Xarm*imiac) - (Rarm*reiac) - Vgrid_RE;
    out(4) = imvac + (Xarm*reiac) - (Rarm*imiac) - Vgrid_IM;
    out(5) = (Vgrid_RE*reiac) + (Vgrid_IM*imiac) - Pgrid;
    out(6) = (Vgrid_IM*reiac) - (Vgrid_RE*imiac) - Qgrid;
end
