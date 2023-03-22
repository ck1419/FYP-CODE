function out = f11_delta(x, Pcon, Xarm, R, Rarm, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid)
    syms revac imvac reiac imiac vdc idc
    state_variables = [revac imvac reiac imiac vdc idc];

    eqn(1) = - (revac*reiac) - (imvac*imiac) + (vdc*idc) - Pcon;
    eqn(2) = vdc + (idc*R) - Vhvdc;
    eqn(3) = revac - (Xarm*imiac) - (Rarm*reiac) - Vgrid_RE;
    eqn(4) = imvac + (Xarm*reiac) - (Rarm*imiac) - Vgrid_IM;
    eqn(5) = (Vgrid_RE*reiac) + (Vgrid_IM*imiac) - Pgrid;
    eqn(6) = (Vgrid_IM*reiac) - (Vgrid_RE*imiac) - Qgrid;

    temp = jacobian(eqn, state_variables);
    out = subs(temp, state_variables, transpose(x));
    out = double(out);
end
