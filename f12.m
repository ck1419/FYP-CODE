function out = f12(in, R, Xl, vhvdc, vgrid, pconu, pconl)
    vdcsum = in(1);
    vdcdif = in(2);
    revacsum = in(3);
    imvacsum = in(4);
    revacdif = in(5);
    imvacdif = in(6);
    reiacac = in(7);        %Iac of AC grid
    imiacac = in(8);
    reiacdc = in(9);        %Iac of DC grid
    imiacdc = in(10);
    idcac = in(11);         %Idc of AC grid
    idcdc = in(12);         %Idc of DC grid

    out = zeros(12,1);

    out(1) = vdcsum + (idcdc * R) - vhvdc;  %REAL
    out(2) = (idcac * Xl) - vdcdif;         %IMAG
    out(3) = (2 * Xarm * imiacdc) - revacsum;       %REAL
    out(4) = (-2 * Xarmm * reiacdc) - imvacsum;     %IMAG
    out(5) = revacdif + (imiacac * (Xl + Xarm/2)) - real(vgrid);    %REAL
    out(6) = imvacdif - (reiacac * (Xl + Xarm/2)) - imag(vgrid);    %IMAG
    out(7) = -idcac;
    out(8) = -reiacdc;
    out(9) =  -imiacdc;
    out(10) = (real(vgrid) * reiacac) + (imag(vgrid) * imacdc) + (vdcdif * idcac);
    out(11) = (vdcsum/4 * (2*idcdc + idcac)) - (revacdif * reiacac/2) + (imvacdif * imiacac/2) - (revacsum * reiacdc) + (imvacdif * imiacdc) - pconu;
    out(12) = (vdcsum/4 * (2*idcdc - idcac)) - (revacdif * reiacac/2) + (imvacdif * imiacac/2) + (revacsum * reiacdc) - (imvacdif * imiacdc) - pconl;
end