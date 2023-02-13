function out = f12(in, R, Xl, Xarm, vhvdc, vgrid, pconu, pconl, sgrid, idcac_ref, reiacdc_ref)

    vdcsum = in(1);
    vdcdif = in(2);         % 0
    revacsum = in(3);       % 0
    imvacsum = in(4);       % 0
    revacdif = in(5);       
    imvacdif = in(6);
    reiacac = in(7);        %Iac of AC grid
    imiacac = in(8);       
    reiacdc = in(9);        %Iac of DC grid - 0
    imiacdc = in(10);       % 0
    idcac = in(11);         %Idc of AC grid - 0
    idcdc = in(12);         %Idc of DC grid

    out = zeros(12,1);

    out(1) = vdcsum + (idcdc * R) - vhvdc;  %REAL
    out(2) = (idcac * Xl) - vdcdif;         %IMAG
    out(3) = (2 * Xarm * imiacdc) - revacsum;       %REAL
    out(4) = (-2 * Xarm * reiacdc) - imvacsum;      %IMAG
    out(5) = revacdif + (imiacac * (Xl + Xarm/2)) - real(vgrid);    %REAL
    out(6) = imvacdif - (reiacac * (Xl + Xarm/2)) - imag(vgrid);    %IMAG
    out(7) = -idcac + idcac_ref;
    out(8) = -reiacdc + reiacdc_ref;
    out(9) = (real(vgrid) * reiacac) + (imag(vgrid) * imiacac) + (vdcdif * idcac) - real(sgrid);
    out(10) = (imag(vgrid) * reiacac) - (real(vgrid) * imiacac) - imag(sgrid); %-imiacdc + imiacdc_ref;
    out(11) = (vdcsum/4 * (2*idcdc + idcac)) - (revacdif * reiacac/2) + (imvacdif * imiacac/2) - (revacsum * reiacdc) + (imvacdif * imiacdc) - pconu;
    out(12) = (vdcsum/4 * (2*idcdc - idcac)) - (revacdif * reiacac/2) + (imvacdif * imiacac/2) + (revacsum * reiacdc) - (imvacdif * imiacdc) - pconl;
end