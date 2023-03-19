function out = f12(in, R, Rl, Xl, Xarm, vhvdc, vgrid, pconu, pconl, sgrid, idcdif_ref, reiacsum_ref)

    vdcsum = in(1);
    vdcdif = in(2);         % 0
    revacsum = in(3);       % 0
    imvacsum = in(4);       % 0
    revacdif = in(5);       
    imvacdif = in(6);
    reiacdif = in(7);        %Iac of AC grid
    imiacdif = in(8);       
    reiacsum = in(9);        %Iac of DC grid - 0
    imiacsum = in(10);       % 0
    idcdif = in(11);         %Idc of AC grid - 0
    idcsum = in(12);         %Idc of DC grid

    out = zeros(12,1);

    out(1) = vdcsum + (idcsum * R) - vhvdc;  %REAL
    out(2) = vdcdif - (idcdif * Rl);       %IMAG
    out(3) = (2 * Xarm * imiacsum) - revacsum;       %REAL
    out(4) = (-2 * Xarm * reiacsum) - imvacsum;      %IMAG
    out(5) = revacdif + (imiacdif * (Xl + Xarm/2)) - (reiacdif*Rl) - real(vgrid);    %REAL
    out(6) = imvacdif - (reiacdif * (Xl + Xarm/2)) - (imiacdif*Rl) - imag(vgrid);    %IMAG
    out(7) = (real(vgrid) * reiacdif) + (imag(vgrid) * imiacdif) + (vdcdif * idcdif) - real(sgrid);
    out(8) = (imag(vgrid) * reiacdif) - (real(vgrid) * imiacdif) - imag(sgrid); %-imiacsum + imiacsum_ref;
    out(9) = ((vdcsum/2) * (idcsum + idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) - (revacsum * reiacsum) - (imvacdif * imiacsum) - pconu;
    out(10) = ((vdcsum/2) * (idcsum + idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) + (revacsum * reiacsum) + (imvacdif * imiacsum) + pconl;

    out(11) = -idcdif + idcdif_ref;
    out(12) = -reiacsum + reiacsum_ref;
%     out(12) = (imvacsum*reiacsum) - (revacsum*imiacsum);
end