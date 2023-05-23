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
    out(8) = (imag(vgrid) * reiacdif) - (real(vgrid) * imiacdif) - imag(sgrid);

    out(9) = (-vdcdif + vdcsum/2)*(idcdif/2 + idcsum) - revacdif*(reiacdif/2+reiacsum) - imvacdif*(imiacdif/2+imiacsum) + 0.5*revacsum*(reiacdif/2+reiacsum) + 0.5*imvacsum*(imiacdif/2+imiacsum) - pconu;
    out(10) = (vdcdif + vdcsum/2)*(-idcdif/2 + idcsum) + revacdif*(-reiacdif/2+reiacsum) + imiacdif*(-imiacdif/2+imiacsum) + 0.5*revacsum*(-revacdif/2+reiacsum) + 0.5*imvacsum*(-imiacdif/2+imiacsum) - pconl;

    out(11) = idcdif - idcdif_ref;
    out(12) = reiacsum - reiacsum_ref;



%     out(9) = ((vdcsum/2) * (idcsum - idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) - (revacsum * reiacsum) - (imvacdif * imiacsum) - pconu;
%     out(10) = ((vdcsum/2) * (idcsum - idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) + (revacsum * reiacsum) + (imvacdif * imiacsum) - pconl;
%     out(12) = (imvacsum*reiacsum) - (revacsum*imiacsum);
end

% function out = f12(in, R, Rl, Xl, Xarm, Vhvdc, Vgrid, Sgrid)
%     vdcu = in(1);
%     vdcl = in(2);
%     idcu = in(3);
%     idcl = in(4);
%     revacu = in(5);
%     imvacu = in(6);
%     revacl = in(7);
%     imvacl = in(8);
%     reiacu = in(9);
%     imiacu = in(10);
%     reiacl = in(11);
%     imiacl = in(12);
% 
%     out = zeros(12,1);
%     
%     out(1) = vdcu + (idcu*R/2) - (Vhvdc/2);
%     out(2) = vdcl + (idcl*R/2) - (Vhvdc/2);
%     out(3) = - revacu + (Xarm*imiacu) + (Xl*imiacu) - (Xl*imiacl) - real(Vgrid);
%     out(4) = - imvacu - (Xarm*reiacu) - (Xl*reiacu) + (Xl*reiacu) - imag(Vgrid);
%     out(5) = revacl - (Xarm*imiacl) + (Xl*imiacu) - (Xl*imiacl) - real(Vgrid);
%     out(6) = imvacl + (Xarm*reiacl) - (Xl*reiacu) + (Xl*reiacl) - imag(Vgrid);
%     out(7) = (real(Vgrid)*(reiacu-reiacl)) + (imag(Vgrid)*(imiacu-imiacl)) - real(Sgrid);
%     out(8) = (imag(Vgrid)*(reiacu-reiacl)) - (real(Vgrid)*(imiacu-imiacl)) - imag(Sgrid);
%     out(9) = (vdcu*idcu) - (revacu*reiacu) - (imvacu*imiacu);
%     out(10) = (vdcl*idcl) - (revacl*reiacl) - (imvacl*imiacl);
%     out(11) = Rl*(idcu-idcl);
%     out(12) = revacu - revacl;
% end

