function out = f12(in, R, Rl, Xl, Xarm, vhvdc, vgrid, pconu, pconl, sgrid, idcdif_ref, imiacsum_ref)

    vdcsum = in(1);
    vdcdif = in(2); 
    idcdif = in(3);  
    idcsum = in(4);  
    revacsum = in(5);   
    imvacsum = in(6); 
    revacdif = in(7);       
    imvacdif = in(8);
    reiacsum = in(9);      
    imiacsum = in(10);    
    reiacdif = in(11);  
    imiacdif = in(12);       

    out = zeros(12,1);

    out(1) = vdcsum + (idcsum * R) - vhvdc;  %REAL
    out(2) = vdcdif - (idcdif * Rl);       %IMAG
    out(3) = revacdif + (imiacdif * (Xl + Xarm/2)) - (reiacdif*Rl) - real(vgrid);    %REAL
    out(4) = imvacdif - (reiacdif * (Xl + Xarm/2)) - (imiacdif*Rl) - imag(vgrid);    %IMAG
    out(5) = (2 * Xarm * imiacsum) - revacsum;       %REAL
    out(6) = (-2 * Xarm * reiacsum) - imvacsum;      %IMAG
    out(7) = (real(vgrid) * reiacdif) + (imag(vgrid) * imiacdif) + (vdcdif * idcdif) - real(sgrid);
    out(8) = (imag(vgrid) * reiacdif) - (real(vgrid) * imiacdif) - imag(sgrid);
    out(9) = (-vdcdif + vdcsum/2)*(idcdif/2 + idcsum) - revacdif*(reiacdif/2+reiacsum) + imvacdif*(imiacdif/2+imiacsum) + 0.5*revacsum*(reiacdif/2+reiacsum) - 0.5*imvacsum*(imiacdif/2+imiacsum) - pconu;
    out(10) = (vdcdif + vdcsum/2)*(-idcdif/2 + idcsum) + revacdif*(-reiacdif/2+reiacsum) - imvacdif*(-imiacdif/2+imiacsum) + 0.5*revacsum*(-reiacdif/2+reiacsum) - 0.5*imvacsum*(-imiacdif/2+imiacsum) - pconl;
    out(11) = idcdif - idcdif_ref;
    out(12) = imiacsum - imiacsum_ref;
end


