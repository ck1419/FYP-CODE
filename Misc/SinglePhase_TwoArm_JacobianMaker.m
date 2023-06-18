clear
clc

syms in R Rl Xl Xarm vhvdc vgrid pconu pconl sgrid idcdif_ref imiacsum_ref vgrid_RE vgrid_IM pgrid qgrid
syms vdcsum vdcdif idcdif idcsum revacsum imvacsum revacdif imvacdif reiacsum imiacsum reiacdif imiacdif 

state_variables = [vdcsum vdcdif idcdif idcsum revacsum imvacsum revacdif imvacdif reiacsum imiacsum reiacdif imiacdif];

fx(1) = vdcsum + (idcsum * R) - vhvdc;
fx(2) = vdcdif - (idcdif * Rl);
fx(3) = revacdif + (imiacdif * (Xl + Xarm/2)) - (reiacdif*Rl) - vgrid_RE;
fx(4) = imvacdif - (reiacdif * (Xl + Xarm/2)) - (imiacdif*Rl) - vgrid_IM;
fx(5) = (2 * Xarm * imiacsum) - revacsum;
fx(6) = (-2 * Xarm * reiacsum) - imvacsum;
fx(7) = (vgrid_RE * reiacdif) + (vgrid_IM * imiacdif) + (vdcdif * idcdif) - pgrid;
fx(8) = (vgrid_IM * reiacdif) - (vgrid_RE * imiacdif) - qgrid;
fx(9) = (-vdcdif + vdcsum/2)*(idcdif/2 + idcsum) - revacdif*(reiacdif/2+reiacsum) - imvacdif*(imiacdif/2+imiacsum) + 0.5*revacsum*(reiacdif/2+reiacsum) + 0.5*imvacsum*(imiacdif/2+imiacsum) - pconu;
fx(10) = (vdcdif + vdcsum/2)*(-idcdif/2 + idcsum) + revacdif*(-reiacdif/2+reiacsum) + imvacdif*(-imiacdif/2+imiacsum) + 0.5*revacsum*(-reiacdif/2+reiacsum) + 0.5*imvacsum*(-imiacdif/2+imiacsum) - pconl;    
fx(11) = idcdif - idcdif_ref;
fx(12) = imiacsum - imiacsum_ref;

out = jacobian(fx, state_variables);

disp(out)