syms in R Rl Xl Xarm vhvdc vgrid pconu pconl sgrid idcdif_ref reiacsum_ref
syms vdcsum vdcdif idcdif idcsum revacsum imvacsum revacdif imvacdif reiacsum imiacsum reiacdif imiacdif 

state_variables = [vdcsum vdcdif idcdif idcsum revacsum imvacsum revacdif imvacdif reiacsum imiacsum reiacdif imiacdif];

eqn(1) = vdcsum + (idcsum * R) - vhvdc;  
eqn(2) = vdcdif - (idcdif * Rl);   
eqn(3) = (2 * Xarm * imiacsum) - revacsum;   
eqn(4) = (-2 * Xarm * reiacsum) - imvacsum;  
eqn(5) = revacdif + (imiacdif * (Xl + Xarm/2)) - (reiacdif*Rl) - real(vgrid);
eqn(6) = imvacdif - (reiacdif * (Xl + Xarm/2)) - (imiacdif*Rl) - imag(vgrid);
eqn(7) = (real(vgrid) * reiacdif) + (imag(vgrid) * imiacdif) + (vdcdif * idcdif) - real(sgrid);
eqn(8) = (imag(vgrid) * reiacdif) - (real(vgrid) * imiacdif) - imag(sgrid);

eqn(9) = (-vdcdif + vdcsum/2)*(idcdif/2 + idcsum) - revacdif*(reiacdif/2+reiacsum) - imvacdif*(imiacdif/2+imiacsum) + 0.5*revacsum*(reiacdif/2+reiacsum) + 0.5*imvacsum*(imiacdif/2+imiacsum) - pconu;
eqn(10) = (vdcdif + vdcsum/2)*(-idcdif/2 + idcsum) + revacdif*(-reiacdif/2+reiacsum) + imiacdif*(-imiacdif/2+imiacsum) + 0.5*revacsum*(-revacdif/2+reiacsum) + 0.5*imvacsum*(-imiacdif/2+imiacsum) - pconl;

eqn(11) = idcdif - idcdif_ref;
eqn(12) = reiacsum - reiacsum_ref;



temp = jacobian(eqn, state_variables);