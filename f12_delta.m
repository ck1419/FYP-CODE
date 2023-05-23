
function out = f12_delta(in, R, Rl, Xl, Xarm, vhvdc, vgrid, pconu, pconl, sgrid, idcdif_ref, reiacsum_ref)

    syms vdcsum vdcdif revacsum imvacsum revacdif imvacdif reiacdif imiacdif reiacsum imiacsum idcdif idcsum
    state_variables = [vdcsum vdcdif revacsum imvacsum revacdif imvacdif reiacdif imiacdif reiacsum imiacsum idcdif idcsum];

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

%     eqn(9) = ((vdcsum/2) * (idcsum - idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) - (revacsum * reiacsum) - (imvacdif * imiacsum) - pconu;
%     eqn(10) = ((vdcsum/2) * (idcsum - idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) + (revacsum * reiacsum) + (imvacdif * imiacsum) - pconl;
%     eqn(12) = (imvacsum*reiacsum) - (revacsum*imiacsum);

    temp = jacobian(eqn, state_variables);
    out = subs(temp, state_variables, transpose(in));
    out = double(out);

end


% function out = f12_delta(in, R, Rl, Xl, Xarm, Vhvdc, Vgrid, Sgrid)
%     syms vdcu vdcl idcu idcl revacu imvacu revacl imvacl reiacu imiacu reiacl imiacl
%     state_variables = [vdcu vdcl idcu idcl revacu imvacu revacl imvacl reiacu imiacu reiacl imiacl];
%     
%     eqn(1) = vdcu + (idcu*R/2) - (Vhvdc/2);
%     eqn(2) = vdcl + (idcl*R/2) - (Vhvdc/2);
%     eqn(3) = - revacu + (Xarm*imiacu) + (Xl*imiacu) - (Xl*imiacl) - real(Vgrid);
%     eqn(4) = - imvacu - (Xarm*reiacu) - (Xl*reiacu) + (Xl*reiacu) - imag(Vgrid);
%     eqn(5) = revacl - (Xarm*imiacl) + (Xl*imiacu) - (Xl*imiacl) - real(Vgrid);
%     eqn(6) = imvacl + (Xarm*reiacl) - (Xl*reiacu) + (Xl*reiacl) - imag(Vgrid);
%     eqn(7) = (real(Vgrid)*(reiacu-reiacl)) + (imag(Vgrid)*(imiacu-imiacl)) - real(Sgrid);
%     eqn(8) = (imag(Vgrid)*(reiacu-reiacl)) - (real(Vgrid)*(imiacu-imiacl)) - imag(Sgrid);
%     eqn(9) = (vdcu*idcu) - (revacu*reiacu) - (imvacu*imiacu);
%     eqn(10) = (vdcl*idcl) - (revacl*reiacl) - (imvacl*imiacl);
%     eqn(11) = Rl*(idcu-idcl);
%     eqn(12) = revacu - revacl;
% 
%     temp = jacobian(eqn, state_variables);
%     out = subs(temp, state_variables, transpose(in));
%     out = double(out);
% end
