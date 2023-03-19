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
    eqn(9) = ((vdcsum/2) * (idcsum + idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) - (revacsum * reiacsum) - (imvacdif * imiacsum) - pconu;
    eqn(10) = ((vdcsum/2) * (idcsum + idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) + (revacsum * reiacsum) + (imvacdif * imiacsum) + pconl;
    
    eqn(11) = -idcdif + idcdif_ref;
    eqn(12) = -reiacsum + reiacsum_ref;
%     eqn(12) = (imvacsum*reiacsum) - (revacsum*imiacsum);

    temp = jacobian(eqn, state_variables);
    out = subs(temp, state_variables, transpose(in));
    out = double(out);

end