function out = f12_delta(in, R, Rl, Xl, Xarm, vhvdc, vgrid, pconu, pconl, sgrid, idcdif_ref, reiacsum_ref)

    syms vdcsum vdcdif revacsum imvacsum revacdif imvacdif reiacdif imiacdif reiacsum imiacsum idcdif idcsum
    state_variables = [vdcsum vdcdif revacsum imvacsum revacdif imvacdif reiacdif imiacdif reiacsum imiacsum idcdif idcsum];

    eqn(1) = vdcsum + (idcsum * R) - vhvdc;  %REAL
    eqn(2) = (idcdif * Rl) - vdcdif;         %IMAG
    eqn(3) = (2 * Xarm * imiacsum) - revacsum;       %REAL
    eqn(4) = (-2 * Xarm * reiacsum) - imvacsum;      %IMAG
    eqn(5) = revacdif + (imiacdif * (Xl + Xarm/2)) - (reiacdif*Rl) - real(vgrid);    %REAL
    eqn(6) = imvacdif - (reiacdif * (Xl + Xarm/2)) - (imiacdif*Rl) - imag(vgrid);    %IMAG
    eqn(7) = -idcdif + idcdif_ref;
    eqn(8) = -reiacsum + reiacsum_ref;
%     eqn(8) = (imvacsum*reiacsum) - (revacsum*imiacsum);
    eqn(9) = (real(vgrid) * reiacdif) + (imag(vgrid) * imiacdif) + (vdcdif * idcdif) - real(sgrid);
    eqn(10) = (imag(vgrid) * reiacdif) - (real(vgrid) * imiacdif) - imag(sgrid); %-imiacsum + imiacsum_ref;
    eqn(11) = ((vdcsum/2) * (idcsum + idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) - (revacsum * reiacsum) - (imvacdif * imiacsum) - pconu;
    eqn(12) = ((vdcsum/2) * (idcsum + idcdif/2)) - (revacdif * reiacdif/2) - (imvacdif * imiacdif/2) + (revacsum * reiacsum) + (imvacdif * imiacsum) + pconl;

    temp = jacobian(eqn, state_variables);
    out = subs(temp, state_variables, transpose(in));
    out = double(out);

end