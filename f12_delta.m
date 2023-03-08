function out = f12_delta(in, R, Rl, Xl, Xarm, vhvdc, vgrid, pconu, pconl, sgrid, idcac_ref, reiacdc_ref)

    syms vdcsum vdcdif revacsum imvacsum revacdif imvacdif reiacac imiacac reiacdc imiacdc idcac idcdc
    state_variables = [vdcsum vdcdif revacsum imvacsum revacdif imvacdif reiacac imiacac reiacdc imiacdc idcac idcdc];

    eqn(1) = vdcsum + (idcdc * R) - vhvdc;  %REAL
    eqn(2) = (idcac * Rl) - vdcdif;         %IMAG
    eqn(3) = (2 * Xarm * imiacdc) - revacsum;       %REAL
    eqn(4) = (-2 * Xarm * reiacdc) - imvacsum;      %IMAG
    eqn(5) = revacdif + (imiacac * (Xl + Xarm/2)) - (reiacac*Rl) - real(vgrid);    %REAL
    eqn(6) = imvacdif - (reiacac * (Xl + Xarm/2)) - (imiacac*Rl) - imag(vgrid);    %IMAG
    eqn(7) = -idcac + idcac_ref;
    eqn(8) = -reiacdc + reiacdc_ref;
    eqn(9) = (real(vgrid) * reiacac) + (imag(vgrid) * imiacac) + (vdcdif * idcac) - real(sgrid);
    eqn(10) = (imag(vgrid) * reiacac) - (real(vgrid) * imiacac) - imag(sgrid); %-imiacdc + imiacdc_ref;
    eqn(11) = ((vdcsum/2) * (idcdc + idcac/2)) - (revacdif * reiacac/2) - (imvacdif * imiacac/2) - (revacsum * reiacdc) - (imvacdif * imiacdc) - pconu;
    eqn(12) = ((vdcsum/2) * (idcdc + idcac/2)) - (revacdif * reiacac/2) - (imvacdif * imiacac/2) + (revacsum * reiacdc) + (imvacdif * imiacdc) + pconl;

    temp = jacobian(eqn, state_variables);
    out = subs(temp, state_variables, transpose(in));
    out = double(out);

end