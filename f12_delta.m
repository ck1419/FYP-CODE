
function out = f12_delta(in, R, Rl, Xl, Xarm, vgrid)
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

    out(1,:) = [                  1,                   0,                     0,                 R,                       0,                       0,                                  0,                       0,                     0,                     0,                       0,                                0];
    out(2,:) = [                  0,                   1,                   -Rl,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,                       0,                                0];
    out(3,:) = [                  0,                   0,                     0,                 0,                       0,                       0,                                  1,                       0,                     0,                     0,                     -Rl,                      Xarm/2 + Xl];
    out(4,:) = [                  0,                   0,                     0,                 0,                       0,                       0,                                  0,                       1,                     0,                     0,           - Xarm/2 - Xl,                              -Rl];
    out(5,:) = [                  0,                   0,                     0,                 0,                      -1,                       0,                                  0,                       0,                     0,                2*Xarm,                       0,                                0];
    out(6,:) = [                  0,                   0,                     0,                 0,                       0,                      -1,                                  0,                       0,               -2*Xarm,                     0,                       0,                                0];
    out(7,:) = [                  0,              idcdif,                vdcdif,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,             real(vgrid),                      imag(vgrid)];
    out(8,:) = [                  0,                   0,                     0,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,             imag(vgrid),                     -real(vgrid)];
    out(9,:) = [idcdif/4 + idcsum/2, - idcdif/2 - idcsum,   vdcsum/4 - vdcdif/2, vdcsum/2 - vdcdif, reiacdif/4 + reiacsum/2, - imiacdif/4 - imiacsum/2,            - reiacdif/2 - reiacsum, imiacdif/2 + imiacsum, revacsum/2 - revacdif,   imvacdif - imvacsum/2, revacsum/4 - revacdif/2,          imvacdif/2 - imvacsum/4];
    out(10,:) = [idcsum/2 - idcdif/4,   idcsum - idcdif/2, - vdcdif/2 - vdcsum/4, vdcdif + vdcsum/2, reiacsum/2 - revacdif/4,   imiacdif/4 - imiacsum/2, reiacsum - reiacdif/2 - revacsum/4,                     0, revacdif + revacsum/2, - imiacdif - imvacsum/2,             -revacdif/2, imiacdif - imiacsum + imvacsum/4];
    out(11,:) = [                  0,                   0,                     1,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,                       0,                                0];
    out(12,:) = [                  0,                   0,                     0,                 0,                       0,                       0,                                  0,                       0,                     1,                     0,                       0,                                0];
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
