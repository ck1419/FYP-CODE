%does the Newton-Raphson process for the single-phase two-arm model
function final = SinglePhase_TwoArm_Calc(in, max_iteration, tolerance, R, Rl, Xl, Xarm, vhvdc, vgrid_RE, vgrid_IM, pconu, pconl, pgrid, qgrid, idcdif_ref, imiacsum_ref)

    %initializes the starting matrices
    x = zeros(12, max_iteration);
    x(:,1) = in;

    %converts 0s for faster calculations
    pconu = pseudo_zero(pconu);
    pconl = pseudo_zero(pconl);
    vgrid_RE = pseudo_zero(vgrid_RE);
    vgrid_IM = pseudo_zero(vgrid_IM);
    pgrid = pseudo_zero(pgrid);
    qgrid = pseudo_zero(qgrid);
    idcdif_ref = pseudo_zero(idcdif_ref);
    imiacsum_ref = pseudo_zero(imiacsum_ref);

    %Starts the iterations
    for n = 2:max_iteration

        %separates variables from x matrix for easier analysis
        vdcsum = x(1,n-1);
        vdcdif = x(2,n-1);
        idcdif = x(3,n-1); 
        idcsum = x(4,n-1); 
        revacsum = x(5,n-1);  
        imvacsum = x(6,n-1);
        revacdif = x(7,n-1);      
        imvacdif = x(8,n-1);
        reiacsum = x(9,n-1);     
        imiacsum = x(10,n-1);  
        reiacdif = x(11,n-1);
        imiacdif = x(12,n-1);     

        %computes f(x)
        fx = zeros(12,1);
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

        %computes the jacobian of f(x)
        jx = zeros(12,12);
        jx(1,:) = [                  1,                   0,                     0,                 R,                       0,                       0,                                  0,                       0,                     0,                     0,                       0,                                0];
        jx(2,:) = [                  0,                   1,                   -Rl,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,                       0,                                0];
        jx(3,:) = [                  0,                   0,                     0,                 0,                       0,                       0,                                  1,                       0,                     0,                     0,                     -Rl,                      Xarm/2 + Xl];
        jx(4,:) = [                  0,                   0,                     0,                 0,                       0,                       0,                                  0,                       1,                     0,                     0,           - Xarm/2 - Xl,                              -Rl];
        jx(5,:) = [                  0,                   0,                     0,                 0,                      -1,                       0,                                  0,                       0,                     0,                2*Xarm,                       0,                                0];
        jx(6,:) = [                  0,                   0,                     0,                 0,                       0,                      -1,                                  0,                       0,               -2*Xarm,                     0,                       0,                                0];
        jx(7,:) = [                  0,              idcdif,                vdcdif,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,             vgrid_RE,                      vgrid_IM];
        jx(8,:) = [                  0,                   0,                     0,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,             vgrid_IM,                     -vgrid_RE];
        jx(9,:) = [idcdif/4 + idcsum/2, - idcdif/2 - idcsum,   vdcsum/4 - vdcdif/2, vdcsum/2 - vdcdif, reiacdif/4 + reiacsum/2, imiacdif/4 + imiacsum/2, - reiacdif/2 - reiacsum, - imiacdif/2 - imiacsum, revacsum/2 - revacdif, imvacsum/2 - imvacdif,   revacsum/4 - revacdif/2,   imvacsum/4 - imvacdif/2];
        jx(10,:) = [idcsum/2 - idcdif/4,   idcsum - idcdif/2, - vdcdif/2 - vdcsum/4, vdcdif + vdcsum/2, reiacsum/2 - reiacdif/4, imiacsum/2 - imiacdif/4,   reiacsum - reiacdif/2,   imiacsum - imiacdif/2, revacdif + revacsum/2, imvacdif + imvacsum/2, - revacdif/2 - revacsum/4, - imvacdif/2 - imvacsum/4];
        jx(11,:) = [                  0,                   0,                     1,                 0,                       0,                       0,                                  0,                       0,                     0,                     0,                       0,                                0];
        jx(12,:) = [                  0,                   0,                     0,                 0,                       0,                         0,                       0,                     0,                     0,                       1,                         0,                       0];

        %computes the current iteration results
        x(:,n) = x(:,n-1) - (jx^-1 * fx);

        %test for whether answer met tolerances
        if all((x(:,n)./x(:,n-1)) <= 1+tolerance) && all((x(:,n)./x(:,n-1)) >= 1-tolerance)
            final = x(:,n);
            break
        elseif n == max_iteration
            final = x(:,n);
        end
    end

    %forces all small values to converge to 0
    for n = 1:12
        if abs(final(n)) <= 1
            final(n) = 0;
        end
    end

end