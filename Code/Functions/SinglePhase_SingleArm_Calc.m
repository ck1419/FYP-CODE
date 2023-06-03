%does the Newton-Raphson process for the single-phase single-arm model
function final = SinglePhase_SingleArm_Calc(in, max_iteration, tolerance, Pcon, Xarm, R, Rarm, Vgrid_RE, Vgrid_IM, Vhvdc, Pgrid, Qgrid)

    %initializes the starting matrices
    x = zeros(6, max_iteration);
    x(:,1) = in;

    %converts 0s for faster calculations
    Pcon = pseudo_zero(Pcon);
    Vgrid_RE = pseudo_zero(Vgrid_RE);
    Vgrid_IM = pseudo_zero(Vgrid_IM);
    Pgrid = pseudo_zero(Pgrid);
    Qgrid = pseudo_zero(Qgrid);

    %Starts the iterations
    for n = 2:max_iteration

        %separates variables from x matrix for easier analysis
        revac = x(1,n-1);
        imvac = x(2,n-1);
        reiac = x(3,n-1);
        imiac = x(4,n-1);
        vdc = x(5,n-1);
        idc = x(6,n-1);

        %approximates idc, helps with converging speed massively
        if n == 3
            idc_candidates = roots([R -Vhvdc (revac*reiac + imvac*imiac)]);
            idc = min(idc_candidates);
        end

        %computes f(x)
        fx = zeros(6,1);
        fx(1) = (revac*reiac) + (imvac*imiac) + (vdc*idc) - Pcon;
        fx(2) = - vdc + (idc*R) - Vhvdc;
        fx(3) = revac + (Xarm*imiac) - (Rarm*reiac) - Vgrid_RE;
        fx(4) = imvac - (Xarm*reiac) - (Rarm*imiac) - Vgrid_IM;
        fx(5) = (Vgrid_RE*reiac) + (Vgrid_IM*imiac) - Pgrid;
        fx(6) = (Vgrid_IM*reiac) - (Vgrid_RE*imiac) - Qgrid;
    
        %computes the jacobian of f(x)
        jx = zeros(6,6);
        jx(1,:) = [reiac, imiac,    revac,     imvac, idc, vdc];
        jx(2,:) = [    0,     0,        0,         0,  -1,   R];
        jx(3,:) = [    1,     0,    -Rarm,      Xarm,   0,   0];
        jx(4,:) = [    0,     1,    -Xarm,     -Rarm,   0,   0];
        jx(5,:) = [    0,     0, Vgrid_RE,  Vgrid_IM,   0,   0];
        jx(6,:) = [    0,     0, Vgrid_IM, -Vgrid_RE,   0,   0];

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
    for n = 1:6
        if abs(final(n)) <= 1
            final(n) = 0;
        end
    end
end