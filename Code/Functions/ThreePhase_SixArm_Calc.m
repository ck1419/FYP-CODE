%does the Newton-Raphson process for the three-phase six-arm model
function final = ThreePhase_SixArm_Calc(in, max_iteration, tolerance, R, Rl, Xl, Xarm, vhvdc, vgrid_RE, vgrid_IM, pconu, pconl, pgrid, qgrid, idcdif_ref, imiacsum_ref)

    %initializes the starting matrices
    x = zeros(36, max_iteration);
    x(:,1) = in;

    %applies phase difference to the Vgrids
    vgrid_RE_a = vgrid_RE;
    vgrid_IM_a = vgrid_IM;
    vgrid_a = vgrid_RE_a + (vgrid_IM_a);
    vgrid_b = vgrid_a * exp(1i*deg2rad(b_phase));
    vgrid_c = vgrid_a * exp(1i*deg2rad(c_phase));
    vgrid_RE_b = real(vgrid_b);
    vgrid_IM_b = imag(vgrid_b);
    vgrid_RE_c = real(vgrid_c);
    vgrid_IM_c = imag(vgrid_c);

    %applies phase difference to the Sgrids 
    pgrid_a = pgrid;
    qgrid_a = qgrid;
    sgrid_a = pgrid_a + (qgrid_a);
    sgrid_b = sgrid_a * exp(1i*deg2rad(b_phase));
    sgrid_c = sgrid_a * exp(1i*deg2rad(c_phase));
    pgrid_b = real(sgrid_b);
    qgrid_b = imag(sgrid_b);
    pgrid_c = real(sgrid_c);
    qgrid_c = imag(sgrid_c);

    %converts 0s for faster calculations for phase A
    pconu_a = pseudo_zero(pconu_a);
    pconl_a = pseudo_zero(pconl_a);
    vgrid_RE_a = pseudo_zero(vgrid_RE_a);
    vgrid_IM_a = pseudo_zero(vgrid_IM_a);
    pgrid_a = pseudo_zero(pgrid_a);
    qgrid_a = pseudo_zero(qgrid_a);
    idcdif_ref_a = pseudo_zero(idcdif_ref_a);
    imiacsum_ref_a = pseudo_zero(imiacsum_ref_a);

    %converts 0s for faster calculations for phase B
    pconu_b = pseudo_zero(pconu_b);
    pconl_b = pseudo_zero(pconl_b);
    vgrid_RE_b = pseudo_zero(vgrid_RE_b);
    vgrid_IM_b = pseudo_zero(vgrid_IM_b);
    pgrid_b = pseudo_zero(pgrid_b);
    qgrid_b = pseudo_zero(qgrid_b);
    idcdif_ref_b = pseudo_zero(idcdif_ref_b);
    imiacsum_ref_b = pseudo_zero(imiacsum_ref_b);

    %converts 0s for faster calculations for phase C
    pconu_c = pseudo_zero(pconu_c);
    pconl_c = pseudo_zero(pconl_c);
    vgrid_RE_c = pseudo_zero(vgrid_RE_c);
    vgrid_IM_c = pseudo_zero(vgrid_IM_c);
    pgrid_c = pseudo_zero(pgrid_c);
    qgrid_c = pseudo_zero(qgrid_c);
    idcdif_ref_c = pseudo_zero(idcdif_ref_c);
    imiacsum_ref_c = pseudo_zero(imiacsum_ref_c);

    %Starts the iterations
    for n = 2:max_iteration

        %separates variables from x matrix for easier analysis (phase A)
        vdcsum_a = x(1,n-1);
        vdcdif_a = x(2,n-1);
        idcdif_a = x(3,n-1); 
        idcsum_a = x(4,n-1); 
        revacsum_a = x(5,n-1);  
        imvacsum_a = x(6,n-1);
        revacdif_a = x(7,n-1);      
        imvacdif_a = x(8,n-1);
        reiacsum_a = x(9,n-1);     
        imiacsum_a = x(10,n-1);  
        reiacdif_a = x(11,n-1);
        imiacdif_a = x(12,n-1);    

        %separates variables from x matrix for easier analysis (phase B)
        vdcsum_b = x(13,n-1);
        vdcdif_b = x(14,n-1);
        idcdif_b = x(15,n-1); 
        idcsum_b = x(16,n-1); 
        revacsum_b = x(17,n-1);  
        imvacsum_b = x(18,n-1);
        revacdif_b = x(19,n-1);      
        imvacdif_b = x(20,n-1);
        reiacsum_b = x(21,n-1);     
        imiacsum_b = x(22,n-1);  
        reiacdif_b = x(23,n-1);
        imiacdif_b = x(24,n-1); 

        %separates variables from x matrix for easier analysis (phase C)
        vdcsum_c = x(25,n-1);
        vdcdif_c = x(26,n-1);
        idcdif_c = x(27,n-1); 
        idcsum_c = x(28,n-1); 
        revacsum_c = x(29,n-1);  
        imvacsum_c = x(30,n-1);
        revacdif_c = x(31,n-1);      
        imvacdif_c = x(32,n-1);
        reiacsum_c = x(33,n-1);     
        imiacsum_c = x(34,n-1);  
        reiacdif_c = x(35,n-1);
        imiacdif_c = x(36,n-1);    

        %computes f(x) for phase A
        fx = zeros(36,1);
        fx(1) = vdcsum_a + (idcsum_a * R) - vhvdc;
        fx(2) = vdcdif_a - (idcdif_a * Rl);
        fx(3) = revacdif_a + (imiacdif_a * (Xl + Xarm/2)) - (reiacdif_a*Rl) - vgrid_RE_a;
        fx(4) = imvacdif_a - (reiacdif_a * (Xl + Xarm/2)) - (imiacdif_a*Rl) - vgrid_IM_a;
        fx(5) = (2 * Xarm * imiacsum_a) - revacsum_a;
        fx(6) = (-2 * Xarm * reiacsum_a) - imvacsum_a;
        fx(7) = (vgrid_RE * reiacdif_a) + (vgrid_IM * imiacdif_a) + (vdcdif_a * idcdif_a) - pgrid_a;
        fx(8) = (vgrid_IM * reiacdif_a) - (vgrid_RE * imiacdif_a) - qgrid_a;
        fx(9) = (-vdcdif_a + vdcsum_a/2)*(idcdif_a/2 + idcsum_a) - revacdif_a*(reiacdif_a/2+reiacsum_a) - imvacdif_a*(imiacdif_a/2+imiacsum_a) + 0.5*revacsum_a*(reiacdif_a/2+reiacsum_a) + 0.5*imvacsum_a*(imiacdif_a/2+imiacsum_a) - pconu_a;
        fx(10) = (vdcdif_a + vdcsum_a/2)*(-idcdif_a/2 + idcsum_a) + revacdif_a*(-reiacdif_a/2+reiacsum_a) + imvacdif_a*(-imiacdif_a/2+imiacsum_a) + 0.5*revacsum_a*(-reiacdif_a/2+reiacsum_a) + 0.5*imvacsum_a*(-imiacdif_a/2+imiacsum_a) - pconl_a;    
        fx(11) = idcdif_a - idcdif_ref_a;
        fx(12) = imiacsum_a - imiacsum_ref_a;

        %computes f(x) for phase B
        fx(13) = vdcsum_b + (idcsum_b * R) - vhvdc;
        fx(14) = vdcdif_b - (idcdif_b * Rl);
        fx(15) = revacdif_b + (imiacdif_b * (Xl + Xarm/2)) - (reiacdif_b*Rl) - vgrid_RE_b;
        fx(16) = imvacdif_b - (reiacdif_b * (Xl + Xarm/2)) - (imiacdif_b*Rl) - vgrid_IM_b;
        fx(17) = (2 * Xarm * imiacsum_b) - revacsum_b;
        fx(18) = (-2 * Xarm * reiacsum_b) - imvacsum_b;
        fx(19) = (vgrid_RE * reiacdif_b) + (vgrid_IM * imiacdif_b) + (vdcdif_b * idcdif_b) - pgrid_b;
        fx(20) = (vgrid_IM * reiacdif_b) - (vgrid_RE * imiacdif_b) - qgrid_b;
        fx(21) = (-vdcdif_b + vdcsum_b/2)*(idcdif_b/2 + idcsum_b) - revacdif_b*(reiacdif_b/2+reiacsum_b) - imvacdif_b*(imiacdif_b/2+imiacsum_b) + 0.5*revacsum_b*(reiacdif_b/2+reiacsum_b) + 0.5*imvacsum_b*(imiacdif_b/2+imiacsum_b) - pconu_b;
        fx(22) = (vdcdif_b + vdcsum_b/2)*(-idcdif_b/2 + idcsum_b) + revacdif_b*(-reiacdif_b/2+reiacsum_b) + imvacdif_b*(-imiacdif_b/2+imiacsum_b) + 0.5*revacsum_b*(-reiacdif_b/2+reiacsum_b) + 0.5*imvacsum_b*(-imiacdif_b/2+imiacsum_b) - pconl_b;    
        fx(23) = idcdif_b - idcdif_ref_b;
        fx(24) = imiacsum_b - imiacsum_ref_b;

        %computes f(x) for phase C
        fx(25) = vdcsum_c + (idcsum_c * R) - vhvdc;
        fx(26) = vdcdif_c - (idcdif_c * Rl);
        fx(27) = revacdif_c + (imiacdif_c * (Xl + Xarm/2)) - (reiacdif_c*Rl) - vgrid_RE_c;
        fx(28) = imvacdif_c - (reiacdif_c * (Xl + Xarm/2)) - (imiacdif_c*Rl) - vgrid_IM_c;
        fx(29) = (2 * Xarm * imiacsum_c) - revacsum_c;
        fx(30) = (-2 * Xarm * reiacsum_c) - imvacsum_c;
        fx(31) = (vgrid_RE * reiacdif_c) + (vgrid_IM * imiacdif_c) + (vdcdif_c * idcdif_c) - pgrid_c;
        fx(32) = (vgrid_IM * reiacdif_c) - (vgrid_RE * imiacdif_c) - qgrid_c;
        fx(33) = (-vdcdif_c + vdcsum_c/2)*(idcdif_c/2 + idcsum_c) - revacdif_c*(reiacdif_c/2+reiacsum_c) - imvacdif_c*(imiacdif_c/2+imiacsum_c) + 0.5*revacsum_c*(reiacdif_c/2+reiacsum_c) + 0.5*imvacsum_c*(imiacdif_c/2+imiacsum_c) - pconu_c;
        fx(34) = (vdcdif_c + vdcsum_c/2)*(-idcdif_c/2 + idcsum_c) + revacdif_c*(-reiacdif_c/2+reiacsum_c) + imvacdif_c*(-imiacdif_c/2+imiacsum_c) + 0.5*revacsum_c*(-reiacdif_c/2+reiacsum_c) + 0.5*imvacsum_c*(-imiacdif_c/2+imiacsum_c) - pconl_c;    
        fx(35) = idcdif_c - idcdif_ref_c;
        fx(36) = imiacsum_c - imiacsum_ref_c;

    end
end