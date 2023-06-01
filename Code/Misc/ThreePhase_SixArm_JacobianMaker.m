clear
clc

syms in R Rl Xl Xarm vhvdc vgrid_RE vgrid_IM pgrid qgrid pconu pconl sgrid idcdif_ref imiacsum_ref
syms vdcsum_a vdcdif_a idcdif_a idcsum_a revacsum_a imvacsum_a revacdif_a imvacdif_a reiacsum_a imiacsum_a reiacdif_a imiacdif_a
syms vdcsum_b vdcdif_b idcdif_b idcsum_b revacsum_b imvacsum_b revacdif_b imvacdif_b reiacsum_b imiacsum_b reiacdif_b imiacdif_b
syms vdcsum_c vdcdif_c idcdif_c idcsum_c revacsum_c imvacsum_c revacdif_c imvacdif_c reiacsum_c imiacsum_c reiacdif_c imiacdif_c

state_variables = [vdcsum_a vdcdif_a idcdif_a idcsum_a revacsum_a imvacsum_a revacdif_a imvacdif_a reiacsum_a imiacsum_a reiacdif_a imiacdif_a vdcdif_b idcdif_b idcsum_b revacsum_b imvacsum_b revacdif_b imvacdif_b reiacsum_b imiacsum_b reiacdif_b imiacdif_b vdcsum_c vdcdif_c idcdif_c idcsum_c revacsum_c imvacsum_c revacdif_c imvacdif_c reiacsum_c imiacsum_c reiacdif_c imiacdif_c];


fx(1) = vdcsum_a + (idcsum_a * R) - vhvdc;
fx(2) = vdcdif_a - (idcdif_a * Rl);
fx(3) = revacdif_a + (imiacdif_a * (Xl + Xarm/2)) - (reiacdif_a*Rl) - vgrid_RE;
fx(4) = imvacdif_a - (reiacdif_a * (Xl + Xarm/2)) - (imiacdif_a*Rl) - vgrid_IM;
fx(5) = (2 * Xarm * imiacsum_a) - revacsum_a;
fx(6) = (-2 * Xarm * reiacsum_a) - imvacsum_a;
fx(7) = (vgrid_RE * reiacdif_a) + (vgrid_IM * imiacdif_a) + (vdcdif_a * idcdif_a) - pgrid;
fx(8) = (vgrid_IM * reiacdif_a) - (vgrid_RE * imiacdif_a) - qgrid;
fx(9) = (-vdcdif_a + vdcsum_a/2)*(idcdif_a/2 + idcsum_a) - revacdif_a*(reiacdif_a/2+reiacsum_a) - imvacdif_a*(imiacdif_a/2+imiacsum_a) + 0.5*revacsum_a*(reiacdif_a/2+reiacsum_a) + 0.5*imvacsum_a*(imiacdif_a/2+imiacsum_a) - pconu;
fx(10) = (vdcdif_a + vdcsum_a/2)*(-idcdif_a/2 + idcsum_a) + revacdif_a*(-reiacdif_a/2+reiacsum_a) + imvacdif_a*(-imiacdif_a/2+imiacsum_a) + 0.5*revacsum_a*(-reiacdif_a/2+reiacsum_a) + 0.5*imvacsum_a*(-imiacdif_a/2+imiacsum_a) - pconl;    
fx(11) = idcdif_a - idcdif_ref;
fx(12) = imiacsum_a - imiacsum_ref;

%computes f(x) for phase B
fx(13) = vdcsum_b + (idcsum_b * R) - vhvdc;
fx(14) = vdcdif_b - (idcdif_b * Rl);
fx(15) = revacdif_b + (imiacdif_b * (Xl + Xarm/2)) - (reiacdif_b*Rl) - vgrid_RE;
fx(16) = imvacdif_b - (reiacdif_b * (Xl + Xarm/2)) - (imiacdif_b*Rl) - vgrid_IM;
fx(17) = (2 * Xarm * imiacsum_b) - revacsum_b;
fx(18) = (-2 * Xarm * reiacsum_b) - imvacsum_b;
fx(19) = (vgrid_RE * reiacdif_b) + (vgrid_IM * imiacdif_b) + (vdcdif_b * idcdif_b) - pgrid;
fx(20) = (vgrid_IM * reiacdif_b) - (vgrid_RE * imiacdif_b) - qgrid;
fx(21) = (-vdcdif_b + vdcsum_b/2)*(idcdif_b/2 + idcsum_b) - revacdif_b*(reiacdif_b/2+reiacsum_b) - imvacdif_b*(imiacdif_b/2+imiacsum_b) + 0.5*revacsum_b*(reiacdif_b/2+reiacsum_b) + 0.5*imvacsum_b*(imiacdif_b/2+imiacsum_b) - pconu;
fx(22) = (vdcdif_b + vdcsum_b/2)*(-idcdif_b/2 + idcsum_b) + revacdif_b*(-reiacdif_b/2+reiacsum_b) + imvacdif_b*(-imiacdif_b/2+imiacsum_b) + 0.5*revacsum_b*(-reiacdif_b/2+reiacsum_b) + 0.5*imvacsum_b*(-imiacdif_b/2+imiacsum_b) - pconl;    
fx(23) = idcdif_b - idcdif_ref;
fx(24) = imiacsum_b - imiacsum_ref;

%computes f(x) for phase C
fx(25) = vdcsum_c + (idcsum_c * R) - vhvdc;
fx(26) = vdcdif_c - (idcdif_c * Rl);
fx(27) = revacdif_c + (imiacdif_c * (Xl + Xarm/2)) - (reiacdif_c*Rl) - vgrid_RE;
fx(28) = imvacdif_c - (reiacdif_c * (Xl + Xarm/2)) - (imiacdif_c*Rl) - vgrid_IM;
fx(29) = (2 * Xarm * imiacsum_c) - revacsum_c;
fx(30) = (-2 * Xarm * reiacsum_c) - imvacsum_c;
fx(31) = (vgrid_RE * reiacdif_c) + (vgrid_IM * imiacdif_c) + (vdcdif_c * idcdif_c) - pgrid;
fx(32) = (vgrid_IM * reiacdif_c) - (vgrid_RE * imiacdif_c) - qgrid;
fx(33) = (-vdcdif_c + vdcsum_c/2)*(idcdif_c/2 + idcsum_c) - revacdif_c*(reiacdif_c/2+reiacsum_c) - imvacdif_c*(imiacdif_c/2+imiacsum_c) + 0.5*revacsum_c*(reiacdif_c/2+reiacsum_c) + 0.5*imvacsum_c*(imiacdif_c/2+imiacsum_c) - pconu;
fx(34) = (vdcdif_c + vdcsum_c/2)*(-idcdif_c/2 + idcsum_c) + revacdif_c*(-reiacdif_c/2+reiacsum_c) + imvacdif_c*(-imiacdif_c/2+imiacsum_c) + 0.5*revacsum_c*(-reiacdif_c/2+reiacsum_c) + 0.5*imvacsum_c*(-imiacdif_c/2+imiacsum_c) - pconl;    
fx(35) = idcdif_c - idcdif_ref;
fx(36) = imiacsum_c - imiacsum_ref;

out = jacobian(fx, state_variables);