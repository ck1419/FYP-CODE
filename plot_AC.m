%Simple AC Phasor Plot
function plot_AC(Vac, Iac)
    figure
    hold on
    grid on
    axis equal
    quiver(0, 0, real(Vac)/1000, imag(Vac)/1000)
    quiver(0, 0, real(Iac)/1000, imag(Iac)/1000)
    xlabel('Real [kV and kA]')
    ylabel('Imaginary [kV and kA]')
    legend('Vac', 'Iac')
    title('Single Arm AC Results')
end


