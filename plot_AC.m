%Simple AC Phasor Plot
function plot_AC(Vac, Iac, msg)
    figure
    hold on
    grid on
    axis equal
    quiver(0, 0, real(Vac)/1e3, imag(Vac)/1e3)
    quiver(0, 0, real(Iac), imag(Iac))
    annotation('textbox', [.131 .131 .795 .795],'String',msg,'FitBoxToText','on');
    xlabel('Real [kV and A]')
    ylabel('Imaginary [kV and A]')
    legend('Vac (kV)', 'Iac (A)')
    title('Single Arm AC Results')
end


