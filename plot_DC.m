%Plots DC components with a max voltage circuit
function plot_DC(Vdc, Idc, Vhvdc)
    %Max voltage
    angles_matrix = linspace(0,2*pi);
    radius = Vhvdc;
    Vmax_x = radius * cos(angles_matrix);
    Vmax_y = radius * sin(angles_matrix);

    %Actual plot
    figure
    hold on
    grid on
    axis equal
    quiver(0, 0, real(Vdc)/1000, imag(Vdc)/1000)
    quiver(0, 0, real(Idc)/1000, imag(Idc)/1000)
    plot(Vmax_x/1000, Vmax_y/1000, '--')
    legend('Vdc', 'Idc', 'Vmax')
    xlabel('Real [kV and kA]')
    ylabel('Imaginary [kV and kA]')
    title('Single Arm DC Results')
end