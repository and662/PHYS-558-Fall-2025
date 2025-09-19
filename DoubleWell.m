% 558 Question 2, Andrew Mckenna 

clear; 
close all

V0 = 1.0; % depth of the well in Ry
xres = 400;
Rres = 50;
initial_well_buffer = 0.2; % how far appart are the well in units of a0

xvals = linspace( -12.0, 12.0, xres );
dx = xvals(2) - xvals(1);

well_inc = ceil(0.5*(1 + initial_well_buffer)/dx);
Rvals = zeros(Rres, 1);

Evals = zeros(Rres, 2);
Evals_TBA = zeros(Rres, 2);

% get the potential of a single well 
Vsingle_well = zeros(xres, 1);
for i = 1:xres % check if we are in the well:
    if -0.5 < xvals(i) && xvals(i) < 0.5
        Vsingle_well(i) = -1.0;
    end
end

% define the kinetic part of the Hamiltonian
KE = (-diag(ones(xres - 1, 1), -1) + 2*eye(xres, xres) - diag(ones(xres - 1, 1), 1))/(dx^2);

for i = 1:Rres

    Vdouble_well = zeros(xres, 1);

    Vleft_well = circshift(Vsingle_well, well_inc);
    Vright_well = circshift(Vsingle_well, -well_inc);
    Vdouble_well = Vleft_well + Vright_well;

    Hdouble_well = KE + 2*V0*diag( Vdouble_well );
    
    R = 2*well_inc*dx;
    Rvals(i) = R;
    
    well_inc = well_inc + 1; % push the wells one step appart 
    
    % diagonalize the Hamiltonian
    [eigenvectors, eigenvalues] = eig( Hdouble_well ); % by default eig does not sort eigenvalues:
    [evs, idx] = sort(diag( eigenvalues ));
    eps_double_well = diag(eigenvalues(idx, idx));
    psi_double_well = eigenvectors(:, idx);
    
    % which energy levels have negative energy?
    bound_states = find( eps_double_well <= 0 );
    Evals(i, :) = eps_double_well(1:2);
    

    % plot everything 
    figure(1);
    t = tiledlayout(2, 1);
    
    nexttile; % plot the wavefunction
    cla(gca, 'reset');
    hold on; grid on;

    bound_states = find( eps_double_well <= 0 );
    for j = 1:length( bound_states )
        
        psi_real = real( psi_double_well(:, bound_states(j)) );
        psi_imag = imag( psi_double_well(:, bound_states(j)) );
        
        % fix gauge:
        sgn = sign( psi_real(1) );
        psi_real = sgn*psi_real;
        psi_imag = sgn*psi_imag;

        plot( xvals, psi_real, 'DisplayName', sprintf( 'state %d', j) );
        plot( xvals, psi_imag, 'HandleVisibility', 'off' );

    end
    
    ylabel('Wavefunction');
    xlim([min(xvals), max(xvals)]);
    legend();
    

    nexttile; % Plot the potential 
    cla(gca, 'reset');
    hold on; grid on;

    plot( xvals, Vdouble_well );
    
    for j = 1:length(bound_states)
        yline( eps_double_well(bound_states(j)), 'r' );
    end

    xlabel('x [Borh radii a0]'); ylabel('E [Rydbergs Ry]'); ylim([-1.2, 0.2]);
    xlim([min(xvals), max(xvals)]);

    hold off; pause(0.2);

end

figure(2);
hold on; grid on;

plot( Rvals, Evals(:, 1) );
plot( Rvals, Evals(:, 2) );
xlabel('R [Borh radii a0]'); ylabel('E [Rydbergs Ry]');
