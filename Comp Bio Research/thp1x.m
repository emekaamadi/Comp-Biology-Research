function [err, timepoints, species_out, observables_out] = thp1x( timepoints, species_init, parameters, suppress_plot )
%THP1X Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'thp1x' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. THP1X returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = thp1x( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 10 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 21 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 6.022e23, 1e-4, 2.5e4, 10, 250, 6e4, 6e3, -11.61, 0.491, -10.85, -0.532, -1.24, -12, -1.38, -12, -1.11, -0.3, 1.17, 100, -1.9, -1.33 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 21  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 21].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 10  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 10].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,20,1000+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'N_A', 'vol', 'aml_init', 'et_ratio', 'a_pm', 'rho3_init', 'rho123_init', 'log_kp3', 'log_km3', 'log_kp123', 'log_km123', 'log_kf', 'log_kr', 'log_kb_tc', 'log_kb_aml', 'log_k_kill', 'log_k_tox1', 'log_k_act', 'nb', 'log_kd3', 'log_ki3' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-8,   ...
               'AbsTol',   0.0001,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
    if(length(timepoints) ~= size(species_out,1))
        exception = MException('ODE15sError:MissingOutput','Not all timepoints output\n');
        throw(exception);
    end
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), 15 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'DEAD_aml', 'CD3_tot', 'CD3_free', 'CD3_bound', 'CD123_tot', 'CD123_free', 'CD123_bound', 'bridge', 'tot_synapse', 'active_synapse', 'Tc_tot', 'Tc_free', 'Tc_active', 'AML_free', 'AML_tot' };

    % construct figure
    plot(timepoints,observables_out);
    title('thp1x observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%

% Define if function to allow nested if statements in user-defined functions
function [val] = if__fun (cond, valT, valF)
% IF__FUN Select between two possible return values depending on the boolean
% variable COND.
    if (cond)
        val = valT;
    else
        val = valF;
    end
end

% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,10);
    species_init(1) = ((params(2)*params(1))*params(5))*1e-12;
    species_init(2) = params(3)*params(7);
    species_init(3) = (params(3)*params(4))*params(6);
    species_init(4) = params(3);
    species_init(5) = params(3)*params(4);
    species_init(6) = 0;
    species_init(7) = 0;
    species_init(8) = 0;
    species_init(9) = 0;
    species_init(10) = 0;

end


% user-defined functions
% function r_synapse
function [val] = r_synapse(expressions, observables)
    val = if__fun( (((observables(11)>1)&&((observables(6)+observables(7))>=1))&&((observables(3)+observables(4))>=1)), (((expressions(7)*observables(2))*((((10^expressions(8))*observables(7))*observables(3))+(((10^expressions(10))*observables(6))*observables(4))))/((observables(11)*(observables(6)+observables(7)))*(observables(3)+observables(4)))), 0);
end

% function r_bridge
function [val] = r_bridge(expressions, observables)
    val = if__fun( (((observables(11)>1)&&((observables(6)+observables(7))>=1))&&((observables(3)+observables(4))>=1)), ((((expressions(19)*observables(12))*observables(2))*expressions(7))/((observables(11)*(observables(4)+observables(3)))*(observables(7)+observables(6)))), 0);
end

% function b3
function [val] = b3(expressions, observables)
    val = ((10^expressions(11))/((10^expressions(9))+(10^expressions(11))));
end

% function r_kill
function [val] = r_kill(expressions, observables)
    val = if__fun( (observables(8)>1), ((expressions(19)*observables(10))/observables(8)), 0);
end

% function r_diss
function [val] = r_diss(expressions, observables)
    val = if__fun( (observables(8)>1), ((expressions(19)*observables(9))/observables(8)), 0);
end

% function Q_Taf_3f
function [val] = Q_Taf_3f(expressions, observables)
    val = if__fun( (observables(3)>1), (observables(13)/observables(3)), 0);
end

% function f_k_tox
function [val] = f_k_tox(expressions, observables)
    val = if__fun( (observables(11)>1), (((observables(2)/observables(11))-expressions(19))/(observables(4)+observables(3))), 0);
end

% function f_k_kill
function [val] = f_k_kill(expressions, observables)
    val = if__fun( (observables(15)>1), (((observables(5)/observables(15))-expressions(19))/(observables(7)+observables(6))), 0);
end

% function rateLaw5
function [val] = rateLaw5(expressions, observables)
    val = ((10^expressions(12))*r_synapse(expressions,observables));
end

% function rateLaw6
function [val] = rateLaw6(expressions, observables)
    val = (((10^(expressions(12)+expressions(8)))*observables(14))*r_bridge(expressions,observables));
end

% function rateLaw7
function [val] = rateLaw7(expressions, observables)
    val = (((10^(expressions(12)+expressions(10)))*observables(14))*r_bridge(expressions,observables));
end

% function rateLaw9
function [val] = rateLaw9(expressions, observables)
    val = (((10^expressions(16))*(1-b3(expressions,observables)))*r_kill(expressions,observables));
end

% function rateLaw10
function [val] = rateLaw10(expressions, observables)
    val = (((10^expressions(16))*b3(expressions,observables))*r_kill(expressions,observables));
end

% function rateLaw11
function [val] = rateLaw11(expressions, observables)
    val = (((10^expressions(17))*observables(10))*f_k_kill(expressions,observables));
end

% function rateLaw12
function [val] = rateLaw12(expressions, observables)
    val = (((10^expressions(17))*observables(10))*f_k_kill(expressions,observables));
end




% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,31);
    expressions(1) = parameters(1);
    expressions(2) = parameters(2);
    expressions(3) = parameters(3);
    expressions(4) = parameters(4);
    expressions(5) = parameters(5);
    expressions(6) = parameters(6);
    expressions(7) = parameters(7);
    expressions(8) = parameters(8);
    expressions(9) = parameters(9);
    expressions(10) = parameters(10);
    expressions(11) = parameters(11);
    expressions(12) = parameters(12);
    expressions(13) = parameters(13);
    expressions(14) = parameters(14);
    expressions(15) = parameters(15);
    expressions(16) = parameters(16);
    expressions(17) = parameters(17);
    expressions(18) = parameters(18);
    expressions(19) = parameters(19);
    expressions(20) = parameters(20);
    expressions(21) = parameters(21);
    expressions(22) = (((expressions(2)*expressions(1))*expressions(5))*1e-12);
    expressions(23) = (expressions(3)*expressions(4));
    expressions(24) = (expressions(23)*expressions(6));
    expressions(25) = (expressions(3)*expressions(7));
    expressions(26) = (expressions(20)+(log(expressions(6))/log(10)));
    expressions(27) = (10^expressions(8));
    expressions(28) = (10^expressions(9));
    expressions(29) = (10^expressions(10));
    expressions(30) = (10^expressions(11));
    expressions(31) = (10^expressions(16));
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,15);
    observables(1) = species(10);
    observables(2) = species(3) +species(6) +species(9);
    observables(3) = species(3);
    observables(4) = species(6);
    observables(5) = species(2) +species(7) +species(9);
    observables(6) = species(2);
    observables(7) = species(7);
    observables(8) = species(9);
    observables(9) = species(8);
    observables(10) = species(8);
    observables(11) = species(5) +species(8);
    observables(12) = species(5);
    observables(13) = species(5);
    observables(14) = species(4);
    observables(15) = species(4) +species(8);

end


% Calculate ratelaws
function [ ratelaws ] = calcratelaws ( species, expressions, observables )

    ratelaws = zeros(1,15);
    ratelaws(1) = (10^expressions(8))*species(1)*species(3);
    ratelaws(2) = (10^expressions(10))*species(1)*species(2);
    ratelaws(3) = rateLaw5(expressions,observables)*species(4)*species(5);
    ratelaws(4) = rateLaw12(expressions,observables)*species(2);
    ratelaws(5) = (10^expressions(9))*species(6);
    ratelaws(6) = (10^expressions(11))*species(7);
    ratelaws(7) = rateLaw6(expressions,observables)*species(3)*species(7);
    ratelaws(8) = rateLaw7(expressions,observables)*species(6)*species(2);
    ratelaws(9) = (10^expressions(16))*species(8);
    ratelaws(10) = rateLaw11(expressions,observables)*species(7);
    ratelaws(11) = rateLaw9(expressions,observables)*species(9);
    ratelaws(12) = rateLaw10(expressions,observables)*species(9);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(10,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calcratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = -ratelaws(1) -ratelaws(2) +ratelaws(5) +ratelaws(6);
    Dspecies(2) = -ratelaws(2) -ratelaws(4) +ratelaws(6) -ratelaws(8);
    Dspecies(3) = -ratelaws(1) +ratelaws(5) -ratelaws(7) +ratelaws(11);
    Dspecies(4) = -ratelaws(3);
    Dspecies(5) = -ratelaws(3) +ratelaws(9);
    Dspecies(6) = ratelaws(1) -ratelaws(5) -ratelaws(8) +ratelaws(12);
    Dspecies(7) = ratelaws(2) -ratelaws(6) -ratelaws(7) -ratelaws(10);
    Dspecies(8) = ratelaws(3) -ratelaws(9);
    Dspecies(9) = ratelaws(7) +ratelaws(8) -ratelaws(11) -ratelaws(12);
    Dspecies(10) = ratelaws(9);

end


end
