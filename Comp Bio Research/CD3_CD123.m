function [err, timepoints, species_out, observables_out] = CD3_CD123( timepoints, species_init, parameters, suppress_plot )
%CD3_CD123 Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'CD3_CD123' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. CD3_CD123 returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = CD3_CD123( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 9 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 12 model parameters.
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
   parameters = [ 6.022e23, 1e-4, 100, 2.5e5, 25000, 60000, 6000, 4.09e4, 3.1032, 2.32e5, 0.29376, 7.2e-6 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 12  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 12].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 9  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 9].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,60,10+1)';
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
param_labels = { 'NA', 'V', 'a_pm', 'T0', 'C0', 'rho3', 'rho123', 'k_on3_Ms', 'k_off3', 'k_on123_Ms', 'k_off123', 'kds' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-6,   ...
               'AbsTol',   1e-6,   ...
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
observables_out = zeros( length(timepoints), 8 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'Syn', 'bound_CD3', 'tot_CD3', 'bound_CD123', 'tot_CD123', 'AML_dead', 'C', 'CD123' };

    % construct figure
    plot(timepoints,observables_out);
    title('CD3_CD123 observables','fontSize',14,'Interpreter','none');
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

    species_init = zeros(1,9);
    species_init(1) = params(4);
    species_init(2) = params(5);
    species_init(3) = params(6)*params(4);
    species_init(4) = params(7)*params(5);
    species_init(5) = ((params(2)*params(1))*params(3))*1e-12;
    species_init(6) = 0;
    species_init(7) = 0;
    species_init(8) = 0;
    species_init(9) = 0;

end


% user-defined functions
% function frac_bound_cd3
function [val] = frac_bound_cd3(expressions, observables)
    val = (observables(2)/observables(3));
end

% function frac_bound_cd123
function [val] = frac_bound_cd123(expressions, observables)
    val = (observables(4)/observables(5));
end

% function dead_AML
function [val] = dead_AML(expressions, observables)
    val = ((observables(6)/(observables(7)+observables(6)))*100);
end

% function k_kill
function [val] = k_kill(expressions, observables)
    val = ((expressions(17)*observables(1))/observables(7));
end

% function k_lose123
function [val] = k_lose123(expressions, observables)
    val = (((expressions(17)*observables(1))*expressions(7))/observables(8));
end




% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,17);
    expressions(1) = parameters(1);
    expressions(2) = parameters(2);
    expressions(3) = parameters(3);
    expressions(4) = parameters(4);
    expressions(5) = parameters(5);
    expressions(6) = parameters(6);
    expressions(7) = parameters(7);
    expressions(8) = (((expressions(2)*expressions(1))*expressions(3))*1e-12);
    expressions(9) = parameters(8);
    expressions(10) = (((expressions(9)/expressions(1))/expressions(2))*3600);
    expressions(11) = parameters(9);
    expressions(12) = parameters(10);
    expressions(13) = (((expressions(12)/expressions(1))/expressions(2))*3600);
    expressions(14) = parameters(11);
    expressions(15) = (expressions(6)*expressions(4));
    expressions(16) = (expressions(7)*expressions(5));
    expressions(17) = parameters(12);
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,8);
    observables(1) = species(9);
    observables(2) = species(6);
    observables(3) = species(3) +species(6) +species(9);
    observables(4) = species(7);
    observables(5) = species(4) +species(7) +species(9);
    observables(6) = species(8);
    observables(7) = species(2);
    observables(8) = species(4) +species(7) +species(9);

end


% Calculate ratelaws
function [ ratelaws ] = calc_ratelaws ( species, expressions, observables )

    ratelaws = zeros(1,8);
    ratelaws(1) = (((expressions(9)/expressions(1))/expressions(2))*3600)*species(3)*species(5);
    ratelaws(2) = (((expressions(12)/expressions(1))/expressions(2))*3600)*species(4)*species(5);
    ratelaws(3) = k_kill(expressions,observables)*species(2);
    ratelaws(4) = k_lose123(expressions,observables)*species(4);
    ratelaws(5) = (((expressions(9)/expressions(1))/expressions(2))*3600)*species(3)*species(7);
    ratelaws(6) = expressions(11)*species(6);
    ratelaws(7) = (((expressions(12)/expressions(1))/expressions(2))*3600)*species(4)*species(6);
    ratelaws(8) = expressions(14)*species(7);
    ratelaws(9) = k_lose123(expressions,observables)*species(7);
    ratelaws(10) = expressions(11)*species(9);
    ratelaws(11) = expressions(14)*species(9);
    ratelaws(12) = expressions(17)*species(9);
    ratelaws(13) = k_lose123(expressions,observables)*species(9);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(9,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calc_ratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = 0.0;
    Dspecies(2) = -ratelaws(3);
    Dspecies(3) = -ratelaws(1) -ratelaws(5) +ratelaws(6) +ratelaws(10);
    Dspecies(4) = -ratelaws(2) -ratelaws(4) -ratelaws(7) +ratelaws(8) +ratelaws(11);
    Dspecies(5) = -ratelaws(1) -ratelaws(2) +ratelaws(6) +ratelaws(8);
    Dspecies(6) = ratelaws(1) -ratelaws(6) -ratelaws(7) +ratelaws(11) +ratelaws(12);
    Dspecies(7) = ratelaws(2) -ratelaws(5) -ratelaws(8) -ratelaws(9) +ratelaws(10);
    Dspecies(8) = ratelaws(3);
    Dspecies(9) = ratelaws(5) +ratelaws(7) -ratelaws(10) -ratelaws(11) -ratelaws(12) -ratelaws(13);

end


end
