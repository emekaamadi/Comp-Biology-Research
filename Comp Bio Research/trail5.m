% this will run thp1x.m
 % function [err, timepoints, species_out, observables_out] = CD3_CD123( timepoints, species_init, parameters, suppress_plot )
 timepoints = linspace(0,20,100)';
 species_init = [];
 % param 3 is drug conc in picomolar
 parameters = [ 6.022e23, 1e-4, 2.5e4, 10, 2.5, 6e4, 6e3, -11.61, 0.491, -10.85, -0.532, 1.24, -12, -1.38, -12, 0, -0.3, 1.17, 100, -1.9, -1.33 ];
 suppress_plot = 1;
[err, timepoints, species_out, observables_out] = thp1x( timepoints, species_init, parameters, suppress_plot );
   % colon in left pos all rows in right you specified 6th column ( above )
   % 
    
    
   a_pm = logspace(-2,10,100); 
    z= zeros(2,100);
    for i = 1:length(a_pm)
        parameters(5) = a_pm(i);
        [err, timepoints, species_out, observables_out] = thp1x( timepoints, species_init, parameters, suppress_plot );
       z(2,i) = observables_out(100,1);
        z(1,i) = a_pm(i);
    end
       
%     parameters(3) =100;
%         [err, timepoints, species_out, observables_out] = CD3_CD123( timepoints, species_init, parameters, suppress_plot );
    plot(z(1,:),z(2,:));
     title('[Apm] x Dead AML','fontSize',14,'Interpreter','none');
    axis([0 a_pm(end) 0 inf]);
    xlabel('[a_pm]','fontSize',12,'Interpreter','none');
    ylabel('Amount of Dead AML','fontSize',12,'Interpreter','none');
    set(gca, 'XScale', 'log');
       
       
        
        
        
    
%     for a_pm = something this needs to be 10-6 to 10^6 on a log scale
%     paramters (3) = a_pm;
%     CD3 function
%     last_val= observables_out(100this is the 10 hour point,6);
%     % } want to save _pm and last_val at each loop
%     make an array that stores the concentration of a_pm and the lastval at the a_pm value 
%     final plot will have x axis on the log scale 
%     x = a_pm and the y will be last_ val 

