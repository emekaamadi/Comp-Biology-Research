% this will run CD_CD123.m
 % function [err, timepoints, species_out, observables_out] = CD3_CD123( timepoints, species_init, parameters, suppress_plot )
 timepoints = linspace(0,20,100)';
 species_init = [];
%  param 3 is drug conc in picomolar
  parameters = [ 6.022e23, 1e-4, 100, 2.5e5, 25000, 60000, 6000, 4.09e4, 3.1032, 2.32e5, 0.29376, .00023645 ];
  % 7.3e-6 is the original kds parameter 12 
 suppress_plot = 1;
[err, timepoints, species_out, observables_out] = CD3_CD123( timepoints, species_init, parameters, suppress_plot );
   % colon in left pos all rows in right you specified 6th column ( above )

actual = zeros(100,2);
ksteps = 100;
z= zeros(ksteps,5);
kon = logspace(8,0,ksteps);
koff3 = 3.1032;
koff123 = 0.29376;
for i = 1:length(kon)
%     z(i,1) = koff(i);
    z(i,2) = kon(i);
    kcd3 = koff3 / kon(i);
    kcd123 = koff123 /kon(i);
    z(i,3) = kcd3;
    z(i,4) = kcd123;
    actual(i,1) = kcd3 / 3600;
    actual(i,2) = kcd123 / 3600;
end
% kd = z(:,3);
kcd3 = z(:,3);
kcd123 = z(:,4);

a_pm = logspace(-1,10,ksteps);
doseresponse = zeros(100,2);
ii = 1;
arr = zeros(100,5);
for l = 1:length(kcd3)
    fprintf("%d\n",l);
    parameters(8) = kon(l);
%     parameters(9) =koff(l);
    for j = 1:length(kcd123)
        parameters(10) = kon(j);
%         parameters(11) =koff(j);
        arr(ii,2) = kcd123(j);
        for i = 1:length(a_pm)
            parameters(3) = a_pm(i);
            [err, timepoints, species_out, observables_out] = CD3_CD123( timepoints, species_init, parameters, suppress_plot );
            doseresponse(i,1) = observables_out(100,6);
            doseresponse(i,2) = a_pm(i);
        end
        
        [width,maxa_pm,height] = get_shapedata(doseresponse);
        arr(ii,1) = kcd3(l);
        arr(ii,3) = width;
        arr(ii,4) = maxa_pm;
        arr(ii,5) = height;
        hillslope = arr(ii,5)/(log(arr(ii,4)));
        ii = ii+1;
        
        
    end
end

%%
%     filename = 'arrfull.mat';
%    load(filename)

%     x = arr(:,1);             % The range of x values.
%     y = arr(:,2);             % The range of y values.
%     [X,Y] = meshgrid (x,y); % This generates the actual grid of x and y values.
%     Z = arr(:,4);            % The the width in this case.
%     figure(1);              % Generating a new window to plot in.
%     % imagesc(x,y,Z)
%     contourf(X,Y,Z)
%     colormap hot
%     colorbar
%     xlabel('kcd3');
%     ylabel('kcd123');
% %     Xtick = 10.^(-12:-4); 
% %   

Z = zeros(100);  
for i =  1:10000
        Z((ceil(i/100)),(mod(i-1,100) +1)) = arr(i,3);
end
 x =[-12 -4];
 %kcd3
 y =[-13 -5];
 %kcd123
imagesc(x,y,Z)
%surf(x,y,Z);
colormap cool;
colorbar;
h = colorbar;
ylabel(h, 'Max Dead AML')
xlabel('KCD3 (logscale)');
ylabel('KCD123 (logscale)'); 
title('Heatmap shows the width of the apm x max dead aml cuve as KCD3 and KCD123 change');


% x = kcd3;            
% y = kcd123;  
% [X,Y] = meshgrid (x,y);
% Z = zeros(100);  
% for i =  1:10000
%         Z(X,Y) = arr(i,3);
% end
% contourf(X,Y,Z);


   
    
   
%     
%     
%     % 2d heatmap, probably better option
%    % https://www.mathworks.com/matlabcentral/answers/59463-plot-heatmap-with-3-variables
%          
%            
           
                    
            
            
            
            
                
                
                
            
       
% 
%     plot(z(1,:),z(2,:));
%      title('CD3_CD123 observables','fontSize',14,'Interpreter','none');
%     axis([0 a_pm(end) 0 inf]);
%     xlabel('a_pm','fontSize',12,'Interpreter','none');
%     ylabel('number or concentration','fontSize',12,'Interpreter','none');
%     set(gca, 'XScale', 'log');
%     
  
             
    
%     for a_pm = something this needs to be 10-6 to 10^6 on a log scale
%     paramters (3) = a_pm;
%     CD3 function
%     last_val= observables_out(100this is the 10 hour point,6);
%     % } want to save _pm and last_val at each loop
%     make an array that stores the concentration of a_pm and the lastval at the a_pm value 
%     final plot will have x axis on the log scale 
%     x = a_pm and the y will be last_ val 

% Find the concentration of A-Pm at the max and find the mid apm on the way up and the way down
% How do the 4 paramters effect the width 1/2 (sd) and the center of the curve with K3 and 123 on and off rates can do the  same thing to find the mean or max value 
% EVery run will give you a plot if you think about it and that will rep a width and a value for KD3 and KD123 those will go into a table that will then be used to make a 3d a graph
% 10^-12 to 10^-4 this is for the for loop
% for the kd value kd = off/ on leave the on and the the off probably between .001 to 100 on logspace
%     if you do this than wont have change in the for loop 
% kd for both and max on one plot as well as kd for both and width
% need to make sure function call is good and working 
% then look at picture on phone and implement that ii stuff 
% look and see if it is flat at the top which max value is actually giving
% you we want the one closest to the center
% for the width why is one in more control than the other vertically
% changes alot but not horizontally. Whats happening in the model think
% about that and change something to test that 
% Subtle diffrences see if thats a real thing or just round off error based off what i calc
% what other paramters beside kd may change the shape. 
% look into a contour plot and surface plot  
% one for loop that does what that double for loop does 


% The reason for changes horzontally may be a result of there being more kcd123 receptors and a greater t cell count as well 
% so changes in kd would have a lower impact as because there are so many more so that makes one axis less dependent then the other.
% That leads to the point that the number of tcells or cancer cells and their receptors will have an impact on shape. 
% This graph makes sense because when one is one extreme and the other is the other extreme the width will be the largest because
% it solves the problem of there being so much drug that it binds to kcd3 and kcd123 receptors seperately. 
% In this way one will have a weak affinity so the width will be large because even with an increase in drug conc one will
% have a weak affinity for the drug and not allow them to be seperatly saturated keeping the dead aml high. 
% 
% test max by seeing if it returns the middle max value at plateau test by
% Anwser = the max returns the index of the first time the max is seen
% making your own array and seeing what index it is returning 
%%% This week
%look at height which is the max dead aml but not the apm at the max dead
%aml which is the max we are looking at
% hill slope approx using the mid value remember the bottom is on a log
% scale. delta effect or height/ delta log apm
%keep one of the receptors constant and change the other and then graph the
%change in the other by the ratio of CD3/CD123 remember to use the bottom
%line 20 and 21
% so look at some spots compared to some other how spots in the middle of
% the picture would have the same max but differ in width description of
% how the paramters of the model control the shape of the dose response
% curve 
% calculate max the way you thought by adding the first and last time of
% the max and dividing by two and making that your max index
% for the slope you can take the point after and before the middle on the
% way up the curve that difference divded by the apm difference will give
% you the slope
% for figures .fig and .png files title and log scale units to the figures 


