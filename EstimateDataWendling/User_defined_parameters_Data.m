% script created by Richard Balson 12/06/13
% This script specifies all the parameters that need to be defined by the
% user in roder to estimate simulated data from the Wendling model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% General parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Simulation_number =20;

filter_simulation =1; % Specify whether or not to filter simulated data

Normalise_data =1; % Specify whether or not data should be normalised.

Max_sim_voltage = 50;%20

Forward_model =0; % Specify whether the forward model should be simulated with the results found from the estimtion procedure.

if filter_simulation
    highcutoff = 2.5; % Specify highcutoff frequency for filter
    
    lowcutoff = 40; % Specify low cutoff frequency for filter
    % Data will have frequency content between highcutoff and lowcutoff
end

if Simulation_number>1
    Decimate = 500; % Specify the distance between corresponding samples for the output matrix when multiple simulation are performed
else Decimate =1;
end

StepbyStepCheck =0; % Check how prediction and correction steps are working on a plot

LimitEst = 0; % If set to 1 estimation results are limited by physiological parameters within the estimation procedure, that is sigma points cannot be aboe the physiological parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Estimation decisions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanCheckTStart = 15; % For multiple simulations determine where to start checking the estimation values

EstStart = 0; % Specify the duration after simulation start when estimation should start. This allows removal of all transients.

% Estimation states
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ds = 8; % Number of differential equations describing model, also the number of fast states to be estiamted

Dp = 3; % Number of parameters to be estimated, also refered to as slow states

Dk =1; %If set to 1 the mean of the stochastic input will be estimated % Note that if Input_mean_variation is not zero than Dk should be set to one to allow tracking of the input mean

Dy =1; % Number of observable outputs from the simulation

% Intialisation of parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

number_of_sigma = 1; % Number of standard deviations from mean. 4 accounts for 99.73 percent of points.

Reinitialise_parameters_attempts = 1; % Specify number of attempts for parameter reinitialisation if results are not physiologucally possible

% Uncertainty parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

num_input_sigma=1;

kappa =0; % Varibale used to define the relative contribution of the mean on the propogation of states, and adjustment of the variance of the sigma points drawn from the Gaussian distribution

% Variable_state_uncertainty = 0;%1e-3; % 1e-3 Uncertianty due to stochastic input

State_uncertainty_adjustment = [1 2 3 40 1 2 3 40]; %[1 1 1 22 1 1 1 22];%[1 1.5 5 20 1 1.5 5 20];%MAy be too far[1 10 20 60 1 10 20 60];% Exponential decrease in uncertainty % All ones good for slow but steady convergence

Base_parameter_uncertainty = 1e-2;%1e-2;%1e-12;%1e-3; % Inherent parameter uncertainty due to model error

Variable_parameter_uncertainty = 0;%1e-3;  % Uncertianty due parameters varying in time

Base_input_uncertainty = 1e-6;%1e-12;%1e-3; % Inherent parameter uncertainty due to model error

Variable_input_uncertainty =0;%1e-3; % Uncertianty due varying input mean, Set to zero if the input mean is not varying

Observation_uncertainty = 1e-1;%1e-12; %1 Specify the uncertainty in observations

uncertainty_adjustment = 1; % Adjuster for model uncertainty

std_adjustment_parameters =1; %2% Variance adjuster for parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Image handling parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_save =1; % Save figures as .fig for future use

printpdf =1;

Image_handling_model_output=[0;0];

plot_uncertainty =1; % Plot covariance of all states

Image_handling_states = [0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0]; % Here a decision is made whether to plot specific states, 
                                                                            % if the value is one the relevant figure is plotted, otherwise it is not.
                                                                            % The columns indicate the state to plot and the rows indicate whether the whole simulation or a zoomed in ploted should be plotted.
                                                                            % Column one corresponds with state 1 and so forth.
Image_handling_inputs = [0 0 0 0;0 0 0 0];  % Here a decision is made whether to plot specific states, 
                                                                            % if the value is one the relevant figure is plotted, otherwise it is not.
                                                                            % The columns indicate the state to plot and the rows indicate whether the whole simulation or a zoomed in ploted should be plotted.
                                                                            % Here column 1-4 are Vp,Ve,Vsi and Vfi respectively.

Image_handling_firing_rates = [0 0 1 0]; % Here a decision is made whether to plot specific states, 
                                                                            % if the value is one the relevant figure is plotted, otherwise it is not.
                                                                            % The columns indicate the firing rate to plot. this is a three image plot where the input potential population firing rate and output potential are plotted.
                                                                            % Here column 1-4 are Vp,Ve,Vsi and Vfi respectively.

Image_handling_multi = [1 1 1;0 0 0];%  % Here a decision is made whether to plot specific states, 
                                                                            % if the value is one the relevant figure is plotted, otherwise it is not.
                                                                            % The columns indicate the figures to be plotted.
                                                                            % Here column 1-4 are for all the model states, all the model states inputs, all the model parameters.
plot_uncertaintyMulti =0;
% Zoom parameters (seconds)
% ~~~~~~~~~~~~~~~~~

tstart =0; % Starting time for zoom

zoomtime = 10; % Duration of zoom



