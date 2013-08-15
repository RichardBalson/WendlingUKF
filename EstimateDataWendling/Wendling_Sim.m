function Wendling_Sim(MVI,fs)

% script created by Richard Balson 21/02/2013

% description
% ~~~~~~~~~~~
% Simulate states and output of the extended neural mass model for
% estimation purposes. Input mean assumed to be 90.

% last edit
% ~~~~~~~~~

MVIT = MVI(:,1);
MVI(:,1:3) = MVI(:,2:4);
MVI(:,4) = MVIT;

% Variables required
% ~~~~~~~~~~~~~~~
addpath(genpath('../../UKFFinal')); % Specify files required for estimation
% stochastic -  Determines level of stochasticity in input
% number_of_sigma_input - Determines the standard deviation of stochastic input
% Parameter_index - Used to decide what parameter should be simulated
% Input_mean_variation - uSed to specifiy whether the mean of the input
% should be different from its described nominal value
SimulationSettings.name = 'VarianceSim'; % Specify name of file to save Wendling model output data to, or to load data from when simulation is not performed
SimulationSettings.fs = fs;
    SimulationSettings.slope_time =2.5; % Specifies the time over which the model gain should be altered
    SimulationSettings.number_of_sigma_input = 1; % Used to determine standard deviation of input if  1: 68.27% of realisations within physiolgical range, 2: 95.45, 3: 99.73 4: 99.994
    SimulationSettings.stochastic = 1; % Used to specifiy the stochastic adjustment on the input 1 is no adjustment. <1 downscalling, >1 upscaling
    SimulationSettings.Parameter_index = 8; % Choose parameters to be simulated: 1 = Seizure Parameter from Wendling 2002;
    SimulationSettings.Input_mean_variation = 0; % If 0 mean stays constant for simulation,
    %if 1 input mean is drawn from a uniform distribution limited by the physiological limits of the input
    % if 2 input mean is drawn from a Gaussian
    % distribution with a mean as per Wendling
    % 2002 and stnadard deviation that satisfies
    % number_of_sigma_input
    
% next edit
% ~~~~~~~~~

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

%%

% Model parameters(Variable parameters are to alter parameter values during
% simulation) Matrix specifies parameter values for a specified time period
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Units mV
% ~~~~~~~~~~~~

% Units Hz
% ~~~~~~~~

sampling_frequency = SimulationSettings.fs; % Sampling frequency for solutions and measurements

%%
% Stationery Parameter defintions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~
% Stationery Parameter defintions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Model parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Units seconds^(-1)
% ~~~~~~~~~~~~

a =100;             %Excitatory time constant
b =30;              %Slow inhibitory time constant original b=50
g =350;             %Fast inhibitory time constant g =500

tcon = [a b g]; % Specify reciprocal of the time constants for simulation


% Units NS
% ~~~~~~~~~

Con = 135; % Connectivity constant, used to specify connectivty between neuronal types

C= [Con; 0.8*Con; 0.25*Con; 0.25*Con; 0.3*Con; 0.1*Con; 0.8*Con]; % Connectivity Constants for all populations

% Input noise limits
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Units Hz
% ~~~~~~~~~~~~~

min_frequency = 30; % Minimum noise firing rate

max_frequency = 150; % Maximum noise input firing rate

frequency_limits = [min_frequency max_frequency];

% Solver Parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dt = 1/sampling_frequency; % Time step between solutions

% Units NS
% ~~~~~~~~~

simulation_time = dt*size(MVI,1);
time = 0:dt:simulation_time;  % Time frame over which solver should solve equations

solver_time = 0:dt:(simulation_time-dt); % Specifiy time parameters for solver loop, ends at simulation_time-dt so that the last solution can be determined

normalise = round((solver_time+dt)*sampling_frequency);% Normalise time for matrix use, used in for loops

no_solutions = sampling_frequency*simulation_time+1; % Specifies the number of solutions required


% Noise Parameters
% ~~~~~~~~~~~~~~~~~~~~
if (SimulationSettings.Input_mean_variation ==0) % Constant input mean
meanf = (min_frequency+max_frequency)/2; % Input mean
elseif (SimulationSettings.Input_mean_variation ==1) % Mean drawn from a uniform distribution
    meanf = min_frequency + (max_frequency-min_frequency)*rand(1);
else % Mean drawn from a Gaussin distribution
    meanf = (min_frequency+max_frequency)/2 + (max_frequency-min_frequency)/(2*number_of_sigma_input)*randn(1);
end
std_deviation = max((max_frequency-meanf)/(SimulationSettings.number_of_sigma_input), (meanf-min_frequency)/(SimulationSettings.number_of_sigma_input)); % Input standard deviation
 % Noise std_deviation


z = zeros(no_solutions,8);% Initialise Solutions for wendling paper model, where rows ...


% e = xxxx*randn(NSamples,1);
% Determine stochastic input to the model
% ~~~~~~~~~~~~~~~~~~~~~~~
gauss = randn(1,length(normalise))*std_deviation*SimulationSettings.stochastic;

%         gauss = randn(1)*std_deviation*stochastic; % Determine random fluctuations in input signal with the specified standard deviation

if size(MVI,2)>3
    normalised_gaussian_input(1,:) = MVI(:,4)' + gauss(1,:); % Determine value of input for all time steps
    normalised_gaussian_inputSDE(1,:) = MVI(:,4)' + gauss(1,:)*sqrt(dt)*sampling_frequency;
else
normalised_gaussian_input(1,:) = meanf + gauss(1,:); % Determine value of input for all time steps
normalised_gaussian_inputSDE(1,:) = meanf + gauss(1,:)*sqrt(dt)*sampling_frequency;
end

%         while (any(normalised_gaussian_input>max_frequency) || any(normalised_gaussian_input <min_frequency))
%
%            gauss(1,k) = randn(1)*std_deviation*stochastic;
%
%             normalised_gaussian_input(k) = meanf + gauss(1,k);
%
%         end

 % Determine value of input used for calculation purposes for all time steps

%%

% Equations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Solver
% ~~~~~~~~~~~~~~~~~~~~~~~
% Iteration to change parameters in function
for k = normalise
    %
    % 8th Order
    % Derivation of equation shown in pdf attached
    
    % Normalised gaussian input
    
    zs1 =  z(k,2)-C(4)*z(k,3)-z(k,4);
    zs2 = C(1)*z(k,1);
    zs3 = C(3)*z(k,1);
    zs4 = C(5)*z(k,1)-C(6)*z(k,3);
    [z(k+1,1) z(k+1,5)] = PSPkernel([z(k,1); z(k,5)],dt,MVI(k,1),a,sigmoid(zs1));
    [z(k+1,2) z(k+1,6)] = PSPkernel([z(k,2); z(k,6)],dt,MVI(k,1),a,normalised_gaussian_inputSDE(k) +C(2)*sigmoid(zs2));
    [z(k+1,3) z(k+1,7)] = PSPkernel([z(k,3); z(k,7)],dt,MVI(k,2),b,sigmoid(zs3));
    [z(k+1,4) z(k+1,8)] = PSPkernel([z(k,4); z(k,8)],dt,MVI(k,3),g,C(7)*sigmoid(zs4));
    
end% End of for loop

% Determine output of model
% ~~~~~~~~~~~~~~~~~~~~~~~~
output8 = z(:,2)-C(4)*z(:,3) -z(:,4);      % 8th order output, Normalised gaussian input

t = linspace(0,1/sampling_frequency*length(output8),length(output8));

fig_settings.label_fontsize = 10;            % point
fig_settings.tick_fontsize = 8;              % point
fig_settings.legend_fontsize = 10;
fig_settings.left_pos = 5;               % cms
fig_settings.bottom_pos = 5;             % cms
fig_settings.font_type = 'Arial';
fig_settings.dirname = 'Results';              % default directory for figure files
fig_settings.scale =0.5;

legendT = {'Forward Model Simulation'};
FMO=state_figure('Wendling Neural Mass Output','Obs',fig_settings,t,output8',legendT,[]);