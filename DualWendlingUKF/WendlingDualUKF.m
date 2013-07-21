function [X Pxx X_Multi Pxx_Multi]= WendlingUKFDual
% script created by Richard Balson 21/02/2013

% description
% ~~~~~~~~~~~
% This script estimates synaptic gains of the extended neural mass
% model using the unscented Kalman filter, it also simualtes the states and
% output of the extended neural mass model. In this script states and
% parameters are realisations with a set mean and a defined standard
% deviation

% last edit
% ~~~~~~~~~

% Altered standard deviation for input, altered uncertainty in input, made
% limit on estimate of input

% next edit
% ~~~~~~~~~

% Work on esstimation of the input

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

% Clear workspace
% ~~~~~~~~~~~
close all
clear
clc

tic

addpath(genpath('../Wendling')); % Specify files required for estimation

system_dependent('setprecision',24); % Set the precision of accuracy in order to reduce the effect of rounding errors.

User_defined_parameters;

dt = 1/fs;

if filter_simulation
    band_coeff = filtercoeff(lowcutoff,highcutoff,fs);
end

% Dynamic variables
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Estimation Procedure Parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dx = Ds+Dp+Dk; % Number of dimensions of augmented state matrix, Note that estimated parameters and inputs are now considered to be 'slow states' in the estimation procedure

% Physiological range of Model gains
% ~~~~~~~~~~~~~~~~~

Max_A =10;
Min_A =0;
Max_B =40;
Min_B =0;
Max_G =40;
Min_G =0;

Max = [Max_A, Max_B, Max_G];
Min = [Min_A, Min_B, Min_G];

for q = 1:Simulation_number
    
    if Random_number_generator(1)
        rng(0);
        else
            rng(cputime);
    end
    if simulate
        SimulationSettings.fs = fs;
        Wendling_Simulation(SimulationSettings); % Simulate states and output of the extended neural mass model
    end
    load([SimulationSettings.name,'.mat'])% Load parameters from simulation
    % output8 contains the simulated
    % model output; normalised_gaussian
    % input contains the input to the
    % model; sampling_frequency
    % contains the sampling frequency
    % used for simulation; stochastic
    % decides whether the variance is
    % applied to the input; MVI
    % contains the model gains; tcon
    % the model time constants; C the
    % connectivity constants and
    % frequency_limits the limits for
    % the input frequency
    
    % Simulated signal data mV
    % ~~~~~~~~~~~~~~~~
    if filter_simulation
        output8 = filtfilt1(band_coeff,1,output8);
    end
    
    % Static Variables
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ~~~~~~~~~~~
    EstimationVariables;
    init =0;
    condition =1;
    
    while condition
        init=init+1
        conditionT = init<Reinitialise_parameters_attempts*(~Parameter_initialisation);
        Initialise_parameters;
        
        %%
        
        % UKF algorithm
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        count(q) = init;
        
        for p =1:Number_of_observations
            [Sigma(:,:,p) err] = Unscented_transform(Ds,Pxx(:,:,p),X(:,p),kappa);
            [SigmaP(:,:,p) err1] =Unscented_transform(Dp,PxxP(:,:,p),Xp(:,p),kappa);
            if (err ==1)
                break
            end
            
            %                 if (LimitEst ==1)
            % %                     for j = 1:Dp
            % %                         if (Sigma(Ds+Dk+j,:,p) <0) % Check if parameter breaches physiological range
            % %                             Sigma(Ds+Dk+j,:,p) =0;
            % %                         end
            % %                     end
            %                 end
            
            if Dp==3
                gain = Xp(:,p);
            elseif Dp ==2
                gain = [X(:,P); MVI(EstStart_Sample+p-1,3)];
            elseif Dp==1
                gain = [X(1,p); MVI(EstStart_Sample+p-1,2);MVI(EstStart_Sample+p-1,3)];
            else
                gain = [MVI(EstStart_Sample+p-1,1); MVI(EstStart_Sample+p-1,2); MVI(EstStart_Sample+p-1,3)];
            end
            
            if Dk ==1
                %                     if (LimitEst ==1)
                %                         if (Sigma(Ds+1,k,p) < 0)
                %                             Sigma(Ds+1,k,p) = 0;
                %                         end
                %                     end
                Input_var = Sigma(Ds+1,:,p);
            else
                Input_var = Input_mean;
            end
            
            [Xout(:,:,p) Yout(:,:,p)] = WNM(Sigma(:,:,p),dt,Input_var, gain, tcon,C);
            
            
            [ExpX(:,p) ExpY(:,p) Pxxn Pxyn Pyyn] = Expectation(Xout(:,:,p), Ds, Yout(:,:,p), 1,kappa);
            
            Pyyn = Pyyn +R;
            %
            Pxxn = Pxxn + Q;
            
            [X(:,p+1) Pxx(:,:,p+1)] = Kalman(ExpX(:,p), ExpY(:,p), Y(p), Pxxn, Pxyn, Pyyn);
            
            % Update parameters
            
%             [Xexp Yexp] = WNM(X(:,p),dt,Input_var, gain, tcon,C);

% Using estimted states as observations

            [XoutP(:,:,p) YoutP(:,:,p)] = WNM1(ExpX(:,p),dt,Input_var, SigmaP(:,:,p), tcon,C);
            
            
            PxxPn = PxxPn+ Qp;
            
            
            PyyPn = PyyPn+ R;
            
            K = PxyPn/PyyPn;
            Xp(:,p+1) = Xp(:,p) + K*(Y(p) - ExpY(1,p)); 
            
            if StepbyStepCheck==1 % perform a step by step check of results
                stepn = 1:p;
                stepC = 0:p;
                figure % Plot the expected value of Y and the observed Y
                plot(stepn,Y(1:p),'k');
                hold on
                plot(stepn, ExpY(:,1:p));
                hold on
                plot(stepC,X(2,1:p+1)-C(4)*X(3,1:p+1)-X(4,1:p+1),'r')
                PxxE(:,:,stepn(end)) = Pxxn;
                for zz = 1:8
                    figure
                    
                    for li = stepn
                        ernE(li) = ExpX(zz,li) - PxxE(zz,zz,li);
                        erpE(li) = ExpX(zz,li) + PxxE(zz,zz,li);
                    end
                    for li = stepC+1
                        ern(li) = X(zz,li) - Pxx(zz,zz,li);
                        erp(li) = X(zz,li) + Pxx(zz,zz,li);
                    end
                    plot(stepn, ExpX(zz,1:p),'b');
                    hold on
                    plot(stepn,z(EstStart_Sample:EstStart_Sample+p-1,zz),'k');
                    hold on
                    plot(stepC, X(zz,1:p+1),'r');
                    hold on
                    plot(stepC, ern(1:p+1),'g+');
                    hold on
                    plot(stepC,  erp(1:p+1),'g+');
                    hold on
                    plot(stepn,  erpE(1:p),'m+');
                    hold on
                    plot(stepn,  ernE(1:p),'m+');
                end
                pause
                close all
            end
            
        end
        if Dp==3
            conditionT1 = ((Xp(1,end) <0) || (Xp(2,end) <0) ||(Xp(3,end) <0));
        else
            conditionT1=1;
        end
        condition = conditionT && conditionT1;
    end
    toc
    % if (((p==length(Y)) && ((X(Ds+Dk+1,1,p) <Min(1)) || (X(Ds+Dk+1,1,p) >Max(1))) || ((X(Ds+Dk+2,1,p) <Min(2)) || (X(Ds+Dk+2,1,p) >Max(2)))|| ((X(Ds+Dk+3,1,p) <Min(3)) || (X(Ds+Dk+3,1,p) >Max(3)))))
    %     q=q-1;
    %     continue;
    % end
    
    if (Simulation_number ~=1)
        if q ==1
            X_Multi = zeros(floor(size(X,2)/Decimate)+1,Dp+Dk,Simulation_number);
            Pxx_Multi = zeros(floor(size(X,2)/Decimate)+1,Dp+Dk,Simulation_number);
        end
        X_Multi(:,:,q) = X(Ds+1:Ds+Dp+Dk,1:500:end)';
        for k =1:Dk+Dp
            Pxx_Multi(:,k,q) = squeeze(Pxx(k,k,1:500:end));
        end
    else
        X_Multi =0;
        Pxx_Multi =0;
    end
end % End Simulation_nuumber loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Post processing and image handling%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % Determne how to plot estimation output
%     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%

if (tstart+zoomtime) >= (SimulationSettings.simulation_time-EstStart)
    zoomtime = SimulationSettings.simulation_time-EstStart-tstart-0.5;
end

Error_calculation;

Generate_figures;

if fig_save
    
    Figure_handling;
    
end
if Simulation_number>1
    Error_calculation_multi;
    Generate_figures_multi;
    if fig_save
        
        Figure_handling_multi;
        
    end
end
%
% load Handel
% sound(y,Fs)
%
% toc
%
%
