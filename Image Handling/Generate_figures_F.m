% script created by Richard Balson 20/01/2013

% description
% ~~~~~~~~~~~
% this script generates figures for UKF estimation results for the Document
% for the Boon group

% last edit
% ~~~~~~~~~


% next edit
% ~~~~~~~~~

% At line 184 editing multi figure plots

% Beginning of script
% ~~~~~~~~~~~~~~~~~~~~~

% Generic Figure handling parameters

fig_settings.label_fontsize = 10;            % point
fig_settings.tick_fontsize = 8;              % point
fig_settings.legend_fontsize = 10;
fig_settings.left_pos = 5;               % cms
fig_settings.bottom_pos = 5;             % cms
fig_settings.font_type = 'Arial';
fig_settings.dirname = 'Results';              % default directory for figure files
fig_settings.scale =0.5;

% Plot of output
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maxlimit = round(sampling_frequency/10):length(check);

% % Plot of state 1
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Error_multiplier = 1;
scale = 0.5;

t = linspace(0,dt*length(check),length(check));
tz = linspace(tstart,tstart+zoom,zoom*sampling_frequency+1);

index = {1:length(check) tstart*sampling_frequency+1:(tstart+zoom)*sampling_frequency+1 ...
    EstStart_Sample:length(output8) round((EstStart+tstart)*sampling_frequency)+1:round((EstStart+tstart+zoom)*sampling_frequency)+1};

if (Image_handling_model_output(1,1))
    
    legendT = {'Noisy Obs.','Estimated Obs.','Obs.','Location'};
    NMO=state_figure('Wendling Neural Mass Output','Obs',fig_settings,t,[Y(index{1})'; X(2,index{1})-C(4)*X(3,index{1})-X(4,index{1}); check(index{1})'],legendT,[]);
    
end

if (Image_handling_model_output(2,1))
    
    legendT = {'Noisy Obs.','Estimated Obs.','Obs.'};
    NMO=state_figure('Wendling Neural Mass Output Zoomed In','Obs',fig_settings,tz,[Y(index{2})'; X(2,index{2})-C(4)*X(3,index{2})-X(4,index{2}); check(index{2})'],legendT,[]);
    
end
legendT = {'Sim. v_{p0}','Est. v_{p0}','Std. Dev.';...
    'Sim. v_{p1}','Est. v_{p1}','Std. Dev.';...
    'Sim. v_{p2}','Est. v_{p2}','Std. Dev.';...
    'Sim. v_{p3}','Est. v_{p3}','Std. Dev.';...
    'Sim. z_{p0}','Est. z_{p0}','Std. Dev.';...
    'Sim. z_{p1}','Est. z_{p1}','Std. Dev.';...
    'Sim. z_{p2}','Est. z_{p2}','Std. Dev.';...
    'Sim. z_{p3}','Est. z_{p3}','Std. Dev.'};

fig_name = {'Neural Mass State1 (v_{p0})';...
    'Neural Mass State2 (v_{p1})';...
    'Neural Mass State3 (v_{p2})';...
    'Neural Mass State4 (v_{p3})';...
    'Neural Mass State5 (z_{p0})';...
    'Neural Mass State6 (z_{p1})';...
    'Neural Mass State7 (z_{p2})';...
    'Neural Mass State8 (z_{p3})'};


for k = 1:Ds
    
    if (Image_handling_states(1,k))
        % Plot of state Not Zoomed in
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        erfn = X(k,index{1})-squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
        erfp = X(k,index{1})+squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
        NMS(k)=state_figure(fig_name{k},'State',fig_settings,t,[z(index{3},k)';X(k,index{1})],legendT(k,:),[erfn;erfp]);
        
    end
end

for k = 1:Ds
    
    if (Image_handling_states(2,k))
        % Plot of state Not Zoomed in
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        erfn = X(k,index{2})-squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
        erfp = X(k,index{2})+squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
        NMSZ(k)=state_figure(fig_name{k},'State',fig_settings,tz,[z(index{4},k)';X(k,index{2})],legendT(k,:),[erfn;erfp]);
        
    end
end

% Plot input state
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~

if Dk
    if (Image_handling_states(1,9))
        
        erfn = X(Ds+Dk,index{1})-squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{1})*Error_multiplier))';
        erfp = X(Ds+Dk,index{1})+squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{1})*Error_multiplier))';
        NMS(Ds+Dk)=state_figureInput('Neural Mass (Input)','State',fig_settings,t,[meanf*ones(1,length(t));X(Ds+Dk,index{1})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        
    end
    
    if (Image_handling_states(2,9))
        
        erfn = X(Ds+Dk,index{2})-squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{2})*Error_multiplier))';
        erfp = X(Ds+Dk,index{2})+squeeze(sqrt(Pxx(Ds+Dk,Ds+Dk,index{2})*Error_multiplier))';
        NMSZ(Ds+Dk)=state_figureInput('Neural Mass (Input)','State',fig_settings,tz,[meanf*ones(1,length(tz));X(Ds+Dk,index{2})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        
    end
    
end

% Determine whether to plot parameter estimates
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~

fig_name = {'Neural Mass Exc. Gain (G_{p})';...
    'Neural Mass Slow Inh. Gain (G_{s})';...
    'Neural Mass Fast Inh. Gain (G_{f})'};
if Dp>0
    for k = 1:Dp
        if (Image_handling_states(1,Ds+1+k))
            
            erfn = X(Ds+Dk+k,index{1})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
            erfp = X(Ds+Dk+k,index{1})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
            NMS(Ds+Dk+k)=state_figure(fig_name{k},'State',fig_settings,t,[MVI(index{3},k)';X(Ds+Dk+k,index{1})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        end
        if (Image_handling_states(2,Ds+1+k))
            erfn = X(Ds+Dk+k,index{2})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
            erfp = X(Ds+Dk+k,index{2})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
            NMSZ(Ds+Dk+k)=state_figure(fig_name{k},'State',fig_settings,tz,[MVI(index{4},k)';X(Ds+Dk+k,index{2})],{'Sim. Input','Est. Input','Std. Dev.'},[erfn;erfp]);
        end
    end
end

%%

% Plot inputs to each neural population
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_name = {'Neural Mass Membrane Potential (v_{p})';...
    'Neural Mass Membrane Potential (v_{e})';...
    'Neural Mass Membrane Potential (v_{s})';...
    'Neural Mass Membrane Potential (v_{f})'};
if any(Image_handling_inputs(1,:))
    Input(:,:,1) = [Y(index{1})'; ExpY(index{1})];
    Input(:,:,2) = [C(1)*z(index{3},1)'; C(1)*X(1,index{1})];
    Input(:,:,3) = [C(3)*z(index{3},1)'; C(3)*X(1,index{1})];
    Input(:,:,4) = [(C(5)*z(index{3},1)-C(6)*z(index{3},3))'; C(5)*X(1,index{1})-C(6)*X(3,index{1})];
end
for k = Ds/2
    if (Image_handling_inputs(1,k))
        NMSI(k)=state_figure(fig_name{k},'State',fig_settings,t,Input(:,:,k),{'Sim. Input','Est. Input'},[]);
    end
end
if any(Image_handling_inputs(2,:))
    clear Input
    Input(:,:,1) = [Y(index{2})'; ExpY(index{2})];
    Input(:,:,2) = [C(1)*z(index{4},1)'; C(1)*X(1,index{2})];
    Input(:,:,3) = [C(3)*z(index{4},1)'; C(3)*X(1,index{2})];
    Input(:,:,4) = [(C(5)*z(index{4},1)-C(6)*z(index{4},3))'; C(5)*X(1,index{2})-C(6)*X(3,index{2})];
end
for k = Ds/2
    if (Image_handling_inputs(2,k))
        NMSIZ(k)=state_figure(fig_name{k},'State',fig_settings,tz,Input(:,:,k),{'Sim. Input','Est. Input'},[]);
    end
end
clear Input

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot all aggregate membrane potentials on a single plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~

if (Image_handling_multi(1,1) ==1) % Plot all model states
    RowP = Ds/2; % Set number of rows for multi plot
    ColP =2; % Set number of columns for multi plot
    fig_name = 'Neural Mass States';
    EEG_Figure_Multi;
    namesS = {'v_{po} (mV)','v_{eo} (mV)','v_{sio} (mV)','v_{fio} (mV)' 'Z_{po} (mV)','Z_{eo} (mV)','Z_{sio} (mV)','Z_{fio} (mV)'};
    NMM(1) =figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    for k = 1:Ds
        subplot(RowP,ColP,k),plot(t,z(index{3},k),SimCol);
        hold on
        plot(t,X(k,index{1}),EstCol)
        hold on
        if (plot_uncertainty ==1)
            erfn = X(k,index{1})-squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
            erfp = X(k,index{1})+squeeze(sqrt(Pxx(k,k,index{1})*Error_multiplier))';
            plot(t,erfn,ErrCol);
            hold on
            plot(t,erfp,ErrCol);
        else
            erfn =X(k,index{1});
            erfp = erfn;
        end
        set(gca,'fontsize',tick_fontsize)
        box off
        xlabel('Time (s)','fontsize',label_fontsize)
        ylabel(namesS(k),'fontsize',label_fontsize)
        %     title('Rate of Change of Vfi','fontsize', label_fontsize)
        minc(1) = min(erfn(maxlimit)); minc(2) = min(z(index{3},k));
        maxc(1) = max(erfp(maxlimit)); maxc(2) = max(z(index{3},k));
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    m=legend('Sim.','Est.','Location',legLoc);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
end

if (Image_handling_multi(2,1) ==1) % Plot all model states
    RowP = Ds/2; % Set number of rows for multi plot
    ColP =2; % Set number of columns for multi plot
    fig_name = 'Neural Mass States Zoomed In';
    EEG_Figure_Multi;
    namesS = {'v_{po} (mV)','v_{eo} (mV)','v_{sio} (mV)','v_{fio} (mV)' 'Z_{po} (mV)','Z_{eo} (mV)','Z_{sio} (mV)','Z_{fio} (mV)'};
    NMM(1) =figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    for k = 1:Ds
        subplot(RowP,ColP,k),plot(tz,z(index{4},k),SimCol);
        hold on
        plot(tz,X(k,index{2}),EstCol)
        hold on
        if (plot_uncertainty ==1)
            erfn = X(k,index{2})-squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
            erfp = X(k,index{2})+squeeze(sqrt(Pxx(k,k,index{2})*Error_multiplier))';
            plot(tz,erfn,ErrCol);
            hold on
            plot(tz,erfp,ErrCol);
        else
            erfn =X(k,index{2});
            erfp = erfn;
        end
        set(gca,'fontsize',tick_fontsize)
        box off
        xlabel('Time (s)','fontsize',label_fontsize)
        ylabel(namesS(k),'fontsize',label_fontsize)
        %     title('Rate of Change of Vfi','fontsize', label_fontsize)
        minc(1) = min(erfn); minc(2) = min(z(index{3},k));
        maxc(1) = max(erfp); maxc(2) = max(z(index{3},k));
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    m=legend('Sim.','Est.','Location',legLoc);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
end

if (Image_handling_multi(1,2) ==1) % Plot all state inputs
    % Plot all aggregate membrane potentials on a single plot
    RowP = 2;
    ColP = 2;
    namesS = {'v_{p} (mV)','v_{e} (mV)','v_{si} (mV)','v_{fi} (mV)'};
    StateIn = [Y(index{1}) C(1)*z(index{3},1) C(3)*z(index{3},1) C(5)*z(index{3},1)-C(6)*z(index{3},3)];
    StateInE = [ExpY(index{1}); C(1)*X(1,index{1}); C(3)*X(1,index{1}); C(5)*X(1,index{1})-C(6)*X(3,index{1})];
    fig_name = 'Neural Mass State Inputs';
    EEG_Figure_Multi;
    NMM(2) =figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    legLocM = [0.35 0.92 0.3 0.1];
    for k =1:size(StateIn,2)
        subplot(RowP,ColP,k),plot(t,StateIn(:,k),SimCol);
        hold on
        plot(t,StateInE(k,:),EstCol);
        set(gca,'fontsize',tick_fontsize)
        box off
        ylabel(namesS{k},'fontsize',label_fontsize)
        %     title('Rate of Change of Vfi','fontsize', label_fontsize)
        minc(1) = min(StateInE(k,:)); minc(2) = min(StateIn(:,k));
        maxc(1) = max(StateInE(k,:)); maxc(2) = max(StateIn(:,k));
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    clear StateIn StatInE namesS
    m=legend('Sim.','Est.','Location',legLocM,'Orientation',legOri);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
    xpos = xlabel('Time (s)','fontsize',label_fontsize)
    %         set(xpos,'Position',[-4 -30 1])
end

if (Image_handling_multi(2,2) ==1) % Plot all state inputs
    % Plot all aggregate membrane potentials on a single plot
    RowP = 2;
    ColP = 2;
    namesS = {'v_{p} (mV)','v_{e} (mV)','v_{si} (mV)','v_{fi} (mV)'};
    StateIn = [Y(index{2}) C(1)*z(index{4},1) C(3)*z(index{4},1) C(5)*z(index{4},1)-C(6)*z(index{4},3)];
    StateInE = [ExpY(index{2}); C(1)*X(1,index{2}); C(3)*X(1,index{2}); C(5)*X(1,index{2})-C(6)*X(3,index{2})];
    fig_name = 'Neural Mass State Inputs Zoomed In';
    EEG_Figure_Multi;
    NMMZ(2) =figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    legLocM = [0.3 0.45 0.5 1];
    for k =1:size(StateIn,2)
        subplot(RowP,ColP,k),plot(tz,StateIn(:,k),SimCol);
        hold on
        plot(tz,StateInE(k,:),EstCol);
        set(gca,'fontsize',tick_fontsize)
        box off
        xlabel('Time (s)','fontsize',label_fontsize)
        ylabel(namesS(k),'fontsize',label_fontsize)
        %     title('Rate of Change of Vfi','fontsize', label_fontsize)
        
        minc(1) = min(StateInE(k,:)); minc(2) = min(StateIn(:,k));
        maxc(1) = max(StateInE(k,:)); maxc(2) = max(StateIn(:,k));
        axis([0 max(tz) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    clear StateIn StatInE namesS
    m=legend('Sim.','Est.','Location',legLocM,'Orientation',legOri);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
end

%%

if (Image_handling_multi(1,3) ==1) % Plot model parameters estimated
    if Dp >0
        RowP =1;
        if Dp>2
            RowP =2;
        end
        if Dp==1
            ColP =1;
        else ColP=2;
        end
        names = {'\theta_{p} (mV)','\theta_{si} (mV)','\theta_{fi} (mV)'};
        fig_name = 'Neural Mass Parameters';
        EEG_Parameter_Figure_Multi;
        NMM(3) =figure('name',fig_name,...
            'units','centimeters',...
            'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
            'papersize',[fig_width fig_height],...
            'filename',fig_dirandname,...
            'PaperPositionMode','auto');
        legLocM = [0.3 0.45 0.5 1];
        for k =1:Dp
            
            subplot(RowP,ColP,k),plot(t,MVI(index{3},k),SimCol);
            hold on
            plot(t,X(Ds+Dk+k,index{1}),EstCol)
            hold on
            if (plot_uncertainty ==1)
                erfn = X(Ds+Dk+k,index{1})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
                erfp = X(Ds+Dk+k,index{1})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{1})*Error_multiplier))';
                plot(t,erfn,ErrCol);
                hold on
                plot(t,erfp,ErrCol);
            else
                erfn =X(Ds+Dk+k,index{1});
                erfp = erfn;
            end
            set(gca,'fontsize',tick_fontsize)
            box off
            ylabel(names(k),'fontsize',label_fontsize)
            minc(1) = min(erfn(maxlimit)); minc(2) = min(MVI(index{3},k));
            maxc(1) = max(erfp(maxlimit)); maxc(2) = max(MVI(index{3},k));
            axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
        end
    end
    m=legend('Sim.','Est.','Std. Dev.','Location',legLocM,'Orientation',legOri);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
    xlabel('Time (s)','fontsize',label_fontsize);
    
end

if (Image_handling_multi(2,3) ==1) % Plot model parameters estimated
    if Dp >0
        RowP =1;
        if Dp>2
            RowP =2;
        end
        if Dp==1
            ColP =1;
        else ColP=2;
        end
        names = {'\theta_{p} (mV)','\theta_{si} (mV)','\theta_{fi} (mV)'};
        fig_name = 'Neural Mass Parameters Zoomed In';
        EEG_Figure_Multi;
        NMMZ(3) =figure('name',fig_name,...
            'units','centimeters',...
            'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
            'papersize',[fig_width fig_height],...
            'filename',fig_dirandname,...
            'PaperPositionMode','auto');
        
        legLocM = [0.3 0.45 0.5 1];
        
        for k =1:Dp
            
            subplot(RowP,ColP,k),plot(tz,MVI(index{4},k),SimCol);
            hold on
            plot(tz,X(Ds+Dk+k,index{2}),EstCol)
            hold on
            if (plot_uncertainty ==1)
                erfn = X(Ds+Dk+k,index{2})-squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
                erfp = X(Ds+Dk+k,index{2})+squeeze(sqrt(Pxx(Ds+Dk+k,Ds+Dk+k,index{2})*Error_multiplier))';
                plot(tz,erfn,ErrCol);
                hold on
                plot(tz,erfp,ErrCol);
            else
                erfn =X(Ds+Dk+k,index{2});
                erfp = erfn;
            end
            set(gca,'fontsize',tick_fontsize)
            box off
            ylabel(names(k),'fontsize',label_fontsize)
            minc(1) = min(erfn); minc(2) = min(MVI(index{4},k));
            maxc(1) = max(erfp); maxc(2) = max(MVI(index{4},k));
            axis([tstart zoomend (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
        end
    end
    m=legend('Actual','Estimated','Location',legLocM,'Orientation',legOri);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
    xlabel('Time (s)','fontsize',label_fontsize);
    
end

if ((Image_handling_multi(1,4) ==1) &&(Dp+Dk >0)) % Plot model parameters estimated
    RowP =1;
    if (Dp+Dk>2)
        RowP =2;
    end
    if (Dp+Dk==1)
        ColP =1;
    else ColP=2;
    end
    Parameters = {MVI(index{3},1) MVI(index{3},2) MVI(index{3},3) meanf*ones(length(t),1)};
    ParametersE = {X(Ds+Dk+1,index{1}) X(Ds+Dk+2,index{1}) X(Ds+Dk+3,index{1}) X(Ds+1,index{1})};
    CoV = {Pxx(Ds+Dk+1,Ds+Dk+1,index{1}) Pxx(Ds+Dk+2,Ds+Dk+2,index{1}) Pxx(Ds+Dk+3,Ds+Dk+3,index{1}) Pxx(Ds+1,Ds+1,index{1})};
    names = {'\theta_{p} (mV)','\theta_{si} (mV)','\theta_{fi} (mV)', 'Input Mean (Hz)'};
    fig_name = 'Neural Mass Parameters and Input Zoomed In';
    EEG_Parameter_Figure_Multi;
    NMM(4) =figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    
    legLocM = [0.3 0.45 0.5 1];
    
    for k =1:(Dp+Dk)
        
        subplot(RowP,ColP,k),plot(t,Parameters{k},SimCol);
        hold on
        plot(t,ParametersE{k},EstCol)
        hold on
        if (plot_uncertainty ==1)
            erfn = ParametersE{k}-squeeze(sqrt(CoV{k}*Error_multiplier))';
            erfp = ParametersE{k}+squeeze(sqrt(CoV{k}*Error_multiplier))';
            plot(t,erfn,ErrCol);
            hold on
            plot(t,erfp,ErrCol);
        else
            erfn =ParametersE{k};
            erfp = erfn;
        end
        set(gca,'fontsize',tick_fontsize)
        box off
        ylabel(names(k),'fontsize',label_fontsize)
        minc(1) = min(erfn(maxlimit)); minc(2) = min(Parameters{k});
        maxc(1) = max(erfp(maxlimit)); maxc(2) = max(Parameters{k});
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    m=legend('Sim.','Est.','Std. Dev.','Location',legLocM,'Orientation',legOri);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
    xlabel('Time (s)','fontsize',label_fontsize);
    
end

if (Image_handling_multi(2,4) ==1) % Plot model parameters estimated
    if (Dp+Dk >0)
        RowP =1;
        if (Dp+Dk>2)
            RowP =2;
        end
        if (Dp+Dk==1)
            ColP =1;
        else ColP=2;
        end
        Parameters = {MVI(index{4},1) MVI(index{4},2) MVI(index{4},3) meanf*ones(length(t),1)};
        ParametersE = {X(Ds+Dk+1,index{2}) X(Ds+Dk+2,index{2}) X(Ds+Dk+3,index{2}) X(Ds+1,index{2})};
        CoV = {Pxx(Ds+Dk+1,Ds+Dk+1,index{2}) Pxx(Ds+Dk+2,Ds+Dk+2,index{2}) Pxx(Ds+Dk+3,Ds+Dk+3,index{2}) Pxx(Ds+1,Ds+1,index{2})};
        names = {'A (mV)','B (mV)','G (mV)' 'Input Mean (Hz)'};
        fig_name = 'Neural Mass Parameters and Input Zoomed In';
        EEG_Figure_Multi;
        NMMZ(4) =figure('name',fig_name,...
            'units','centimeters',...
            'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
            'papersize',[fig_width fig_height],...
            'filename',fig_dirandname,...
            'PaperPositionMode','auto');
        
        
        for k =1:(Dp+Dk)
            
            subplot(RowP,ColP,k),plot(tz,Parameters{k},SimCol);
            hold on
            plot(tz,ParametersE{k},EstCol)
            hold on
            if (plot_uncertainty ==1)
                erfn = ParametersE{k}-squeeze(sqrt(CoV{k}*Error_multiplier))';
                erfp = ParametersE{k}+squeeze(sqrt(CoV{k}*Error_multiplier))';
                plot(tz,erfn,ErrCol);
                hold on
                plot(tz,erfp,ErrCol);
            else
                erfn =ParametersE{k};
                erfp = erfn;
            end
            set(gca,'fontsize',tick_fontsize)
            box off
            ylabel(names(k),'fontsize',label_fontsize)
            minc(1) = min(erfn); minc(2) = min(Parameters{k});
            maxc(1) = max(erfp); maxc(2) = max(Parameters{k});
            axis([tstart zoomend (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
        end
    end
    m=legend('Actual','Estimated','Location',legLocM,'Orientation',legOri);
    legend(m,'boxoff');
    set(m,'fontsize',legend_fontsize);
    xlabel('Time (s)','fontsize',label_fontsize);
    
end


%%

if (Image_handling_firing_rates(1) ==1)
    % Plot of state 1
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RowP =3;
    fig_name = 'Neural Mass Pyramidal Neuron Input, Output and Firing Rate';
    EEG_Figure_Multi;
    MassP = {Y(index{1}); Sigmoid(Y(index{1})); z(index{3},1)};
    MassPE = {ExpY(index{1}); Sigmoid(ExpY(index{1})); X(1,index{1})};
    ylab = {'Pyramidal Input (mV)' 'Pyramidal Firing Rate (Hz)' 'Pyramidal Output (mV)'};
    NMFR(1) = figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    for k =1:3
        subplot(RowP,1,k),plot(t,MassP{k},SimCol);
        hold on
        plot(t,MassPE{k},EstCol);
        set(gca,'fontsize',tick_fontsize)
        box off
        ylabel(ylab(k),'fontsize',label_fontsize)
        minc(1) = min(MassPE{k}); minc(2) = min(MassP{k});
        maxc(1) = max(MassPE{k}); maxc(2) = max(MassP{k});
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    k=legend('Sim.','Est.','Location',legLoc);
    legend(k,'boxoff');
    set(k,'fontsize',legend_fontsize);
    % Print image with resolution 600 to pdf test.
end

if (Image_handling_firing_rates(2) ==1)
    
    % Plot of state 2
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RowP =3;
    fig_name = 'Neural Mass Excitatory Neuron Input, Output and Firing Rate';
    EEG_Figure_Multi;
    MassP = {C(1)*z(index{3},1); Sigmoid(C(1)*z(index{3},1)); z(index{3},2)};
    MassPE = {C(1)*X(1,index{3}); Sigmoid(C(1)*X(1,index{3})); X(2,index{3})};
    ylab = {'Excitatory Input (mV)' 'Excitatory Firing Rate (Hz)' 'Excitatory Output (mV)'};
    NMFR(2) = figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    for k =1:3
        subplot(RowP,1,k),plot(t,MassP{k},SimCol);
        hold on
        plot(t,MassPE{k},EstCol);
        set(gca,'fontsize',tick_fontsize)
        box off
        ylabel(ylab(k),'fontsize',label_fontsize)
        minc(1) = min(MassPE{k}); minc(2) = min(MassP{k});
        maxc(1) = max(MassPE{k}); maxc(2) = max(MassP{k});
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    k=legend('Sim.','Est.','Location',legLoc);
    legend(k,'boxoff');
    set(k,'fontsize',legend_fontsize);
    % Print image with resolution 600 to pdf test.
end

if (Image_handling_firing_rates(3) ==1)
    
    % Plot of state 3
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RowP =3;
    fig_name = 'Neural Mass Excitatory Neuron Input, Output and Firing Rate';
    EEG_Figure_Multi;
    MassP = {C(3)*z(index{3},1); Sigmoid(C(3)*z(index{3},1)); z(index{3},3)};
    MassPE = {C(3)*X(1,index{3}); Sigmoid(C(3)*X(1,index{3})); X(3,index{3})};
    ylab = {'Slow Inh. Input (mV)' 'Slow Inh. Firing Rate (Hz)' 'Slow Inh. Output (mV)'};
    NMFR(3) = figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    for k =1:3
        subplot(RowP,1,k),plot(t,MassP{k},SimCol);
        hold on
        plot(t,MassPE{k},EstCol);
        set(gca,'fontsize',tick_fontsize)
        box off
        ylabel(ylab(k),'fontsize',label_fontsize)
        minc(1) = min(MassPE{k}); minc(2) = min(MassP{k});
        maxc(1) = max(MassPE{k}); maxc(2) = max(MassP{k});
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    k=legend('Sim.','Est.','Location',legLoc);
    legend(k,'boxoff');
    set(k,'fontsize',legend_fontsize);
    % Print image with resolution 600 to pdf test.
end


if (Image_handling_firing_rates(4) ==1)
    
    % Plot of state 4
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RowP =3;
    fig_name = 'Neural Mass Fast Inh. Neuron Input, Output and Firing Rate';
    EEG_Figure_Multi;
    MassP = {C(5)*z(index{3},1)-C(6)*z(index{3},3); Sigmoid(C(5)*z(index{3},1)-C(6)*z(index{3},3)); z(index{3},4)};
    MassPE = {C(5)*X(1,index{3})-C(6)*X(3,index{3}); Sigmoid(C(5)*X(1,index{3})-C(6)*X(3,index{3})); X(4,index{3})};
    ylab = {'Fast Inh. Input (mV)' 'Fast Inh. Firing Rate (Hz)' 'Fast Inh. Output (mV)'};
    NMFR(4) = figure('name',fig_name,...
        'units','centimeters',...
        'position',[fig_left_pos fig_bottom_pos fig_width fig_height],...
        'papersize',[fig_width fig_height],...
        'filename',fig_dirandname,...
        'PaperPositionMode','auto');
    for k =1:3
        subplot(RowP,1,k),plot(t,MassP{k},SimCol);
        hold on
        plot(t,MassPE{k},EstCol);
        set(gca,'fontsize',tick_fontsize)
        box off
        ylabel(ylab(k),'fontsize',label_fontsize)
        minc(1) = min(MassPE{k}); minc(2) = min(MassPE{k});
        maxc(1) = max(MassPE{k}); maxc(2) = max(MassP{k});
        axis([0 max(t) (min(minc)-abs(min(minc))*scale) (max(maxc)+abs(max(maxc))*scale)]);
    end
    k=legend('Sim.','Est.','Location',legLoc);
    legend(k,'boxoff');
    set(k,'fontsize',legend_fontsize);
    % Print image with resolution 600 to pdf test.
end

