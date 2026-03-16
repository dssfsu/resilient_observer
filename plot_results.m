%% Creating visualization and plotting results for the simulation

% Olugbenga Moses Anubi, Jan 2020


%% Extracting values
time_vec  = out.logsout.getElement('x').Values.Time;

% State vectors
x         = out.logsout.getElement('x').Values.Data; % true states
x_hat_LO  = out.logsout.getElement('x_hat_LO').Values.Data; % Observed states (Luenberger)
x_hat_L1O = out.logsout.getElement('x_hat_L1O').Values.Data; % Observed states (Unconstrained l1)
x_hat_MMO = out.logsout.getElement('x_hat_MMO').Values.Data; % Observed states (Constrained MMO)

% Auxiliary model and FDIA performance
% y_obsv     = out.logsout.getElement('y_obsv').Values.Data; % measured y for observer (y_obsv = y-D*u)
% mu_k       = out.logsout.getElement('mu_k').Values.Data; % mean value for auxiliary variables
aux_model_res  = out.logsout.getElement('y_obsv').Values.Data; % measured y for observer (y_obsv = y-D*u)
BDD_res        = out.logsout.getElement('BDD_res').Values.Data; % Bad data residue


%% Ploting/Visualization/Metric tables

%% 1. delta estimates
n_delta = size(x,2)/2;  % number of generator angles

figure(1) % Luenberger Observer
figure(2) % Unconstrained l_1 based minimization
figure(3) % Multimodel

LW = 1.3;  % linewidth
FS = 15;   % font size

if n_delta > 5
    % Plot all angles together

    time_select_mask = time_vec >= T_start_opt; % ignoring the optimizer idle time

    % Luenberger
    figure(1)
    subplot(3,1,1)
    plot(time_vec,x,'LineWidth',LW)
    ylabel('\delta','FontWeight','bold')
    title('Actual Angles')

    subplot(3,1,2)
    plot(time_vec,x_hat_LO,'LineWidth',LW)
    ylabel('\delta_{est}','FontWeight','bold')
    title('Estimated Angles')

    subplot(3,1,3)
    plot(time_vec,(x - x_hat_LO).*time_select_mask,'LineWidth',LW)
    ylabel('\delta - \delta_{est}','FontWeight','bold')
    title('Estimation Error') 
    sgtitle('Luenberger Observer')

    % L1 unconstrained
    figure(2)
    subplot(3,1,1)
    plot(time_vec,x,'LineWidth',LW)
    ylabel('\delta','FontWeight','bold')
    title('Actual Angles')

    subplot(3,1,2)
    plot(time_vec,x_hat_L1O,'LineWidth',LW)
    ylabel('\delta_{est}','FontWeight','bold')
    title('Estimated Angles')

    subplot(3,1,3)
    plot(time_vec,(x - x_hat_L1O).*time_select_mask,'LineWidth',LW)
    ylabel('\delta - \delta_{est}','FontWeight','bold')
    title('Estimation Error') 
    sgtitle('Unconstrained L1 Observer')

    % MMO
    figure(3)
    subplot(3,1,1)
    plot(time_vec,x,'LineWidth',LW)
    ylabel('\delta','FontWeight','bold')
    title('Actual Angles')

    subplot(3,1,2)
    plot(time_vec,x_hat_MMO,'LineWidth',LW)
    ylabel('\delta_{est}','FontWeight','bold')
    title('Estimated Angles')

    subplot(3,1,3)
    plot(time_vec,(x - x_hat_MMO).*time_select_mask,'LineWidth',LW)
    ylabel('\delta - \delta_{est}','FontWeight','bold')
    title('Estimation Error') 
    sgtitle('Multi-Model Observer')



else
    % Plot all angles separately
    for iter = 1:n_delta

        % Luenberger
        figure(1)
        subplot(n_delta,2,2*(iter-1)+1)
        plot(time_vec,x(:,iter),'k','LineWidth',LW)
        ylabel(['\delta_{' num2str(iter) '}'],'FontWeight','bold')
        if(iter == 1)
            title('Actual')
        end

        subplot(n_delta,2,2*iter)
        plot(time_vec,x_hat_LO(:,iter),'k','LineWidth',LW)
        ylabel(['\delta_{' num2str(iter) '}'],'FontWeight','bold')

        if(iter == 1)
            title('Luenberger Observer')
        end

        % L1 unconstrained
        figure(2)
        subplot(n_delta,2,2*(iter-1)+1)
        plot(time_vec,x(:,iter),'k','LineWidth',LW)
        ylabel(['\delta_{' num2str(iter) '}'],'FontWeight','bold')
        if(iter == 1)
            title('Actual')
        end

        subplot(n_delta,2,2*iter)
        plot(time_vec,x_hat_L1O(:,iter),'k','LineWidth',LW)
        ylabel(['\delta_{' num2str(iter) '}'],'FontWeight','bold')
        if(iter == 1)
            title('Unconstrained l_1-based Observer')
        end


        % MMO
        figure(3)
        subplot(n_delta,2,2*(iter-1)+1)
        plot(time_vec,x(:,iter),'k','LineWidth',LW)
        ylabel(['\delta_{' num2str(iter) '}'],'FontWeight','bold')
        if(iter == 1)
            title('Actual')
        end

        subplot(n_delta,2,2*iter)
        plot(time_vec,x_hat_MMO(:,iter),'k','LineWidth',LW)
        ylabel(['\delta_{' num2str(iter) '}'],'FontWeight','bold')
        if(iter == 1)
            title('Multi-Model Observer')
        end

    end
end

%% BDD_res
figure,
plot(time_vec,BDD_res,'k','LineWidth',LW), hold on
plot(time_vec,BDD_thresh*ones(length(time_vec),1),'r--','LineWidth',2*LW)
ylabel('Bad Data Detection Residual','FontWeight','bold')



%% Error table
error_LO = x - x_hat_LO;
error_LO_delta = error_LO(:,1:n_delta).'; % NOTE THE TRANSPOSE

error_L1O = x - x_hat_L1O;
error_L1O_delta = error_L1O(:,1:n_delta).';% NOTE THE TRANSPOSE

error_MMO = x - x_hat_MMO;
error_MMO_delta = error_MMO(:,1:n_delta).';% NOTE THE TRANSPOSE

metric_rms = @(x,dim) sqrt(sum(x.^2,dim)/size(x,dim));
metric_max = @(x,dim) max(abs(x),[],dim);
error_table = [metric_rms(error_LO_delta,2) metric_rms(error_L1O_delta,2) metric_rms(error_MMO_delta,2) ...
               metric_max(error_LO_delta,2) metric_max(error_L1O_delta,2) metric_max(error_MMO_delta,2)];


