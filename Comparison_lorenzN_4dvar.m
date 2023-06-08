%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program to perform 4D-Var data assimilation on the
% Lorenz N problem - error in parameters as model error
% as in 2008 by M. A. Freitag
% 
% Comparison with Time-limited balanced truncation (TLBT) model order
% reduction -- 2023 J. Koenig
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables:
%
%
%    force, N:   parameters
%
%    time_assi:  Length of assimilation window
%    time_fore:  Length of forecast window
%    h:          Time step for numerical scheme
%    freq:       Frequency of observations in time steps
%    tstep_assi: Number of time steps for assimilation
%    tstep_fore: Total number of time steps for assimilation + forecast
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------------------------------------------
%   1. Set up the problem
%---------------------------------------------------------------------
clear all
close all
format long

% set random number generator
rng('default')
rng(1)

% F and system size
force   = 8;
N       = 500;

% initial conditions; collect all the intial conditions in a vector
guessXtruth = randn(N,1);

% choice of observed points
choice = menu('Choose observations','In all N variables','every s variables');
if choice == 1
    % observation matrix is constant and identity
    H           = speye(N);
    s           = 1;
    obs_vars    = N;
    last_ind    = N;
else
    s   = input('Take observations every ... variables. ');
    obs_vars    = ceil(N/s);
    H           = zeros(obs_vars,N);
    for i = 1:obs_vars
        H(i,(i-1)*s+1)  = 1;
    end
    last_ind    = (obs_vars-1)*s+1;
    H           = sparse(H);
end

% Choice of prior covariance matrix
choiceB = menu('Covariance matrix','Scaled Identity','Gaussian exponential');
sb = 0.1;
if choiceB == 1
    % Identity
    %sb = input('Multiple of Identity: ');
    B = sb*speye(N);
else
    %sb = input('Multiple of Gaussian Exponential: ');
    Ls = 1;
    for i = 1:N
        for j = 1:N
        B(i,j) = sb*exp((-abs(i-j))/(2*Ls^2));
        end
    end
end

%---------------------------------------------------------------------
%  2. Input true parameters and calculate truth trajectory
%---------------------------------------------------------------------

% Inputs 
tstep_assi  = 50;   % assimilation window (data collection)
tstep_fore  = 100;  % forecast window (no observable data)
h           = 0.01; % time step between observation times
freq        = 1;    % frequency of observations

tstep_truth = tstep_assi +  tstep_fore; % total time window
endtime     = tstep_assi * h;           % end time of observations

% Find truth trajectory for the whole length of time
[Xtruth]    = rk4(tstep_truth,h,guessXtruth,force);
 
% plot true solution
figure(1)

xvals = 0:tstep_truth;
clf;
subplot(2,1,1)
plot(xvals,Xtruth(1,:),'g:','LineWidth',2)
hold on
xlabel('time step')
ylabel('x_1')
title('x_1')
legend('Truth x_1')
%
subplot(2,1,2)
plot(xvals,Xtruth(last_ind,:),'g:','LineWidth',2)
hold on
xlabel('time step')
ylabel('x_{last}')
title('x_{last}')
legend('Truth x_{last}')

% initial conditions error
guessX = guessXtruth + sqrt(B)*randn(N,1);

% calculate trajectory for imperfect model
[Xtruthimperfect] = rk4(tstep_truth,h,guessX,force);

%---------------------------------------------------------------------
% 4D-Var
%--------------------------------------------------------------------- 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for simplicity we use observation frequency 1;

% number of observations; observations only in data assimilations window 
n_obs       = ceil(tstep_assi/freq); 

% covariance of observational noise
R           = 0.01*diag(ones(obs_vars,1));  
stand_dev   = sqrt(R);

RX  = zeros(obs_vars,n_obs); % vector of observational noise     
for n = 1:obs_vars
    noise   = randn(1,n_obs);
    RX(n,:) = stand_dev(n,n)*noise;
end
% generate observations
obX = H * Xtruth(:,2:freq:tstep_assi+1) + RX; % vector of observations  

% Plot truth and observations
vec = 1:freq:tstep_assi;

% temporary arrays for plotting qx1_plot, qx2_plot
% which coordinates are really observed??? fix this all throughout!
qx1_plot=obX(1,:);
qxlast_plot=obX(end,:);

figure(1)
subplot(2,1,1)
hold on
plot(vec,qx1_plot,'om')
hold on
xlabel('time step')
ylabel('x_1')
title('x_1')
legend('Truth x_1','Observations for 4DVar')

%
subplot(2,1,2)
hold on
plot(vec,qxlast_plot,'om')
hold on
xlabel('time step')
ylabel('x_{last}')
title('x_{last}')
legend('Truth x_{last}','Observations for 4DVar')

disp('* ... observations for 4DVar done')

% finally start 4D-Var (at step i = 1)
% initial guess for 4D-Var

currentX = Xtruthimperfect;
backgroundX = currentX;

% Plot first guess
figure(1)
xvals = 0:tstep_truth;
subplot(2,1,1)
hold on
plot(xvals,backgroundX(1,:),'Color','#7E2F8E','LineStyle','--')
xlabel('time step')
ylabel('x_1')
title('x_1')
legend('Truth x_1','Observations for 4DVar','Initial guess 4DVar','Location','best')
subplot(2,1,2)
hold on
plot(xvals,backgroundX(last_ind,:),'Color','#7E2F8E','LineStyle','--')
xlabel('time step')
ylabel('x_{last}')
title('x_{last}')
legend('Truth x_{last}','Observations for 4DVar','Initial guess 4DVar','Location','best')

disp('* ... initial guess for 4DVar done')

% plot rms error in background over time 
figure(2)
Z = backgroundX(:,:)-Xtruth(:,:);
semilogy(xvals,sqrt(sum(Z.*Z,1)),'Color','#7E2F8E','LineStyle','--','LineWidth',2);
hold on
legend('before assimilation','Location','best')
xlabel('time step')
ylabel('RMS error to truth')

currentX0 = currentX(:,1);
firstX_MOR = currentX;

% set parameters for Gauss-Newton method
maxit       = N;
its         = 0;
tolerance   = 10^-8;

%% 4D-Var cost function minimisation for the full model
tic
% calculate cost function, Jacobian and Hessian for initial state
[J,Jgrad,JHess] = costfunction_fixed_H(tstep_assi,h, ...
    backgroundX(:,1),currentX(:,1:tstep_assi+1),obX,B,R,H,freq,force);

% collect cost function and gradient values
allJ        = J;
allJgrad    = norm(Jgrad);

% perform minimisation using Gauss-Newton method

while norm(Jgrad) >= tolerance && its < maxit
    
    % solve system
    increment   = - JHess\Jgrad;
       
    % update initial guess X_0
    currentX0   = currentX0 + increment;

    clear currentX
    
    % perform forecast on the new initial guess to obtain a current state
    [currentX]  = rk4(tstep_assi,h,currentX0,force);
    
    % calculate cost function, Jacobian and Hessian for iterate
    [J,Jgrad,JHess] = costfunction_fixed_H(tstep_assi,h,backgroundX(:,1), ...
        currentX(:,1:tstep_assi+1),obX,B,R,H,freq,force);

    % just for monitoring purposes of cost function J and its gradient
    its = its+1
    J
    norm(Jgrad)
    
    % collect cost function and gradient values
    allJ        = [allJ;  J];
    allJgrad    = [allJgrad; norm(Jgrad)];
    
end
runtime_total_full = toc;
% end minimisation

%---------------------------------------------------------------------
% 6. Finally, run forecast from final analysis and plot it
%---------------------------------------------------------------------

[Xfinal] = rk4(tstep_assi+tstep_fore,h,currentX0,force);

% calculate final x_1 error
disp('final x_1 error at end of assimilation window')
disp('full model')
finalx1 = norm(Xfinal(1,tstep_assi+1)-Xtruth(1,tstep_assi+1))

% calculate final x_{last} error
disp('final x_{last} error at end of assimilation window')
disp('full model')
finalx2 = norm(Xfinal(last_ind,tstep_assi+1)-Xtruth(last_ind,tstep_assi+1))

disp('final error in norm at end of assimilation window')
disp('full model')
finalerror = norm(Xfinal(:,tstep_assi+1)-Xtruth(:,tstep_assi+1))

disp('RMS error over the whole assimilation window')
disp('full model')
assi_error = norm(Xfinal(:,1:tstep_assi+1)-Xtruth(:,1:tstep_assi+1),'fro')

%% 4D-Var cost function minimisation for the TLBT-reduced model
% rank of reduced model
r = 3*ceil(N/4);

allJ_MOR            = zeros(length(r),maxit);
allJgrad_MOR        = zeros(length(r),maxit);
allJ_MOR_a          = zeros(length(r),maxit);
allJgrad_MOR_a      = zeros(length(r),maxit);

for k = 1:length(r)

    % set staring values
    currentX0_MOR   = firstX_MOR(:,1);
    currentX_MOR    = firstX_MOR;
    its             = 0;

    tic
    % compute constant system matrix (middle of observation window)
    M   = JacobianRK4(h*freq,currentX_MOR(:,(ceil(n_obs/2)-1)*freq+1),force);

    % compute the reduced order model via TLBT
    tic
    [JHess,B_r,H_r,XB_r,Sr_TL,G_BTQ_TL] = Simple_TLBT_Quantities_lorenz(r(k), ...
        N-1, n_obs,endtime,B,R,M,H,backgroundX(:,1));
    time_build_TLBT = toc;

    tic
    % calculate cost function, Jacobian and Hessian for initial state
    [J,Jgradient_temp] = costfunction_Gradient_reduced(n_obs,freq, ...
        currentX_MOR(:,1:tstep_assi+1),obX,B_r,XB_r,Sr_TL,R,H_r,G_BTQ_TL);
    % add initial state correction to the gradient
    Jgrad = (B\(currentX_MOR(:,1)-backgroundX(:,1))) - Jgradient_temp;

    % collect cost function and gradient values
    allJ_MOR(k,its+1)        = J;
    allJgrad_MOR(k,its+1)    = norm(Jgrad);

    % perform minimisation using Gauss-Newton method (constant Hessian)

    while norm(Jgrad) >= tolerance && its < maxit
    
        % solve system
        increment       = - JHess\Jgrad;
       
        % update initial guess X_0
        currentX0_MOR   = currentX0_MOR + increment;

        clear currentX_MOR
    
        % perform forecast on the new initial guess to obtain a current state
        [currentX_MOR]  = rk4(tstep_assi,h,currentX0_MOR,force);
    
        % calculate cost function, Jacobian and Hessian for iterate
        [J,Jgradient_temp]  = costfunction_Gradient_reduced(n_obs, ...
        freq,currentX_MOR(:,1:tstep_assi+1),obX,B_r,XB_r,Sr_TL,R,H_r,G_BTQ_TL);
        % add initial state correction to the gradient
        Jgrad = (B\(currentX_MOR(:,1)-backgroundX(:,1))) - Jgradient_temp;

        % just for monitoring purposes of cost function J and its gradient
        its = its+1
        J
        norm(Jgrad)
    
        % collect cost function and gradient values
        allJ_MOR(k,its+1)        = J;
        allJgrad_MOR(k,its+1)    = norm(Jgrad);
    end
    time_reduction_TLBT     = toc;
    % end minimisation
    
    %---------------------------------------------------------------------
    % 6. Finally, run forecast from final analysis and plot it
    %---------------------------------------------------------------------
    [Xfinal_MOR] = rk4(tstep_assi+tstep_fore,h,currentX0_MOR,force);

    %% alpha-bounded balancing
    % set staring values
    currentX0_MOR   = firstX_MOR(:,1);
    currentX_MOR    = firstX_MOR;
    its             = 0;

    tic
    % compute the reduced order model via alpha-BT
    [JHess,B_r,H_r,XB_r,Sr_TL,G_BTQ_TL] = alpha_bounded_BT_Quantities_lorenz(r(k), ...
        N-1, n_obs,endtime,B,R,M,H,backgroundX(:,1));
    time_build_BT_a = toc;

    tic
    % calculate cost function, Jacobian and Hessian for initial state
    [J,Jgradient_temp] = costfunction_Gradient_reduced(n_obs,freq, ...
        currentX_MOR(:,1:tstep_assi+1),obX,B_r,XB_r,Sr_TL,R,H_r,G_BTQ_TL);
    % add initial state correction to the gradient
    Jgrad = (B\(currentX_MOR(:,1)-backgroundX(:,1))) - Jgradient_temp;

    % collect cost function and gradient values
    allJ_MOR_a(k,its+1)        = J;
    allJgrad_MOR_a(k,its+1)    = norm(Jgrad);

    % perform minimisation using Gauss-Newton method (constant Hessian)

    while norm(Jgrad) >= tolerance && its < maxit
    
        % solve system
        increment       = - JHess\Jgrad;
       
        % update initial guess X_0
        currentX0_MOR   = currentX0_MOR + increment;

        clear currentX_MOR
    
        % perform forecast on the new initial guess to obtain a current state
        [currentX_MOR]  = rk4(tstep_assi,h,currentX0_MOR,force);
    
        % calculate cost function, Jacobian and Hessian for iterate
        [J,Jgradient_temp]  = costfunction_Gradient_reduced(n_obs, ...
        freq,currentX_MOR(:,1:tstep_assi+1),obX,B_r,XB_r,Sr_TL,R,H_r,G_BTQ_TL);
        % add initial state correction to the gradient
        Jgrad = (B\(currentX_MOR(:,1)-backgroundX(:,1))) - Jgradient_temp;

        % just for monitoring purposes of cost function J and its gradient
        its = its+1
        J
        norm(Jgrad)
    
        % collect cost function and gradient values
        allJ_MOR_a(k,its+1)        = J;
        allJgrad_MOR_a(k,its+1)    = norm(Jgrad);
    end
    time_reduction_BT_a     = toc;
    % end minimisation
    
    %---------------------------------------------------------------------
    % 6. Finally, run forecast from final analysis and plot it
    %---------------------------------------------------------------------

    [Xfinal_MOR_a] = rk4(tstep_assi+tstep_fore,h,currentX0_MOR,force);

    % Plot
    figure(1)
    xvals = 0:tstep_truth;
    subplot(2,1,1)
    hold on
    plot(xvals,Xfinal_MOR(1,:),'color',[0,0.6*(k-1)/length(r), 0.6/k],'DisplayName',['TLBT-4DVar analysis with r = ', num2str(r(k))]);
    plot(xvals,Xfinal_MOR_a(1,:),'color',[0.6*(k-1)/length(r), 0.6/k,0],'DisplayName',['alpha-BT-4DVar analysis with r = ', num2str(r(k))]);
    xlabel('time step')
    ylabel('x_1')
    title('x_1')
    % plot vertical line at beginning of forecast
    if k == length(r)
        plot(xvals,Xfinal(1,:),'color',[0.7,0,0],'LineStyle','--','DisplayName','4DVar analysis')
        yL = get(gca,'YLim');
        line([tstep_assi tstep_assi],yL,'Color','k','DisplayName','end of assimilation window');
    end

    subplot(2,1,2)
    hold on
    plot(xvals,Xfinal_MOR(last_ind,:),'color',[0,0.6*(k-1)/length(r), 0.6/k],'DisplayName',['TLBT-4DVar analysis with r = ', num2str(r(k))]);
    plot(xvals,Xfinal_MOR_a(last_ind,:),'color',[0.6*(k-1)/length(r), 0.6/k,0],'DisplayName',['alpha-BT-4DVar analysis with r = ', num2str(r(k))]);
    xlabel('time step')
    ylabel('x_{last}')
    title('x_{last}')
    % plot vertical line at beginning of forecast
    if k == length(r)
        plot(xvals,Xfinal(last_ind,:),'color',[0.7,0,0],'LineStyle','--','DisplayName','4DVar analysis');
        yL = get(gca,'YLim');
        line([tstep_assi tstep_assi],yL,'Color','k','DisplayName','end of assimilation window');
    end

    % plot rms error in analysis over time
    figure(2)
    clear Z Z_MOR Z1 Z2 Z1_MOR Z2_MOR
    Z = Xfinal(:,:)-Xtruth(:,:);
    Z_MOR = Xfinal_MOR(:,:)-Xtruth(:,:);
    semilogy(xvals,sqrt(sum(Z_MOR.*Z_MOR,1)),'color',[0,0.6*(k-1)/length(r), 0.6/k],'LineWidth',2,'DisplayName',['after assim., TLBT, r = ', num2str(r(k))]);
    Z_MOR_a = Xfinal_MOR_a(:,:)-Xtruth(:,:);
    semilogy(xvals,sqrt(sum(Z_MOR_a.*Z_MOR_a,1)),'color',[0.6*(k-1)/length(r), 0.6/k,0],'Marker','x','DisplayName',['after assim., alpha-BT, r = ', num2str(r(k))]);
    if k == length(r)
        semilogy(xvals,sqrt(sum(Z.*Z,1)),'color',[0.7,0,0],'LineStyle','-.','LineWidth',2,'DisplayName','after assim.');
        % plot vertical line at beginning of forecast
        yL = get(gca,'YLim');
        line([tstep_assi tstep_assi],yL,'Color','k','DisplayName','end of assimilation window');
        xlabel('time step')
        ylabel('RMS error to truth')
        title('RMS error over all variables over time')
    end

    
    figure(3)
    ZS = [Z(1:s:end,:)]; % observed variables
    ZS_MOR = [Z_MOR(1:s:end,:)];
    ZS_MOR_a = [Z_MOR_a(1:s:end,:)];
    vecdiff = setdiff(1:N,1:s:N);
    ZP = [Z(vecdiff,:)]; % unobserved variable
    ZP_MOR = [Z_MOR(vecdiff,:)];
    ZP_MOR_a = [Z_MOR_a(vecdiff,:)];
    semilogy(xvals,sqrt(sum(ZS.*ZS,1)),'color',[0.7,0,0],'LineStyle','-.','LineWidth',2);
    hold on
    legend('observed, full',Location='best')
    semilogy(xvals,sqrt(sum(ZS_MOR_a.*ZS_MOR_a,1)),'color',[0.6*(k-1)/length(r), 0.6/k,0],'Marker','x','DisplayName',['observed, alpha-BT, r = ', num2str(r(k))]);
    semilogy(xvals,sqrt(sum(ZS_MOR.*ZS_MOR,1)),'color',[0,0.6*(k-1)/length(r), 0.6/k],'LineWidth',2,'DisplayName',['observed, TLBT, r = ', num2str(r(k))]);
    semilogy(xvals,sqrt(sum(ZP.*ZP,1)),'color',[0.7,0,0],'LineStyle',':' ,'LineWidth',2,'DisplayName','unobserved, full');
    semilogy(xvals,sqrt(sum(ZP_MOR_a.*ZP_MOR_a,1)),'color',[0.6*(k-1)/length(r), 0.6/k,0],'Marker','o','LineStyle','--','DisplayName',['unobserved, alpha-BT, r = ', num2str(r(k))]);
    semilogy(xvals,sqrt(sum(ZP_MOR.*ZP_MOR,1)),'color',[0,0.6*(k-1)/length(r), 0.6/k],'LineStyle','--','LineWidth',2,'DisplayName',['unobserved, TLBT, r = ', num2str(r(k))]);
    title('RMS errors in observed and unobserved variables')

    % calculate final x_1 error
    disp('final x_1 error at end of assimilation window')
    disp('with TLBT')
    finalx1_MOR(k) = norm(Xfinal_MOR(1,tstep_assi+1)-Xtruth(1,tstep_assi+1))
    disp('with alpha-BT')
    finalx1_MOR_a(k) = norm(Xfinal_MOR_a(1,tstep_assi+1)-Xtruth(1,tstep_assi+1))

    % calculate final x_last_ind error
    disp('final x_{last} error at end of assimilation window')
    disp('with TLBT')
    finalx2_MOR(k) = norm(Xfinal_MOR(last_ind,tstep_assi+1)-Xtruth(last_ind,tstep_assi+1))
    disp('with alpha-BT')
    finalx2_MOR_a(k) = norm(Xfinal_MOR_a(last_ind,tstep_assi+1)-Xtruth(last_ind,tstep_assi+1))

    disp('final error in norm at end of assimilation window')
    disp('with TLBT')
    finalerror_MOR(k) = norm(Xfinal_MOR(:,tstep_assi+1)-Xtruth(:,tstep_assi+1))
    disp('with alpha-BT')
    finalerror_MOR_a(k) = norm(Xfinal_MOR_a(:,tstep_assi+1)-Xtruth(:,tstep_assi+1))

    disp('RMS error over the whole assimilation window')
    disp('with TLBT')
    assi_error_MOR(k) = norm(Xfinal_MOR(:,1:tstep_assi+1)-Xtruth(:,1:tstep_assi+1),'fro')
    disp('with alpha-BT')
    assi_error_MOR_a(k) = norm(Xfinal_MOR_a(:,1:tstep_assi+1)-Xtruth(:,1:tstep_assi+1),'fro')
    
end    

% plot cost function and its gradient
figure(4)
semilogy(allJ,'-r','LineWidth',2);
legend('full model 4DVar')
hold all
for k = 1:length(r)
    semilogy(allJ_MOR(k,:),'color',[0,0.6*(k-1)/length(r), 0.6/k],'LineWidth',2,'DisplayName',['TLBT-4DVar with r = ', num2str(r(k))]);
    semilogy(allJ_MOR_a(k,:),'color',[0.6*(k-1)/length(r), 0.6/k,0],'LineWidth',2,'DisplayName',['alpha-BT-4DVar with r = ', num2str(r(k))]);
end
title('cost function')
xlabel('iterations')
ylabel('cost function')

figure(5)
semilogy(allJgrad,'-r','LineWidth',2);
hold all
legend('full model 4DVar')
for k = 1:length(r)
   semilogy(allJgrad_MOR(k,:),'color',[0,0.6*(k-1)/length(r), 0.6/k],'LineWidth',2,'DisplayName',['TLBT-4DVar with r = ', num2str(r(k))]);
   semilogy(allJgrad_MOR_a(k,:),'color',[0.6*(k-1)/length(r), 0.6/k,0],'LineWidth',2,'DisplayName',['alpha-BT-4DVar with r = ', num2str(r(k))]);
end
title('cost function gradient')
xlabel('iterations')
ylabel('cost function gradient')
