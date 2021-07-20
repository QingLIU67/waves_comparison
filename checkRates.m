function checkRates(time,Q,R,D,kappaFun,lambdaFun,kappa,lambda)
% function checkRates(time,Q,D,R,kappaFun,lambdaFun) compares the fitted
% and calcualted death and recovered ratios. The idea is to check whether
% the approximation of these ratios is appropriate
%
% Inputs
%   time: datetime: [1xN]: time array
%   Q: double [1xN]: Time histories of the quarantined/active cases :
%   D: double [1xN]: Time histories of the deceased cases :
%   R: double [1xN]: Time histories of the recovered cases :
%   kappaFun: anonymous function approximating the death rate
%   lambdaFun: anonymous function approximating the recovery rate
%
% Outputs:
% None
%
% Author: E. Cheynet - UiB - last modified 07-05-2020
%
% see also SEIQRDP.m fit_SEIQRDP.m

%% Compute the rate of deceased and recovered cases

Q = Q(:);
R = R(:);
D = D(:);
time = time(:);

rateD = (diff(D)./diff(datenum(time-time(1))))./Q(2:end);
rateD(abs(rateD)>3) = nan; % remove obvious outliers
rateD(abs(rateD)<0) = nan; % This is not a zombie simulation

rateD_rshp = reshape(rateD, 7, []);
% [sum(rateD_rshp)]

if ~isempty(R)
    rateR = (diff(R)./diff(datenum(time-time(1))))./Q(2:end);
    rateR(abs(rateR)>3) = nan;
    
    rateR_rshp = reshape(rateR, 7, []);
    
end
%% Define the time
x = datenum(time(2:end)-time(1));
x1 = x(1):1/24:x(end);

%% Compare the fitted and effective rates

if ~isempty(R)
    figure;
    subplot(221)
    title('Death rate')
    plot(x,rateD,'k*',x1,kappaFun(kappa,x1),'r')
    xlabel('Time (days)')
    ylabel('Death rate (day^{-1})')
    axis tight
    legend('Measured','Fitted','location','northwest')
    
    
    subplot(222)
    title('Recovery rate')
    
    plot(x,rateR,'b*',x1,lambdaFun(lambda,x1),'r')
    axis tight
    set(gcf,'color','w')
    xlabel('Time (days)')
    ylabel('Recovery rate (day^{-1})')
    legend('Measured','Fitted')
    
    % test draw
    subplot(223)
    w = [1:1:length(x(1:7:end))]
    % plot(x(1:7:end), [sum(rateD_rshp)])
    plot(w, [sum(rateD_rshp)]/7)
    xlabel('Time (weeks)')
    ylabel('Death rate (week^{-1})')

    subplot(224)
    plot(w, [sum(rateR_rshp)]/7)
    xlabel('Time (weeks)')
    ylabel('Recovery rate (week^{-1})')
else
    figure;
    plot(x,rateD,'k*',x1,kappaFun(kappa,x1),'r')
    xlabel('Time (days)')
    ylabel('Pseudo-death rate (day^{-1})')
    axis tight
    legend('Measured','Fitted')
end
end

