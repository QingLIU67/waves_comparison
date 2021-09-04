function [S,E,I,Q,R,D,P] = SEIQRDP(alpha,beta,gamma,delta,lambda0,kappa0,tau,time_delay,NPI_matrix,NPIs,Npop,E0,I0,Q0,R0,D0,C0,t,lambdaFun,kappaFun) % 1st wave
%function [S,E,I,Q,R,D,P] = SEIQRDP(alpha,beta,gamma,delta,lambda0,kappa0,tau,time_delay,NPI_matrix,Npop,E0,I0,Q0,R0,D0,C0,t,lambdaFun,kappaFun) % 2nd wave
% [S,E,I,Q,R,D,P] = SEIQRDP(alpha,beta,gamma,delta,lambda,kappa,Npop,E0,I0,R0,D0,t,lambdaFun)
% simulate the time-histories of an epidemic outbreak using a generalized
% SEIR model.
%
% Input
%
%   alpha: scalar [1x1]: fitted protection rate
%   beta: scalar [1x1]: fitted  infection rate
%   gamma: scalar [1x1]: fitted  Inverse of the average latent time
%   delta: scalar [1x1]: fitted  rate at which people enter in quarantine
%   lambda: scalar [1x1]: fitted  cure rate
%   kappa: scalar [1x1]: fitted  mortality rate
%   Npop: scalar: Total population of the sample
%   E0: scalar [1x1]: Initial number of exposed cases
%   I0: scalar [1x1]: Initial number of infectious cases
%   Q0: scalar [1x1]: Initial number of quarantined cases
%   R0: scalar [1x1]: Initial number of recovered cases
%   D0: scalar [1x1]: Initial number of dead cases
%   t: vector [1xN] of time (double; it cannot be a datetime)
%   lambdaFun: anonymous function giving the time-dependant recovery rate
%   kappaFun: anonymous function giving the time-dependant death rate
% 
% Output
%   S: vector [1xN] of the target time-histories of the susceptible cases
%   E: vector [1xN] of the target time-histories of the exposed cases
%   I: vector [1xN] of the target time-histories of the infectious cases
%   Q: vector [1xN] of the target time-histories of the quarantined cases
%   R: vector [1xN] of the target time-histories of the recovered cases
%   D: vector [1xN] of the target time-histories of the dead cases
%   P: vector [1xN] of the target time-histories of the insusceptible cases
%
% Author: E. Cheynet - UiB - last modified 27-04-2020
%
% see also fit_SEIQRDP.m

%% Initial conditions
N = numel(t);
Y = zeros(7,N);
Y(1,1) = Npop-Q0-E0-R0-D0-I0-C0;
Y(2,1) = E0;
Y(3,1) = I0;
Y(4,1) = Q0;
Y(5,1) = R0;
Y(6,1) = D0;
Y(7,1) = C0;

if round(sum(Y(:,1))-Npop)~=0
    error(['the sum must be zero because the total population',...
        ' (including the deads) is assumed constant']);
end
%% Computes the seven states
modelFun = @(Y,A,F) A*Y + F;
dt = median(diff(t));

lambda = lambdaFun(lambda0,t);
kappa = kappaFun(kappa0,t);

% ODE resolution
time_delay = round(time_delay);
NPI_ind = 1;
for ii=1:N-1
    % update alpha and tau with NPIs for 1st wave
    if ii > time_delay && ii < 92
        if NPIs(ii-time_delay, 1) == 1
            alpha = alpha + NPI_matrix(NPI_ind);
            NPI_ind = NPI_ind + 1;
        elseif NPIs(ii-time_delay, 1) == 2
            tau = tau + NPI_matrix(NPI_ind);
            NPI_ind = NPI_ind + 1;
        end
    end
    
    A = getA(alpha,gamma,delta,lambda(ii),kappa(ii),tau);
    SI = Y(1,ii)*Y(3,ii);
    F = zeros(7,1);
    F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
    Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
end

% Y = round(Y);
%% Write the outputs
S = Y(1,1:N);
E = Y(2,1:N);
I = Y(3,1:N);
Q = Y(4,1:N);
R = Y(5,1:N);
D = Y(6,1:N);
P = Y(7,1:N);


%% Nested functions
    function [A] = getA(alpha,gamma,delta,lambda,kappa,tau)
        %  [A] = getA(alpha,gamma,delta,lambda,kappa) computes the matrix A
        %  that is found in: dY/dt = A*Y + F
        %
        %   Inputs:
        %   alpha: scalar [1x1]: protection rate
        %   beta: scalar [1x1]: infection rate
        %   gamma: scalar [1x1]: Inverse of the average latent time
        %   delta: scalar [1x1]: rate of people entering in quarantine
        %   lambda: scalar [1x1]: cure rate
        %   kappa: scalar [1x1]: mortality rate
        %   Output:
        %   A: matrix: [7x7]
        A = zeros(7);
        % S
        A(1,1) = -alpha;
        A(1,7) = tau;
        % E
        A(2,2) = -gamma;
        % I
        A(3,2:3) = [gamma,-delta];
        % Q
        A(4,3:4) = [delta,-kappa-lambda];
        % R
        A(5,4) = lambda;
        % D
        A(6,4) = kappa;
        % P
        A(7,1) = alpha;
        A(7,7) = -tau;
    end
    function [Y] = RK4(Fun,Y,A,F,dt)
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*dt,A,F);
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end
end


