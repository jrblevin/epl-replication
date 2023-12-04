% Gfunc.m -- Evaluates equilibrium conditions for dynamic entry/exit model
%
% FORMAT:
% [result] = Gfunc(zval, ptranz, theta, v, bdisc, nplayer, naction)
%
% INPUTS:
%     zval    - (numz x 1) vector with values of market characteristics
%     ptranz  - (numz x numz) matrix of transition probabilities
%               of market characteristics.
%     theta   - (kparam x 1) vector of structural parameters
%     v       - (numx*nplayer*naction x 1) value function vector
%     bdisc   - discount factor
%     nplayer - number of players
%     naction - number of actions
%
% OUTPUT:
%
%     result  - (numx*nplayer*naction x 1) vector of equilibrium condition
%               values (i.e., the residual G(v) = v - \Gamma(v) of the fixed
%               point mapping evaluated at the input v).

function [ result ] = Gfunc(zval, ptranz, theta, v, bdisc, nplayer, naction)

% ---------------
% Some constants
% ---------------
eulerc = 0.5772;
numa = 2^nplayer;
numz = size(zval, 1);
numx = numz * numa;
kparam = nplayer + 3;
sigmaeps = 1.0;

% fprintf('\n\n');
% fprintf('*****************************************************************************************\n');
% fprintf('   EVALUATING EPL EQUILIBRIUM CONDITIONS\n');
% fprintf('*****************************************************************************************\n');
% fprintf('\n');
% fprintf('----------------------------------------------------------------------------------------\n');
% fprintf('       Values of the structural parameters\n');
% fprintf('\n');
% for i = 1:nplayer
%     fprintf('                       Fixed cost firm %d   = %12.4f\n', i, theta(i));
% end
% fprintf('       Parameter of market size (theta_rs) = %12.4f\n', theta(nplayer+1));
% fprintf('Parameter of competition effect (theta_rn) = %12.4f\n', theta(nplayer+2));
% fprintf('                     Entry cost (theta_ec) = %12.4f\n', theta(nplayer+3));
% fprintf('                       Discount factor     = %12.4f\n', bdisc);
% fprintf('                    Std. Dev. epsilons     = %12.4f\n', 1.0);
% fprintf('\n');
% fprintf('----------------------------------------------------------------------------------------\n');
% fprintf('\n');

% ----------------------------------------------
% Matrix with values of states of (s[t],a[t-1])
% ----------------------------------------------
aval = zeros(numa, nplayer);
for i = 1:nplayer
    aval(:,i) = kron(ones(2^(i-1),1), kron([ 0; 1 ], ones(2^(nplayer-i),1)));
end
mstate = zeros(numx, nplayer+1);
mstate(:,1) = kron(zval, ones(numa, 1));
mstate(:,2:nplayer+1) = kron(ones(numz, 1), aval);

% Reshape vector v into a multi-dimensional array for clarity
v = reshape(v, [numx, nplayer, naction]);

% ------------------------------------------------
% Matrices of choice probabilities,
% social surplus, and
% derivatives of choice probabilities.
% ------------------------------------------------
pchoice = zeros(numx, nplayer);
ssurplus = zeros(numx, nplayer);
for i = 1:nplayer
    vdiff = v(:,i,2) - v(:,i,1);
    expvdiff = exp(vdiff);
    pchoice(:,i) = expvdiff ./ (1 + expvdiff);
    ssurplus(:,i) = log(expvdiff + 1) + v(:,i,1) + eulerc;
end

% ----------------------
% Equilibrium conditions
% ----------------------
u = zeros(numx, kparam, nplayer, naction);
fv = zeros(numx, nplayer, naction);
G = zeros(numx, nplayer, naction);

% -----------------------------------------------------------
% Matrix of transition probabilities Pr(a[t]|s[t],a[t-1])
% -----------------------------------------------------------
ptrana = ones(numx, numa);
for i = 1:nplayer
    mi = aval(:,i)';
    ppi = pchoice(:,i);
    ppi1 = repmat(ppi, 1, numa);
    ppi0 = 1 - ppi1;
    mi1 = repmat(mi, numx, 1);
    mi0 = 1 - mi1;
    ptrana = ptrana .* (ppi1 .^ mi1) .* (ppi0 .^ mi0);
end

for i = 1:nplayer
    % --------------------------------------------
    % Matrices Pr(a[t] | s[t], a[t-1], ai[t])
    % --------------------------------------------
    mi = aval(:,i)';
    ppi = pchoice(:,i);
    ppi1 = repmat(ppi, 1, numa);
    ppi0 = 1 - ppi1;
    mi1 = repmat(mi, numx, 1);
    mi0 = 1 - mi1;
    ptrani = ((ppi1 .^ mi1) .* (ppi0 .^ mi0));
    ptranai0 = ptrana .* (mi0 ./ ptrani);
    ptranai1 = ptrana .* (mi1 ./ ptrani);
    clear mi;

    % --------------------------------------------
    % Transition probability matrix Pr(s[t+1] | s[t], a[t-1], ai[t])
    % --------------------------------------------
    Fi0 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai0);
    Fi1 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai1);

    % ------------------------------------
    % Computing hi = E(ln(Sum(aj)+1))
    % ------------------------------------
    hi = aval;
    hi(:,i) = ones(numa, 1);
    hi = ptranai1 * log(sum(hi, 2));

    % ---------------------------
    % Creating U and FV
    % ---------------------------

    % Player i chooses 0
    %u(:,:,i,1) = zeros(numx, kparam); % u is initialized with zeros
    fv(:,i,1) = bdisc * Fi0 * ssurplus(:,i);      % Player i chooses 0

    % Player i chooses 1
    u1i = zeros(1, nplayer);
    u1i(1,i) = 1; % Dummy for FC for player i
    u(:,:,i,2) = [ repmat(u1i, numx, 1), mstate(:,1), (-hi), (mstate(:,i+1)-1) ];
    fv(:,i,2) = bdisc * Fi1 * ssurplus(:,i);      % Player i chooses 0

    % Equilibrium conditions for player i
    G(:,i,1) = v(:,i,1) - (u(:,:,i,1) * theta + bdisc * Fi0 * ssurplus(:,i));
    G(:,i,2) = v(:,i,2) - (u(:,:,i,2) * theta + bdisc * Fi1 * ssurplus(:,i));
end

% Reshape equilibrium conditions G as vector and return
result = G(:);

% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Profit functions (a = 1):')
% profits = [ u(:,:,1,2) * theta, u(:,:,2,2) * theta, u(:,:,3,2) * theta, u(:,:,4,2) * theta, u(:,:,5,2) * theta ];
% disp(profits)
% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Choice probabilities (a = 1):')
% disp(pchoice)
% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Social surplus functions:')
% disp(ssurplus)
% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Norm of equilibrium conditions:')
% disp(norm(result))
end

