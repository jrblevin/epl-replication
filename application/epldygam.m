% ----------------------------------------------------------------------------
% EPLDYGAM       Procedure that estimates the structural parameters
%                of dynamic game of firms' entry/exit
%                using the Efficient Pseudo-Likelihood (EPL) algorithm
%
% by Jason Blevins and Adam Dearing
%
% ----------------------------------------------------------------------------
%
% FORMAT:
% [bepl,varb,llike,iter,err,times,pchoice,v] = epldygam(aobs,zobs,aobs_1,zval,ptranz,theta,v,bdisc,miniter,maxiter,verbose)
%
% INPUTS:
%     aobs    - (nobs x nplayer) matrix with observations of
%               firms' activity decisions (1=active ; 0=no active)
%     zobs    - (nobs x 1) vector with observations of market
%               exogenous characteristics.
%     aobs_1  - (nobs x nplayer) matrix with observations of
%               firms' initial states (1=incumbent; 0=potential entrant)
%     zval    - (numz x 1) vector with values of market characteristics
%     ptranz  - (numz x numz) matrix of transition probabilities
%               of market characteristics.
%     theta   - (kparam x 1) vector of parameters used to initialize the
%               procedure.
%     v       - ((numz*2^nplayer) x 2 x nplayer) matrix of players'
%               choice-specific value functions used to initialize EPL.
%     bdisc   - Discount factor
%     miniter - Minimum number of EPL iterations
%     maxiter - Maximum number of EPL iterations
%     verbose - 0 to suppress output, 1 to enable output
%
% OUTPUTS:
%     bepl    - (kparam x iter) matrix with 1-EPL estimates as rows:
%               (theta_fc_1; theta_fc_2; ... ; theta_fc_n; theta_rs; theta_rn; theta_ec).
%     varb    - (kparam x kparam) matrix with estimated
%               covariance matrices of the converged EPL estimates.
%     llike   - Final log likelihood value
%     iter    - Number of NPL iterations completed.
%     err     - Error flag indicator
%     times   - Vector of computational time (sec.) per iteration
%     pchoice - Choice probabilities
%     v       - Value functions
%
% ----------------------------------------------------------------------------

function [ bepl, varb, llike, iter, err, times, pchoice, v ] = epldygam(aobs, zobs, aobs_1, zval, ptranz, theta, v, bdisc, miniter, maxiter, verbose)
% ---------------
% Some constants
% ---------------
eulerc = 0.5772;
myzero = 1e-16;
nobs = size(aobs, 1);
nplayer = size(aobs, 2);
naction = 2;
numa = 2^nplayer;
numz = size(zval, 1);
numx = numz*numa;
kparam = nplayer + 3;
critconv = (1e-2)*(1/kparam);
bcepl = zeros(kparam,1);
bepl = [];
times = [];
varb = zeros(kparam, kparam);
llike = -1e300;

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

% ------------------------------------------------
% Matrix with observed indexes of state variables
% ------------------------------------------------
indzobs = (repmat(zobs, 1,numz)==repmat(zval', nobs, 1))*[ 1:numz ]';
twop = kron(ones(nobs, 1), 2.^(nplayer-[1:nplayer]));
indobs = sum(aobs_1.*twop, 2);
indobs = (indzobs-1).*(2^nplayer) + indobs + 1;
indobs_orig = indobs;
for i = 2:nplayer
    indobs = [ indobs; indobs_orig + (i-1)*numx ];
end

% -------------
% EPL Algorithm
% -------------
aobs = 1 + reshape(aobs, [nobs*nplayer,1]);

% ----------------------------------------------------------------
% Initialize storage for flow utilities, continuation values, etc.
% ----------------------------------------------------------------
u = zeros(numx, kparam, nplayer, naction);
fv = zeros(numx, nplayer, naction);
pchoice = zeros(numx, nplayer);
pchoice_est = zeros(numx, nplayer);
ssurplus = zeros(numx, nplayer);
pchoice_dv = zeros(numx, naction, naction, nplayer); % Derivative of Lambda(x,a,i) wrt v(x,a',i)

% --------------------------
% Initialize Jacobian matrix
% --------------------------

% Jacobian dimension: sparse (jacdim x jacdim)
jacdim = numx * naction * nplayer;

% Number of nonzeros in Jacobian
% There are (naction * nplayer) x (naction * nplayer) diagonal blocks of dimension (numx x numx).
jacnz = numx*(naction*nplayer)^2 - 4*numx*nplayer + nplayer*4*numx^2;

% Triplets for sparse jacobian
I = zeros(jacnz, 1);
J = zeros(jacnz, 1);
X = zeros(jacnz, 1);

% Calculate initial choice probabilities
for i = 1:nplayer
    vdiff = v(:,i,2) - v(:,i,1);
    expvdiff = exp(vdiff);
    pchoice(:,i) = expvdiff ./ (1 + expvdiff);
end

% Iterate until convergence criterion met
criter = 1000;
iter = 1;
while ((criter > critconv) && (iter <= maxiter))
    if (verbose > 0)
        %fprintf('\n');
        fprintf('-----------------------------------------------------\n');
        fprintf('EPL ESTIMATOR: K = %d\n', iter);
        fprintf('-----------------------------------------------------\n');
        fprintf('\n');
    end
    tic;

    % ------------------------------------------------
    % Matrices of choice probabilities,
    % social surplus, and
    % derivatives of choice probabilities.
    % ------------------------------------------------
    for i = 1:nplayer
        vdiff = v(:,i,2) - v(:,i,1);
        expvdiff = exp(vdiff);
        ssurplus(:,i) = log(expvdiff + 1) + v(:,i,1) + eulerc;

        % Derivative of Lambda_^i(:,a^i=2) with respect to v^i(:,a^i=2)
        pchoice_dv(:,2,2,i) = pchoice(:,i) .* (1 - pchoice(:,i));
        % Derivative of Lambda_^i(:,a^i=2) with respect to v^i(:,a^i=1)
        pchoice_dv(:,2,1,i) = -pchoice_dv(:,2,2,i);
        % Derivative of Lambda_^i(:,a^i=1) with respect to v^i(:,a^i=1)
        pchoice_dv(:,1,1,i) = pchoice_dv(:,2,2,i);
        % Derivative of Lambda_^i(:,a^i=1) with respect to v^i(:,a^i=2)
        pchoice_dv(:,1,2,i) = -pchoice_dv(:,1,1,i);
    end

    % -----------------------------------------------------------
    % Matrix of transition probabilities Pr(a[t]|s[t],a[t-1])
    % -----------------------------------------------------------
    ptrana = ones(numx, numa);
    for i = 1:nplayer;
        mi = aval(:,i)';
        ppi = pchoice(:,i);
        ppi1 = repmat(ppi, 1, numa);
        ppi0 = 1 - ppi1;
        mi1 = repmat(mi, numx, 1);
        mi0 = 1 - mi1;
        ptrana = ptrana .* (ppi1 .^ mi1) .* (ppi0 .^ mi0);
    end

    % Reset Jacobian index counter
    idx = 0;
    for i = 1:nplayer
        % --------------------------------------------
        % Matrices Pr(a[t] | s[t], a[t-1], ai[t])
        % --------------------------------------------
        mi = aval(:,i)';
        ppi = pchoice(:,i);
        ppi= (ppi>=myzero).*(ppi<=(1-myzero)).*ppi ...
             + (ppi<myzero).*myzero ...
             + (ppi>(1-myzero)).*(1-myzero);
        ppi1 = repmat(ppi, 1, numa);
        ppi0 = 1 - ppi1;
        mi1 = repmat(mi, numx, 1);
        mi0 = 1 - mi1;
        ptrani = ((ppi1 .^ mi1) .* (ppi0 .^ mi0));
        ptranai0 = ptrana .* (mi0 ./ ptrani);
        ptranai1 = ptrana .* (mi1 ./ ptrani);

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
        u(:,:,i,1) = zeros(numx, kparam); % u is initialized with zeros
        fv(:,i,1) = bdisc * Fi0 * ssurplus(:,i);      % Player i chooses 0

        % Player i chooses 1
        u1i = zeros(1, nplayer);
        u1i(1,i) = 1; % Dummy for FC for player i
        u(:,:,i,2) = [ repmat(u1i, numx, 1), mstate(:,1), (-hi), (mstate(:,i+1)-1) ];
        fv(:,i,2) = bdisc * Fi1 * ssurplus(:,i);      % Player i chooses 0

        % ------------------------------
        % Calculate the Jacobian
        % ------------------------------
        %
        % The Jacobian is a sparse matrix with the following block structure.
        % Each (ai,aj) cell represents a matrix of dimension (numx x numx).
        %
        % [  (a1=0,a1=0)  (a1=0,a2=0)  ...  (a1=0,aN=0) |  (a1=0,a1=1)  (a1=0,a2=1)  ...  (a1=0,aN=1) ]
        % [  (a2=0,a1=0)  (a2=0,a2=0)  ...  (a2=0,aN=0) |  (a2=0,a1=1)  (a2=0,a2=1)  ...  (a2=0,aN=1) ]
        % [      ...          ...      ...      ...     |      ...          ...      ...      ...     ]
        % [  (aN=0,a1=0)  (aN=0,a2=0)  ...  (aN=0,aN=0) |  (aN=0,a1=1)  (aN=0,a2=1)  ...  (aN=0,aN=1) ]
        % [---------------------------------------------|---------------------------------------------]
        % [  (a1=1,a1=0)  (a1=1,a2=0)  ...  (a1=1,aN=0) |  (a1=1,a1=1)  (a1=1,a2=1)  ...  (a1=1,aN=1) ]
        % [  (a2=1,a1=0)  (a2=1,a2=0)  ...  (a2=1,aN=0) |  (a2=1,a1=1)  (a2=1,a2=1)  ...  (a2=1,aN=1) ]
        % [      ...          ...      ...      ...     |      ...          ...      ...      ...     ]
        % [  (aN=1,a1=0)  (aN=1,a2=0)  ...  (aN=1,aN=0) |  (aN=1,a1=1)  (aN=1,a2=1)  ...  (aN=1,aN=1) ]

        for j = 1:nplayer
            if j == i
                % ------------------------------
                % Diagonal Jacobian blocks
                % ------------------------------

                % First, store diagonal elements of these blocks as vectors (for clarity)
                jac_i0_i0 = eye(numx,numx) - bdisc * Fi0 .* repmat((1-ppi)', numx, 1);
                jac_i1_i1 = eye(numx,numx) - bdisc * Fi1 .* repmat(ppi', numx, 1);
                jac_i0_i1 = - bdisc * Fi0 .* repmat(ppi', numx, 1);
                jac_i1_i0 = - bdisc * Fi1 .* repmat((1-ppi)', numx, 1);

                % Calculate row and column indices
                jac_i0_row = repmat([(i-1)*numx+1:(i-1)*numx+numx]', 1, numx);
                jac_i0_col = repmat([(i-1)*numx+1:(i-1)*numx+numx], numx, 1);
                jac_i1_row = jac_i0_row + jacdim/2;
                jac_i1_col = jac_i0_col + jacdim/2;

                % Store values for sparse Jacobian
                I(idx+1:idx+4*numx^2) = [ jac_i0_row(:); jac_i0_row(:); jac_i1_row(:); jac_i1_row(:); ];
                J(idx+1:idx+4*numx^2) = [ jac_i0_col(:); jac_i1_col(:); jac_i0_col(:); jac_i1_col(:); ];
                X(idx+1:idx+4*numx^2) = [ jac_i0_i0(:);  jac_i0_i1(:);  jac_i1_i0(:);  jac_i1_i1(:);  ];
                idx = idx + 4*numx^2;
            else
                % --------------------------------------------
                % Matrices Pr(a[t] | s[t], a[t-1], ai[t], aj[t])
                % --------------------------------------------
                mj = aval(:,j)';
                ppj = pchoice(:,j);
                ppj= (ppj>=myzero).*(ppj<=(1-myzero)).*ppj ...
                     + (ppj<myzero).*myzero ...
                     + (ppj>(1-myzero)).*(1-myzero);
                ppj1 = repmat(ppj, 1, numa);
                ppj0 = 1 - ppj1;
                mj1 = repmat(mj, numx, 1);
                mj0 = 1 - mj1;
                ptranj = ((ppj1 .^ mj1) .* (ppj0 .^ mj0));

                % --------------------------------------------
                % Derivatives of transition probability Pr(s[t+1] | s[t], a[t-1], ai[t])
                % --------------------------------------------
                ppj1_dvj0 = repmat(pchoice_dv(:,2,1,j), 1, numa);
                ppj1_dvj1 = repmat(pchoice_dv(:,2,2,j), 1, numa);
                ppj0_dvj0 = repmat(pchoice_dv(:,1,1,j), 1, numa);
                ppj0_dvj1 = repmat(pchoice_dv(:,1,2,j), 1, numa);

                ptranj_dvj0 = ((ppj1_dvj0 .^ mj1) .* (ppj0_dvj0 .^ mj0));
                ptranj_dvj1 = ((ppj1_dvj1 .^ mj1) .* (ppj0_dvj1 .^ mj0));

                ptranai0_dvj0 = ptranai0 ./ ptranj .* ptranj_dvj0;
                ptranai0_dvj1 = ptranai0 ./ ptranj .* ptranj_dvj1;
                ptranai1_dvj0 = ptranai1 ./ ptranj .* ptranj_dvj0;
                ptranai1_dvj1 = ptranai1 ./ ptranj .* ptranj_dvj1;

                Fi0_dvj0 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai0_dvj0);
                Fi0_dvj1 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai0_dvj1);
                Fi1_dvj0 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai1_dvj0);
                Fi1_dvj1 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai1_dvj1);

                % ------------------------------
                % Off-diagonal Jacobian blocks
                % ------------------------------

                % Choice ai = 0 (flow utility is zero, future value term only)
                jac_i0_j0_val = -bdisc*Fi0_dvj0*ssurplus(:,i);
                jac_i0_j1_val = -bdisc*Fi0_dvj1*ssurplus(:,i);

                % Used for calculating derivative of expected log number of rivals given ai=1
                hitmp = aval;
                hitmp(:,i) = ones(numa, 1);

                utmp = -ptranai1_dvj0 * log(sum(hitmp, 2)); % Appears in u with a negative sign
                jac_i1_j0_val = -utmp.*bcepl(nplayer+2) - bdisc*Fi1_dvj0*ssurplus(:,i);

                utmp = -ptranai1_dvj1 * log(sum(hitmp, 2)); % Appears in u with a negative sign
                jac_i1_j1_val = -utmp.*bcepl(nplayer+2) - bdisc*Fi1_dvj1*ssurplus(:,i);

                % Calculate row and column offsets
                jac_i0_row = [(i-1)*numx+1:(i-1)*numx+numx]';
                jac_i1_row = jac_i0_row + jacdim/2;
                jac_j0_col = [(j-1)*numx+1:(j-1)*numx+numx]';
                jac_j1_col = jac_j0_col + jacdim/2;

                % Store values for sparse Jacobian
                I(idx+1:idx+4*numx) = [ jac_i0_row;    jac_i0_row;    jac_i1_row;    jac_i1_row;    ];
                J(idx+1:idx+4*numx) = [ jac_j0_col;    jac_j1_col;    jac_j0_col;    jac_j1_col;    ];
                X(idx+1:idx+4*numx) = [ jac_i0_j0_val; jac_i0_j1_val; jac_i1_j0_val; jac_i1_j1_val; ];
                idx = idx + 4*numx;
            end
        end
    end

    % -------------------
    % Combine Jacobian blocks into a single sparse matrix
    % -------------------
    jac = sparse(I, J, X, jacdim, jacdim, jacnz);

    vvec = v(:);
    eqcond = @(v) Gfunc(zval, ptranz, theta, v, bdisc, nplayer, naction);

    % -------------------
    % Reshape utility function (u), future values (fv), and value function (v) to conform
    % -------------------
    vvec = v(:);
    fvvec = fv(:);
    uperm = permute(u, [ 1, 3, 4, 2 ]);
    uvec = reshape(uperm, [numx * nplayer * naction, kparam]);

    % -------------------
    % Creating ww
    % -------------------

    % Solve jac * ww = [ u, fv - v ] for ww.
    rhs = [ uvec, fvvec - vvec ];
    ww = jac \ rhs; % MLDIVIDE uses sparsity.

    % -------------------------------------------
    % Separate vectors by choice a = 0, a = 1
    % -------------------------------------------
    a0ind = [1:jacdim/2]';
    a1ind = a0ind + jacdim/2;
    x0 = ww(a0ind,1:kparam);
    x1 = ww(a1ind,1:kparam);
    r0 = vvec(a0ind) + ww(a0ind,kparam + 1);
    r1 = vvec(a1ind) + ww(a1ind,kparam + 1);

    % -------------------------------------------
    % Creating observations xobs (_x_ vars) and robs (_r_est)
    % -------------------------------------------
    xobs0 = x0(indobs,:);
    xobs1 = x1(indobs,:);
    robs0 = r0(indobs,:);
    robs1 = r1(indobs,:);

    % ------------------------------------------
    % Pseudo-Maximum Likelihood Estimation
    % ------------------------------------------
    [ thetaest, varb, llike, cl_iter, err ] = clogit(aobs, [xobs0, xobs1], [robs0, robs1], 100, bcepl);
    if err > 0
        % Undo final iteration increment
        iter = iter - 1;
        return;
    end

    % ----------------------------
    % Update the value function
    % ----------------------------
    vvec = vvec + ww(:,1:kparam)*thetaest + ww(:,kparam+1);

    % Reshape v for next iteration
    v = reshape(vvec, [numx, nplayer, naction]);

    % Update choice probabilities
    for i = 1:nplayer
        vdiff = v(:,i,2) - v(:,i,1);
        expvdiff = exp(vdiff);
        pchoice_est(:,i) = expvdiff ./ (1 + expvdiff);
    end

    % Check for convergence after minimum number of iterations
    if (iter > miniter)
        criter1 = max(abs(thetaest - bcepl));
        criter2 = max(max((abs(pchoice_est - pchoice))));
        criter = max(criter1, criter2);
    end

    % Save values for next iteration
    bcepl = thetaest;
    pchoice = pchoice_est;

    % Collect estimates
    bepl = [ bepl; thetaest' ];

    % Collect times
    times = [ times; toc ];

    if (verbose > 0)
        disp('theta = ')
        disp(bcepl')
        disp('llike = ')
        disp(llike)
    end

    % Proceed to the next iteration
    iter = iter + 1;
end

% Undo final iteration increment
iter = iter - 1;

end
