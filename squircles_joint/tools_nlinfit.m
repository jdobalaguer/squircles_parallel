function [beta,r,J] = tools_nlinfit(X,y,model,beta,options)
%NLINFIT Nonlinear least-squares regression.
%   BETA = NLINFIT(X,Y,MODELFUN,BETA0) estimates the coefficients of a
%   nonlinear regression function, using least squares estimation.  Y is a
%   vector of response (dependent variable) values.  Typically, X is a
%   design matrix of predictor (independent variable) values, with one row
%   for each value in Y and one column for each coefficient.  However, X
%   may be any array that MODELFUN is prepared to accept.  MODELFUN is a
%   function, specified using @, that accepts two arguments, a coefficient
%   vector and the array X, and returns a vector of fitted Y values.  BETA0
%   is a vector containing initial values for the coefficients.
%
%   [BETA,R,J] = NLINFIT(X,Y,MODELFUN,BETA0) returns the fitted
%   coefficients BETA, the residuals R, and the Jacobian J of MODELFUN,
%   evaluated at BETA. You can use these outputs with NLPREDCI to produce
%   confidence intervals for predictions, and with NLPARCI to produce
%   confidence intervals for the estimated coefficients.  sum(R.^2)/(N-P),
%   where [N,P] = size(X), is an estimate of the population variance, and
%   inv(J'*J)*sum(R.^2)/(N-P) is an estimate of the covariance matrix of
%   the estimates in BETA.
%
%   [...] = NLINFIT(X,Y,MODELFUN,BETA0,OPTIONS) specifies control parameters
%   for the algorithm used in NLINFIT.  This argument can be created by a
%   call to STATSET.  Applicable STATSET parameters are:
%
%      'MaxIter'     - Maximum number of iterations allowed.  Defaults to 100.
%      'TolFun'      - Termination tolerance on the residual sum of squares.
%                      Defaults to 1e-8.
%      'TolX'        - Termination tolerance on the estimated coefficients
%                      BETA.  Defaults to 1e-8.
%      'Display'     - Level of display output during estimation.  Choices
%                      are 'off' (the default), 'iter', or 'final'.
%      'DerivStep'   - Relative difference used in finite difference gradient
%                      calculation.  May be a scalar, or the same size as
%                      the parameter vector BETA.  Defaults to EPS^(1/3).
%      'FunValCheck' - Check for invalid values, such as NaN or Inf, from
%                      the objective function [ 'off' | 'on' (default) ].
%
%
%   NLINFIT treats NaNs in Y or MODELFUN(BETA0,X) as missing data, and
%   ignores the corresponding observations.
%
%   Examples:
%
%      Use @ to specify MODELFUN:
%         load reaction;
%         beta = nlinfit(reactants,rate,@mymodel,beta);
%
%      where MYMODEL is a MATLAB function such as:
%         function yhat = mymodel(beta, x)
%         yhat = (beta(1)*x(:,2) - x(:,3)/beta(5)) ./ ...
%                        (1+beta(2)*x(:,1)+beta(3)*x(:,2)+beta(4)*x(:,3));
%   
%   See also NLPARCI, NLPREDCI, NLINTOOL, STATSET.

%   References:
%      [1] Seber, G.A.F, and Wild, C.J. (1989) Nonlinear Regression, Wiley.

if nargin < 4
    error('stats:nlinfit:TooFewInputs','NLINFIT requires four input arguments.');
elseif ~isvector(y)
    error('stats:nlinfit:NonVectorY','Requires a vector second input argument.');
end
if nargin < 5
    options = tools_statset('nlinfit');
else
    options = tools_statset(tools_statset('nlinfit'), options);
end

% Check sizes of the model function's outputs while initializing the fitted
% values, residuals, and SSE at the given starting coefficient values.
model = fcnchk(model);
try
    yfit = model(beta,X);
catch
    [errMsg,errID] = lasterr;
    if isa(model, 'inline')
        error('stats:nlinfit:ModelFunctionError',...
             ['The inline model function generated the following ', ...
              'error:\n%s'], errMsg);
    elseif strcmp('MATLAB:UndefinedFunction', errID) ...
                && ~isempty(strfind(errMsg, func2str(model)))
        error('stats:nlinfit:ModelFunctionNotFound',...
              'The model function ''%s'' was not found.', func2str(model));
    else
        error('stats:nlinfit:ModelFunctionError',...
             ['The model function ''%s'' generated the following ', ...
              'error:\n%s'], func2str(model),errMsg);
    end
end
if ~isequal(size(yfit), size(y))
    error('stats:nlinfit:WrongSizeFunOutput', ...
          'MODELFUN should return a vector of fitted values the same length as Y.');
end

% Set up convergence tolerances from options.
maxiter = options.MaxIter;
betatol = options.TolX;
rtol = options.TolFun;
fdiffstep = options.DerivStep;
funValCheck = strcmp(options.FunValCheck, 'on');

% Find NaNs in either the responses or in the fitted values at the starting
% point.  Since X is allowed to be anything, we can't just check for rows
% with NaNs, so checking yhat is the appropriate thing.  Those positions in
% the fit will be ignored as missing values.  NaNs that show up anywhere
% else during iteration will be treated as bad values.
nans = (isnan(y(:)) | isnan(yfit(:))); % a col vector
n = sum(~nans);
p = numel(beta);

sqrteps = sqrt(eps(class(beta)));
J = zeros(n,p,class(yfit));
r = y(:) - yfit(:);
r(nans) = [];
sse = r'*r;
if funValCheck && ~isfinite(sse), checkFunVals(r); end

zbeta = zeros(size(beta),class(beta));
zerosp = zeros(p,1,class(r));

% Set initial weight for LM algorithm.
lambda = .01;

% Set the level of display
switch options.Display
case 'off',    verbose = 0;
case 'notify', verbose = 1;
case 'final',  verbose = 2;
case 'iter',   verbose = 3;
end

if verbose > 2 % iter
    disp(' ');
    disp('                                     Norm of         Norm of');
    disp('   Iteration             SSE        Gradient           Step ');
    disp('  -----------------------------------------------------------');
    disp(sprintf('      %6d    %12g',0,sse));
end

iter = 0;
breakOut = false;
while iter < maxiter
    iter = iter + 1;
    betaold = beta;
    sseold = sse;

    % Compute a finite difference approximation to the Jacobian
    for k = 1:p
        delta = zbeta;
        if (beta(k) == 0)
            nb = sqrt(norm(beta));
            delta(k) = fdiffstep * (nb + (nb==0));
        else
            delta(k) = fdiffstep*beta(k);
        end
        yplus = model(beta+delta,X);
        dy = yplus(:) - yfit(:);
        dy(nans) = [];
        J(:,k) = dy/delta(k);
    end

    % Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r
    diagJtJ = sum(abs(J).^2);
    if funValCheck && ~all(isfinite(diagJtJ)), checkFunVals(J(:)); end
    Jplus = [J; diag(sqrt(lambda*diagJtJ))];
    rplus = [r; zerosp];
    step = Jplus \ rplus;
    beta(:) = beta(:) + step;

    % Evaluate the fitted values at the new coefficients and
    % compute the residuals and the SSE.
    yfit = model(beta,X);
    r = y(:) - yfit(:);
    r(nans) = [];
    sse = r'*r;
    if funValCheck && ~isfinite(sse), checkFunVals(r); end

    % If the LM step decreased the SSE, decrease lambda to downweight the
    % steepest descent direction.
    if sse < sseold
        lambda = 0.1*lambda;

    % If the LM step increased the SSE, repeatedly increase lambda to
    % upweight the steepest descent direction and decrease the step size
    % until we get a step that does decrease SSE.
    else
        while sse > sseold
            lambda = 10*lambda;
            if lambda > 1e16
                warning('stats:nlinfit:UnableToDecreaseSSE', ...
                        'Unable to find a step that will decrease SSE.  Returning results from last iteration.');
                breakOut = true;
                break
            end
            Jplus = [J; diag(sqrt(lambda*sum(J.^2)))];
            step = Jplus \ rplus;
            beta(:) = betaold(:) + step;
            yfit = model(beta,X);
            r = y(:) - yfit(:);
            r(nans) = [];
            sse = r'*r;
            if funValCheck && ~isfinite(sse), checkFunVals(r); end
        end
    end
    
    if verbose > 2 % iter
        disp(sprintf('      %6d    %12g    %12g    %12g', ...
                     iter,sse,norm(2*r'*J),norm(step)));
    end

    % Check step size and change in SSE for convergence.
    if norm(step) < betatol*(sqrteps+norm(beta))
        if verbose > 1 % 'final' or 'iter'
            disp('Iterations terminated: relative norm of the current step is less than OPTIONS.TolX');
        end
        break
    elseif abs(sse-sseold) <= rtol*sse
        if verbose > 1 % 'final' or 'iter'
            disp('Iterations terminated: relative change in SSE less than OPTIONS.TolFun');
        end
        break
    elseif breakOut
        break
    end
end

if (iter >= maxiter)
    warning('stats:nlinfit:IterationLimitExceeded', ...
            'Iteration limit exceeded.  Returning results from final iteration.');
end

% If the Jacobian is ill-conditioned, then two parameters are probably
% aliased and the estimates will be highly correlated.  Prediction at new x
% values not in the same column space is dubious.  NLPARCI will have
% trouble computing CIs because the inverse of J'*J is difficult to get
% accurately.  NLPREDCI will have the same difficulty, and in addition,
% will in effect end up taking the difference of two very large, but nearly
% equal, variance and covariance terms, lose precision, and so the
% prediction bands will be erratic.
[Q,R] = qr(J,0);
if condest(R) > 1/(eps(class(beta)))^(1/2)
    warning('stats:nlinfit:IllConditionedJacobian', ...
            ['The Jacobian at the solution is ill-conditioned, and some\n' ...
             'model parameters may not be estimated well (they are not ' ...
             'identifiable).\nUse caution in making predictions.']);
end
if nargout > 1
    % Return residuals and Jacobian that have missing values where needed.
    r = y - yfit;
    JJ(~nans,:) = J;
    JJ(nans,:) = NaN;
    J = JJ;
    
    % We could estimate the population variance and the covariance matrix
    % for beta here as
    % mse = sum(abs(r).^2)/(n-p);
    % Rinv = inv(R);
    % Sigma = Rinv*Rinv'*mse;
end


function checkFunVals(v)
if any(~isfinite(v))
    error('stats:nlinfit:NonFiniteFunOutput', ...
          'MODELFUN has returned Inf or NaN values.');
end