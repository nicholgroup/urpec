function [beta1,r,j,COVB,mse,err,se] = fitwrap(ctrl, x, y, beta0, model, mask)
%function [beta1,r,j,COVB,mse,err] = fitwrap(ctrl, x, y, beta0, model, mask)
% beta1 = fitwrap(ctrl, x, y, beta0, model, mask)
% ctrl: plinit, plfit, woff, nofit, pause, samefig, fine, robust
% y = model(beta, x);     
%plinit plots the initial guess
%plfit plots the fit
%woff turns off warnings
%no fit does not fit the data
%pause will pause after each plot
%samefig uses current figure. not 500.
%fine: fitting option
%robust: fitting option
% beta0 can be a vector, array, cell array with initialization function handles 
% or initial value vectors, or a struct array with fields fn and args, and optionally vals.
% Vinite vals entries override the return of the initialization function.
% Initialization functions are called as beta0 = fn(x, y, args).
% data should be in columns.
% beta1 are the fitted parameters
% r are the residuals
% j is the jacobian
% COVB is the variance-covariance matrix
% mse contains details about the error model
% err are the upper and lower error bounds for the predictions. These are
% 95% confidence bounds.
% SE are the standard error bounds
% is val=beta(i,j) - err(i,j,1) + err(i,j,1)
% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.


% This allows the model to be something like '@(p,x) x*p(1)', which
% overcomes issues with saving functions to disk.
if ischar(model)
    model=str2func(model);
end

n = size(y, 1);
if size(y,2) == 1
    fprintf('X is %d x %d, Y is %d x %d\n',size(x,1),size(x,2),size(y,1),size(y,2));
    error('It is unlikely you wanted to fit a single data point.  Try transposing Y');
end

if size(x,2) ~= size(y,2)
    fprintf('X is %d x %d, Y is %d x %d\n',size(x,1),size(x,2),size(y,1),size(y,2));
    warning('This is probably not what you want.');
end

if size(x, 1) == 1
    x = repmat(x, n, 1);
end

if isa(beta0, 'function_handle')
    beta0 = {beta0};
end

if isreal(beta0) && size(beta0, 1) == 1 || length(beta0) == 1
    beta0 = repmat(beta0, n, 1);
end

if nargin >= 6 
    mask = logical(mask);
end

if strfind(ctrl, 'fine')
    options=statset('TolX',1e-20,'TolFun',1e-20);
else
    options=statset();
end
if strfind(ctrl, 'robust')
    options=statset(options,'Robust','on');
end
if strfind(ctrl, 'woff')
    ws(1) = warning('query', 'stats:nlinfit:IllConditionedJacobian');
    ws(2) = warning('query', 'stats:nlinfit:IterationLimitExceeded');
    ws2 = ws;
    [ws2.state] = deal('off');
    warning(ws2);
end

for i = 1:n
    if ~isempty(strfind(ctrl, 'pl')) && isempty(strfind(ctrl,'samefig'))
        figure(500);
        clf;
        hold on;
    end

    if iscell(beta0) 
        if isreal(beta0{i})
            beta2 = beta0{i};
        else 
            beta2 = beta0{i}(x(i, :), y(i, :));
        end
    elseif isstruct(beta0)      
        if ~isempty(beta0(i).fn)
            beta2 = beta0(i).fn(x(i, :), y(i, :), beta0(i).args{:});
            if isfield(beta0, 'vals');
                beta2(isfinite(beta0(i).vals)) = beta0(i).vals(isfinite(beta0(i).vals));
            end
        else
            beta2 = beta0(i).vals;
        end                    
    else
        beta2 = beta0(i, :);
    end
    
    if i == 1
        nfp = length(beta2);
        
        if nargin < 6 || isempty(mask)
            mask = true(1, nfp);
        end
        beta1 = zeros(n, nfp);
    end
    
    beta1(i, :) = beta2;
        
    if ~isempty(strfind(ctrl, 'plinit'))
        plot(x(i, :), y(i, :), '.-', x(i, :), model(beta1(i, :), x(i, :)), 'r--');
    end
    
    if isempty(strfind(ctrl, 'nofit'))
      [betaT,rT,JT,COVBT,mseT] = nlinfit(x(i, :), y(i, :), @fitfn, beta1(i, mask),options);
      beta1(i, mask) = betaT;
      r(i,:)=rT(:);
      j(i,:)=JT(:);
      COVB(i,:)=COVBT(:);
      mse(i,:)=mseT(:);
      try
          %compute standard error JMN 2020/03/04
          alpha=.05; %95 confidence interval
          ci = nlparci(betaT,rT,JT,alpha);
          t = tinv(1-alpha/2,length(x(i,:))-length(betaT));
          nlinfit_se = (ci(:,2)-ci(:,1)) ./ (2*t);
          se(i,mask)=nlinfit_se;
          
          %compute 95% confidence interval
          err=nlparci(betaT,rT,'Jacobian',JT);
          err(i,mask,1:2)=abs(repmat(beta1(i,mask)',1,2)-nlparci(betaT,rT,'Jacobian',JT));
      catch
         warning('could not propagate errors properly')
         err(i,mask,1:2) = nan(sum(mask),2);
         se(i,mask,1:2) = nan(sum(mask),1);
      end
    end
     
    if ~isempty(strfind(ctrl, 'plfit'))
        plot(x(i, :), y(i, :), '.-', x(i, :), model(beta1(i, :), x(i, :)), 'k');
    end
    
    if ~isempty(strfind(ctrl, 'pause')) && (i < n)
        pause
    end
    
    if ~isempty(strfind(ctrl,'resid'))
      f=gcf;
      figure(501);
      if isempty(strfind(ctrl,'samefig'))
        clf;
      else
        hold on;
      end
      plot(x(i,:),y(i,:)-model(beta1(i,:),x(i,:)),'rx-');
      figure(f);
    end
end

   
if strfind(ctrl, 'woff')
    warning(ws);
end


    function y = fitfn(beta, x)
        beta([find(mask), find(~mask)]) = [beta, beta2(~mask)];
        y = model(beta, x);        
    end

end


