function [x,f, exitflag,costdata] = ism_steepDesc(f,x0,tol,maxit)
%
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,histout,costdata] = steep(x0,f,tol,maxit)
%
% steepest descent, simple linesearch 
% 
%
% Input: x0 = initial iterate
%        f = objective function,
%            the calling sequence for f should be
%            [fout,gout]=f(x) where fout=f(x) is a scalar
%              and gout = grad f(x) is a COLUMN vector
%        tol = termination criterion norm(grad) < tol
%              optional, default = 1.d-6
%        maxit = maximum iterations (optional) default = 1000
%
% Output: x = solution
%         histout = iteration history   
%             Each row of histout is
%            [norm(grad), f, number of step length reductions, iteration count]
%         costdata = [num f, num grad, num hess] (for steep, num hess=0)
%

%
% linesearch parms

if nargin < 4                   %Maximum number of iterations
maxit=15; 
end
if nargin < 3                   %Tolerance, as fraction (xc - xn)/xc
tol=0.1;
end

itc=1; xc=x0;                   %Initialize                    
[fc,gc]=feval(f,xc);            %cost,gradient
costdata = zeros(maxit,1);
costdata(1) = fc;
chng = Inf;
exitflag = 0;


tau = 0.5;                     %Line stepping parameters
lambda = 10;

disp(['Initial Cost: ', num2str(fc)])


while(chng > tol && itc <= maxit)
    disp(['Iteration: ', num2str(itc)])    
    xn=xc-lambda*gc;
    fn = feval(f,xn); 
    
    iarm = 1;
    fp = Inf;
    
	while(((fn > fc) && (fp > fc)) || fn < fp)                      %line search
		fp = fn;
        xp = xn;
        
        lambda = lambda*tau;
		xn=xc-lambda*gc;
		fn=feval(f,xn); 
		
        if (iarm > 15) 
		disp('Error in steepest descent ')
        x = xc;
        f = fn;
        exitflag = 1;
		return; end
                
    iarm=iarm+1;
    end
    
    lambda = lambda*(1/tau);
    
    chng = (fc-fp)/fc;          %Iteration information
    itc = itc + 1;    
    
    disp(['Updated Cost: ', num2str(fn)])
    disp(['% Change: ', num2str(100*chng)])
    
	xc=xp; 
    fc = fp;
    if chng > tol, [fc,gc]=feval(f,xc); end; %Update Values
    costdata(itc) = fc;
end

x=xc;
f= fc;


end
