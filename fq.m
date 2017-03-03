function I = fq(x)
global t T N a b sigma eta rho cashflow tenor x0 y0
ti=tenor;
NT1 = N(1);
psi=@(t) 0.072;
ipsi=@(t,T) psi(t).*(T-t);
MxT=@(s,t) (sigma.^2./a.^2  + rho.*sigma.*eta./(a.*b)) .*(1 - exp(-a.*(t-s))) - (sigma.^2)./(2.*a.^2).*(exp(-a.*(T-t)) - exp(-a.*(T+t-2.*s)))-rho.*sigma.*eta./(b.*(a+b)).*...
    (exp(-b.*(T -t)) - exp(-b.*T -a.*t+(a+b).*s));
MyT=@(s,t) (eta.^2./b.^2  + rho.*sigma.*eta./(a.*b)) .*(1 - exp(-b.*(t-s))) - eta.^2./(2.*b.^2).*(exp(-b.*(T-t)) - exp(-b.*(T+t-2.*s)))-rho.*sigma.*eta./(a.*(a+b)).*...
    (exp(-a.*(T -t)) - exp(-a.*T -b.*t+(a+b).*s));
B=@(z,t,T) (1 - exp(-z.*(T-t)))./z;
V = @(t,T) (sigma.^2./a.^2).*((T-t) + (2./a).*exp(-a.*(T-t)) - (1./(2.*a)).*exp(-2.*a.*(T-t))- 3./(2.*a))...
    +(eta.^2 ./ b.^2).*((T-t) +   (2./b).*exp(-b.*(T-t)) - (1./(2.*b)).*exp(-2.*b.*(T-t))   -3./(2.*b))...
    +(2.*rho.*sigma.*eta./(a.*b)).*((T-t) +(1./a).*(exp(-a.*(T-t)) -1) + (1./b).*(exp(-b.*(T-t))-1) -(1./(a+b)).*(exp(-(a+b).*(T-t)) -1));
A =@(t,T) exp(0.5.*V(t,T) - ipsi(t,T));
P2 = @(t,T,x,y) repmat(A(t,T),length(x),1).*exp(-repmat(B(a,t,T),length(x),1).*repmat(x,1,length(T)) - repmat(B(b,t,T),length(y),1).*repmat(y,1,length(T)));



mux = -MxT(t,T);
muy = -MyT(t,T);
sigmax = sigma.*sqrt((1 - exp(-2.*a.*(T-t)))./(2.*a));
sigmay = eta.*sqrt((1 - exp(-2.*b.*(T-t)))./(2.*b));
rhoxy = rho.*sigma.*eta./((a+b).*sigmax.*sigmay) .*(1 - exp(-(a+b).*(T-t)));
h1=@(x,ybar) rhoxy.*(x - mux)./(sigmax.*sqrt(1 - rhoxy.^2)) - (ybar - muy)./(sigmay.*sqrt(1 - rhoxy.^2));
h2=@(x,ybar) h1(x,ybar) - B(b,T,ti).*sigmay.*sqrt(1 - rhoxy.^2);
lambdai=@(x) cashflow.*A(T,ti).*exp(-B(a,T,ti).*x);
kappai=@(x) -B(b,T,ti).*(muy - 0.5.*(1 - rhoxy.^2).*sigmay.^2 .* B(b,T,ti) + rhoxy.*sigmay.*((x - mux)./sigmax));

func = @(y) sum(cashflow.*A(T,ti).*exp(-B(a,T,ti).*x - B(b,T,ti).*y)) - 1; %solve by fzero
ybar = fzero(func,0);
I = P2(t,T,x0,y0).*((exp(-0.5.*((x-mux)./sigmax).^2)./(sigmax.*sqrt(2.*pi))) .* (sum(lambdai(x).*exp(kappai(x)).*normcdf(h2(x,ybar)),2)-...
        NT1.*normcdf(h1(x,ybar))));
end