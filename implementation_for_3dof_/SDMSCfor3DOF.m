function xk = SDMSCfor3DOF(rd, x0)
% the function is the implementation of SDMSC algorithm for 3oF model
% inputs: rd:---is the computed Lissajous target
%	  x0:---is the angular vector at a particular time step, t_i
%
% output: xk:---is solution obtained by the SDMSC algorithm from solving equ. 4.2

% initializations
xk = x0;
gam=1e-5;
rho = 0.5;
Q0=1;
iter = 0; tmin = 1e-16;

rk = [cos(xk(1)) + cos(xk(1) + xk(2)) + cos(xk(1) + xk(2) + xk(3)); sin(xk(1)) + sin(xk(1) + xk(2)) + sin(xk(1) + xk(2) + xk(3))] - rd;
fkz = 1/2*norm(rk,2)^2;
C0 = fkz;
% Initial Jacobian matrix
Jk = [-sin(xk(1))-sin(xk(1)+xk(2))-sin(xk(1)+xk(2)+xk(3)) -sin(xk(1)+xk(2))-sin(xk(1)+xk(2)+xk(3)) -sin(xk(1)+xk(2)+xk(3)); cos(xk(1))+cos(xk(1)+xk(2))+cos(xk(1)+xk(2)+xk(3)) cos(xk(1)+xk(2))+cos(xk(1)+xk(2)+xk(3)) cos(xk(1)+xk(2)+xk(3))];
% Gradient function gk = J'F;
gk = Jk'*rk;
% Initiate b
diag = ones(length(x0), 1);
%Initial Direction
dk = -diag.*gk;
normg = norm(gk); %norm(gk, inf);
while (normg > 10^(-6)&&iter<=1000)      
    tk = 1;
    xkn = xk + tk*dk;
    rkn = [cos(xkn(1)) + cos(xkn(1) + xkn(2)) + cos(xkn(1) + xkn(2) + xkn(3)); sin(xkn(1)) + sin(xkn(1) + xkn(2)) + sin(xkn(1) + xkn(2) + xkn(3))] - rd;
    fkzn = 0.5*norm(rkn,2)^2;
    dd = dk'*gk;
    while(fkzn > C0 + gam*tk*(dd) && tk > tmin)
%             while (fkzn > fkz + sigma*tk*gk'*dk)
%         if tk <= 0.01
%            tk = rho*tk; % simple backtrack
        if (iter>=5 || iter<=15) 
           t1=rho*tk;
        else
          % quadratic interpolation backtrack only
            t1=-((tk)^2*(dd))/(2*(fkzn-C0-tk*(dd)));
         % safeguard
%             if (t1 < 0.1*tk || t1 > 0.9*tk)
%                t1 = tk/2;
%             end
        end
     tk = t1;
     xkn = xk + tk*dk;
     rkn = [cos(xkn(1)) + cos(xkn(1) + xkn(2)) + cos(xkn(1) + xkn(2) + xkn(3)); sin(xkn(1)) + sin(xkn(1) + xkn(2)) + sin(xkn(1) + xkn(2) + xkn(3))] - rd;
     fkzn = 0.5*(rkn'*rkn); 
    end
    xkn = xk + tk*dk;
    Jkn = [-sin(xkn(1))-sin(xkn(1)+xkn(2))-sin(xkn(1)+xkn(2)+xkn(3)) -sin(xkn(1)+xkn(2))-sin(xkn(1)+xkn(2)+xkn(3)) -sin(xkn(1)+xkn(2)+xkn(3)); cos(xkn(1))+cos(xkn(1)+xkn(2))+cos(xkn(1)+xkn(2)+xkn(3)) cos(xkn(1)+xkn(2))+cos(xkn(1)+xkn(2)+xkn(3)) cos(xkn(1)+xkn(2)+xkn(3))];
    rkn = [cos(xkn(1)) + cos(xkn(1) + xkn(2)) + cos(xkn(1) + xkn(2) + xkn(3)); sin(xkn(1)) + sin(xkn(1) + xkn(2)) + sin(xkn(1) + xkn(2) + xkn(3))] - rd;
%     xkn = xk + alph*dir;
    sk = xkn - xk;
    gkn = Jkn'*rkn;
    c1 = Jkn'*Jkn*sk;
    %   gradient  d1=Jk^Trk
    d1  = Jkn'*rkn;  
%   gradient  d2=J_{k-1}rk
    d2 = Jk'*rkn;
    %d=gt-d2; % d3 = ybar = Jk^Trk - J_{k-1}^Trk
    d3  = d1 - d2;
 
    beta1 = c1+ d3; % short beta
    % the entries of the correction matrix is computed using this formulation
    cdiag = (beta1 - (diag.*sk))./(sk); % newly written formulation (06/05/2022)
    tol1 = 1e-4;
    diagk = diag+cdiag;
%     if diagk >= tol1
%     dk = -gkn./diagk;
%     else
%     dk = -gkn;
%     end
     % updating the values of xc, gc, f's, Q & C
    eta = 1e+30;
    diag = min( max(diagk, tol1), eta);
    dk = -gkn./diag;
    fkz = fkzn;
%     etak= 0.75*exp(-((iter/45))^2)+0.1; % nonmonotone line search
    etak=0; %%Armijo linesearch
    Q1=etak*Q0+1;
    C1=(etak*Q0*C0+fkz)/Q1;
    C0=C1;
    xk = xkn;
    gk = gkn;
    iter = iter +1;  
    normg=norm(gk); %norm(gk,inf);
end
