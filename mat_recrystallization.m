function [sig, C_T, para, his2, tpvpf, plotv, OutputIntegerVariables] = ...
    mat_recrystallization  (eps, para, T, Tdot, his1, time, FEA_FILE, tpvpf, IntegerVariables)

% model for dynamic recrystallization

% this subroutine calculates nonlinear hardening during dynamic
% recrystallization based on the continuum model proposed by brown.

% it accounts for single peak and multiple peak dynamic recrystallization
% for (nrx) continuous cycles of recrystallization

% input:                strain (eps), temperature (T), material parameters
% output:               stress (sig)
% internal variables:   hardening variables/cycle of rx (kappa),
%                       misorientation variables/cycle of rx (zeta),
%                       averaged hardening variable (kappaG)

% - radial return solves for norm_sig_dev
% - the nonlinear evolutions of the internal variables are updated with an
%   explicit algorithm
%
% this routine can be considered as semi implicit

isw = IntegerVariables(1);
fid = IntegerVariables(5);
fids = IntegerVariables(6);
iptyp = IntegerVariables(9);
nstr = IntegerVariables(3);
dt = IntegerVariables(8);

switch isw
    
    case 1
        para = str2double(strsplit(fgetl(fid)));
        npar = size(para,2);
        
        % Initialization
        sig = 0; % why is it initialized as a scalar
        C_T = 0; % why is it initialized as a scalar
        plotv = 0; % why is it initialized as a scalar
        nhis = 40; % nstr + 1; has to work dynamically in the end
        his2 = zeros(nhis,1);
        
        % define plot variables
        [tpvpf,npvpf5] = mat_recrystallization_pp (iptyp);
        OutputIntegerVariables = horzcat(npar, nhis, npvpf5);
        
    case 2
                
        % todo
        % history fields with cell arrays
        
        %% material parameters
        nu = para(4);
        c1 = para(43);
        c2 = para(44);
        c3 = para(45);
        c4 = para(46);
        c5 = para(47);
        c6 = para(48);
        c7 = para(49);
        c8 = para(50);
        r = para(51);
        hc = para(52);
        ctheta = para(53);
        B = para(54);
        m = para(55);
        ckappa = para(56);
        czeta = para(57);
        a = para(58);
        b = para(59);
        c = para(60);
        nrx = para(61);
        
        shear  = 5.47e4;
        dshear = -34.1;
        mu = shear + dshear*T;
        dmu = dshear;
        
        K0 = (2*mu*(1+nu))/(3*(1-2*nu));
        
        %% numerical parameters
        nmax = 200;               % maximum number of iterations
        tol = 1e-8;               % tolerance for newtons scheme
        
        %% history variables out
        epspn   = his1(1:nstr,1);
        Xn      = his1(nstr+1:nstr+(nrx+2),1);
        kappan  = his1(nstr+(nrx+2)+1:nstr+(nrx+2)+(nrx+1),1);
        zetan   = his1(nstr+(nrx+2)+(nrx+1)+1:nstr+(nrx+2)+2*(nrx+1),1);
        kappaGn = his1(nstr+(nrx+2)+2*(nrx+1)+1,1);
        epspvn  = his1(nstr+(nrx+2)+2*(nrx+1)+2,1);
        
        npar   = IntegerVariables(4);
        nhis   = IntegerVariables(7);
        npvpf5 = IntegerVariables(10);
        
        %% initializing variables
        error_local = zeros(nmax,1);
        X      = zeros(nrx+2,1);
        dX     = zeros(nrx+2,1);
        dkappa = zeros(nrx+1,1);
        kappa  = zeros(nrx+1,1);
        dzeta  = zeros(nrx+1,1);
        zeta   = zeros(nrx+1,1);
        kappaX = zeros(nrx+1,1);
        
        %% initial values
        kappa(:,1) = 1e-6;
        zeta(:,1) = 1e-6;
        kappan(:,1) = max(kappan(:,1),1e-6);
        zetan(:,1) = max(zetan(:,1),1e-6);
        
        X(1,1) = 1;
        
        if Xn(1) == 0
            Xn(1) = 1;
            Xn(2) = 1e-6;
            for irx = 2:nrx+1
                Xn(irx+1,1) = 1e-6*Xn(irx,1);
            end
        end
        
        for irx = 1:nrx+1
            X(irx+1,1) = 1e-6*X(irx,1);
        end
        
        %% preliminary preparation
        id2 = [1 1 1 0 0 0]';
        id4 = eye(6);        
        id4_vol = id2*id2'/3;
        id4_dev = id4-id4_vol;
        
        fth = c1*exp(-c2/T);
        nth = c3+c4/T;
        Yth = c5*mu;
        
        %% calculation of plastic evolution
        eps_dev = eps -(eps(1)+eps(2)+eps(3))/3*id2;
        sig_tr_dev = 2*mu*(eps_dev-epspn);
        norm_sig_tr_dev = norm(sig_tr_dev);
        n = sig_tr_dev/norm_sig_tr_dev;
        
        % calculation of dlambda with newtons scheme (local iteration)
%         f = @(x) dt*fth*(sinh(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*x)/(kappaGn+Yth)-1,0)))^nth-x;
%         df = @(x) dt*fth*nth*(sinh(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*x)/(kappaGn+Yth)-1,0)))^(nth-1)*cosh(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*x)/(kappaGn+Yth)-1,0))*(-3*mu/(kappaGn+Yth))-1;

        % calculation of dlambda with newtons scheme (local iteration) WITHOUT SINH()!
        f = @(x) dt*fth*(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*x)/(kappaGn+Yth)-1,0) )^nth-x;
        df = @(x) dt*fth*nth*(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*x)/(kappaGn+Yth)-1,0))^(nth-1)*(-3*mu/(kappaGn+Yth))-1;
        
        x0 = 2;
        [temp,err] = newton(f, df, x0, tol, nmax);
        z = size(err',1);
        error_local(1:z) = err';
        dlambda = temp(end);
        
        % updating plastic multiplier and plastic strain
        epsp = epspn + sqrt(3/2)*dlambda*n;
        epspv = epspvn + dlambda;
        
        % calculation of stresses
        sig = sig_tr_dev - 2*mu*dlambda*sqrt(3/2)*n + K0*(id2*id2')'*eps;
        sigv = sqrt(0.5*((sig(1)-sig(2))^2 + (sig(2)-sig(3))^2 + (sig(3)-sig(1))^2+6*(sig(4)^2+sig(5)^2+sig(6)^2)));
 
        % tangent modulus
%         dres_deps = dt*fth*nth*(sinh(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*dlambda)/(kappaGn+Yth)-1,0)))^(nth-1)*cosh(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*dlambda)/(kappaGn+Yth)-1,0))*sqrt(3/2)*2*mu*n/(kappaGn+Yth);
        dres_deps = dt*fth*nth*(max(sqrt(3/2)*(norm_sig_tr_dev-2*mu*sqrt(3/2)*dlambda)/(kappaGn+Yth)-1,0))^(nth-1)*sqrt(3/2)*2*mu*n/(kappaGn+Yth);
        C_T = K0*(id2*id2') + 2*mu*id4_dev - sqrt(3/2)*dlambda*(2*mu)^2/norm_sig_tr_dev*(id4_dev-n*n')+2*mu/df(dlambda)*sqrt(3/2)*n*dres_deps';
        
        %% calculation of internal hardening variables for each recrystallization cycle
        for irx = 1:nrx+1
            dX(irx+1,1) = 1/mu*exp(-ctheta/T)*(1-exp(-B*(zetan(irx,1)/mu)^m))*(ckappa*kappan(irx,1)^2+czeta*zetan(irx,1)^2)*Xn(irx,1)*(Xn(irx+1,1)/Xn(irx,1))^a*(1-Xn(irx+1,1)/Xn(irx,1))^b*(1+c*(1-Xn(irx,1)));
            X(irx+1,1)  = Xn(irx+1,1) + dt*dX(irx+1,1);
            
            dzeta(irx,1) = zetan(irx,1)/mu*dmu*Tdot + hc*mu*(zetan(irx,1)/mu)^(1-1/r)*abs(dlambda/dt) - zetan(irx,1)*dX(irx,1)/(Xn(irx,1)-Xn(irx+1,1));
            zeta(irx,1) = zetan(irx,1) + dt*dzeta(irx,1);
            
            dkappa(irx,1) = kappan(irx,1)/mu*dmu*Tdot + (c8*mu*(1+zetan(irx,1)/kappan(irx,1)) - c6*exp(-c7/T)*kappan(irx,1))*abs(dlambda/dt) - kappan(irx,1)*dX(irx,1)/(Xn(irx,1)-Xn(irx+1,1));
            kappa(irx,1) = kappan(irx,1) + dt*dkappa(irx,1);
        end
        
        % rule of mixture for averaging the final hardening variable
        for irx = 1:nrx
            kappaX(irx,1) = kappa(irx,1)*(X(irx,1)-X(irx+1,1));
        end
        kappaG = sum(kappaX(1:nrx,1)) + kappa(nrx+1,1)*X(nrx+1,1);
        
        %% history field     
        his2 = vertcat(epsp,X,kappa,zeta,kappaG,epspv);
        
        % insert  plot variables
        [plotv] = mat_recrystallization_ip(sig,sigv,epspv,kappa,zeta,X,kappaX,kappaG);
        OutputIntegerVariables = horzcat(npar, nhis, npvpf5);
end

    function [tpvpf,npvpf5] = mat_recrystallization_pp(iptyp)
        if iptyp == 3
            tpvpf(1) = {'SRR          '};
            tpvpf(2) = {'STT          '};
            tpvpf(3) = {'SZZ          '};
            tpvpf(4) = {'SRZ          '};
        elseif iptyp == 4
            tpvpf(1) = {'SXX              '};
            tpvpf(2) = {'SYY              '};
            tpvpf(3) = {'SZZ              '};
            tpvpf(4) = {'SXY              '};
            tpvpf(5) = {'SYZ              '};
            tpvpf(6) = {'SXZ              '};
        else
            tpvpf(1) = {'SXX          '};
            tpvpf(2) = {'SYY          '};
            tpvpf(3) = {'SZZ          '};
            tpvpf(4) = {'SXY          '};
        end
        
        if iptyp < 4
            tpvpf(7)  = {'SIGV             '};
        else
            tpvpf(7)  = {'SIGV             '};
            tpvpf(8)  = {'EPSV             '};
            tpvpf(9)  = {'K1               '};
            tpvpf(10) = {'K2               '};
            tpvpf(11) = {'K3               '};
            tpvpf(12) = {'K4               '};
            tpvpf(13) = {'K5               '};
            tpvpf(14) = {'K6               '};
            tpvpf(15) = {'K7               '};
            tpvpf(16) = {'K8               '};
            tpvpf(17) = {'K9               '};
            tpvpf(18) = {'Z1               '};
            tpvpf(19) = {'Z2               '};
            tpvpf(20) = {'Z3               '};
            tpvpf(21) = {'Z4               '};
            tpvpf(22) = {'Z5               '};
            tpvpf(23) = {'Z6               '};
            tpvpf(24) = {'Z7               '};
            tpvpf(25) = {'Z8               '};
            tpvpf(26) = {'Z9               '};
            tpvpf(27) = {'X1               '};
            tpvpf(28) = {'X2               '};
            tpvpf(29) = {'X3               '};
            tpvpf(30) = {'X4               '};
            tpvpf(31) = {'X5               '};
            tpvpf(32) = {'X6               '};
            tpvpf(33) = {'X7               '};
            tpvpf(34) = {'X8               '};
            tpvpf(35) = {'X9               '};
            tpvpf(36) = {'X10              '};
            tpvpf(37) = {'KX1              '};
            tpvpf(38) = {'KX2              '};
            tpvpf(39) = {'KX3              '};
            tpvpf(40) = {'KX4              '};
            tpvpf(41) = {'KX5              '};
            tpvpf(42) = {'KX6              '};
            tpvpf(43) = {'KX7              '};
            tpvpf(44) = {'KX8              '};
            tpvpf(45) = {'KX9              '};
            tpvpf(46) = {'KG               '};
        end
        
        npvpf5 = size(tpvpf,2);
        
        if npvpf5 ~= 46
            error('uma_mises: nur iptyp = 4 moeglich')
        end
    end

    function [plotv] = mat_recrystallization_ip (sig,sigv,epspv,kappa,zeta,X,kappaX,kappaG)
        plotv = vertcat(sig,sigv,epspv,kappa,zeta,X,kappaX,kappaG);
    end

end