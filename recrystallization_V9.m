% model for dynamic recrystallization V9
% mode - dependent parameters 
% this subroutine calculates nonlinear hardening during dynamic
% recrystallization based on the continuum model proposed by brown.

% it accounts for single peak and multiple peak dynamic recrystallization
% for (nrx) continuous cycles of recrystallization

% input:                strain (eps), temperature (T), material parameters
% output:               stress (sig)
% internal variables:   hardening variables/cycle of rx (kappa), 
%                       misorientation variables/cycle of rx (zeta), 
%                       averaged hardening variable (kappaG)

% - radial return solves for norm_sig_dev (newton)
% - the nonlinear evolutions of the internal variables are updated with an
%   explicit algorithm
%
%  therefore this routine can be considered as semi implicit

clear; clc; close all; tic;

% temp = constant

%% incrementation
tfac = 2;
tmax = 1500*tfac;              % total time
nt =  8000*tfac;               % number of increments

t = linspace(0,tmax,nt);
dt = tmax/nt;
nmax = 200;               % maximum number of iterations
tol = 1e-8;              % tolerance for newtons scheme
nrx = 8;                  % number of concurrent recrystallization cycles        % variable input parameter

nstr = 6;                 % number of stress components

mode = [1,1,0,0,0,0];     % on/off switch for the different models
nmod = nnz(mode);         % number of different recrystallization models

ntd = nt; %   nt;         % dimension of time dependent arrays


% Temperatur
temp = 337 + 273.15;  


%% Material parameters for stress - based formulation, mode 2 (original Bammann)

shear  = 5.47e4;  % RM 
dshear = - 34.1;  % RM

B = 1.5e9;
m = 2.88;
a = 0.225;
b = 1.31;
c = 25;
%czeta  = 8.7e11;
ckappa = 8.7e11;
ctheta = 1.8e4;
c1 = 0.214;
c2 = 5.22;
c3 = 16.4;
c4 = 0;
c5 = 4.64e-4;
c6 = 8.8;
c7 = 241;
c8 = 0.0137;

hc = 2.17e-5;
r  = 0.657;
be = 0.25;        
ceps   = 5e-10; % (2*c8)/(be*k1);
calpha = 5e-10;    
clg = 1e-8;


% Define 
k1(1) = 1; 
% not needed
%k2(1) = 1; 
%k3(1) = 1;

czeta  = ckappa/k1(1); 

%% Material parameters for strain - based formulation, mode 2

k1(2)     = 1; 
k2(2)     = 1;   % ohne Einfluss ?
k3(2)     = 2*4;   

c6(2)     = 1.5*c6(1);
c8(2)     = 2*c8(1);

ckappa(2) = ckappa(1);          
czeta(2)  = ckappa(2)/k1(2);

calpha(2) = calpha(1);   

ceps(2)   = c8(2)/be; 
calpha(2) = k1(2)*ceps(2); 
clg(2)    = hc*(calpha(2)*be*k2(2)/k1(2))^(-1/r);

% Kontrolle
hch  = clg(2)*(calpha(2)*be*k2(2)/k1(2))^(1/r);  

%% Material parameters for own stress - based formulation, mode 3 
 
factor2 = 2;
c6(3)     = c6(1)*factor2; 
c8(3)     = c8(1)*factor2;
ckappa(3) = ckappa(1);

k1(3)     = 1;
czeta(3)  = ckappa(3)/k1(3);
  
% k2(3)     = k2(1); 
k3(3)     = k3(1);   

%% initializing variables
T         = zeros(ntd,1);
eps       = zeros(ntd,nstr);
eps_dev   = zeros(ntd,nstr);
epsp      = zeros(ntd,nstr);
epsel     = zeros(ntd,nstr);

sig             = zeros(ntd,nstr);
sig_tr_dev      = zeros(ntd,nstr);
norm_sig_tr_dev = zeros(ntd,1);
norm_sig_dev    = zeros(ntd,1);
mises           = zeros(ntd,1);
n               = zeros(ntd,nstr);
dlambda         = zeros(ntd,1);
epsv            = zeros(ntd,1);

error_local     = zeros(nmax,ntd);
error           = zeros(nmax,1);

mu  = zeros(ntd,1);
dmu = zeros(ntd,1);

fth = zeros(ntd,1);
nth = zeros(ntd,1);
Yth = zeros(ntd,1);

X      = zeros(ntd,nrx+2,nmod);
dX     = zeros(ntd,nrx+2,nmod);
dkappa = zeros(ntd,nrx+1,nmod);
kappa  = zeros(1,nrx+1,nmod);
dzeta  = zeros(ntd,nrx+1,nmod);
zeta   = zeros(ntd,nrx+1,nmod);

kappaX = zeros(ntd,nrx+1,nmod);
kappaG = zeros(ntd,1,nmod);

rho    = zeros(ntd,nrx+1,nmod);
frho   = zeros(ntd,nrx+1,nmod);
invLg  = zeros(ntd,nrx+1,nmod);
epss   = zeros(ntd,nrx+1,nmod);

alphag = zeros(ntd,nrx+1,nmod);

V    = zeros(ntd,nrx+1,nmod);
dV   = zeros(ntd,nrx+1,nmod);
VG   = zeros(ntd,1,nmod);
QG   = zeros(ntd,1,nmod);

Rrx  = zeros(ntd,nrx+2,nmod);
Qrx1 = zeros(ntd,nrx+2,nmod);
Qrx2 = zeros(ntd,nrx+2,nmod);
grx  = zeros(ntd,nrx+2,nmod);

%% strain and temperature
Tdot = 0; % -0.002;
epsdot = [0.0004 0 0 0 0 0];    % variable input parameter
id = [1 1 1 0 0 0];

epsdot(2) = -epsdot(1)/2;
epsdot(3) =  epsdot(2);

%% initial values
T(1) = temp;       % temperature 

kappa(1,:,:) = 1e-6;
zeta(1,:,:)  = 1e-6;
X(:,1,:)     = 1;

for irx = 1:nrx+1
    X(1,irx+1,:) = 1e-6*X(1,irx,:);  % ??????
end

p = zeros(nrx+1,1);
index = 0;

invLg(1,:,:) = 1e-6;
dV(:,:,:)    = 1e-6;
rho(1,:,:)   = 1e-6;
 
%% time discretization
for i = 1:nt-1 
    %% calculation of strain and temperature
    T(i+1) = T(i) + dt*Tdot;
    eps(i+1,:) = eps(i,:) + dt*epsdot;
    
    mu(i+1)  =  shear + dshear*T(i+1);
    dmu(i+1) =  dshear; % RM
    
    fth(i+1) = c1*exp(-c2/T(i+1));
    nth(i+1) = c3+c4/T(i+1);
    Yth(i+1) = c5*mu(i+1);

    %% calculation of internal hardening variables and stresses
    
    eps_dev(i+1,:) = eps(i+1,:)-(eps(i+1,1)+eps(i+1,2)+eps(i+1,3))/3*id;
  
    sig_tr_dev(i+1,:) = 2*mu(i+1)*(eps_dev(i+1,:)-epsp(i,:));
    norm_sig_tr_dev(i+1) = norm(sig_tr_dev(i+1,:));
    n(i+1,:) = sig_tr_dev(i+1,:)/norm_sig_tr_dev(i+1);
    
    % calculation of norm_sig_dev with newtons scheme (local iteration)
       
    kmod = 1;   % Mode fuer Fliessfunktion 
%     A = @(x) dt*fth(i+1)*(sinh(max((sqrt(3/2)*x)/(kappaG(i,kmod)+Yth(i+1))-1,0)))^nth(i+1);
%     f = @(x) x-norm_sig_tr_dev(i+1)+2*mu(i+1)*sqrt(3/2)*A(x);
%     dA = @(x) dt*fth(i+1)*nth(i+1)*(sinh(max((sqrt(3/2)*x)/(kappaG(i,kmod)+ Yth(i+1))-1,0)))^(nth(i+1)-1)*cosh(max((sqrt(3/2)*x)/(kappaG(i,kmod)+Yth(i+1))-1,0))*sqrt(3/2)/(kappaG(i,kmod) + Yth(i+1));
%     df = @(x) 1+2*mu(i+1)*sqrt(3/2)*dA(x);
    
    A = @(x) dt*fth(i+1)*(max((sqrt(3/2)*x)/(kappaG(i,kmod)+Yth(i+1))-1,0))^nth(i+1);
    f = @(x) x-norm_sig_tr_dev(i+1)+2*mu(i+1)*sqrt(3/2)*A(x);
    dA = @(x) dt*fth(i+1)*nth(i+1)*(max((sqrt(3/2)*x)/(kappaG(i,kmod)+ Yth(i+1))-1,0))^(nth(i+1)-1)*sqrt(3/2)/(kappaG(i,kmod) + Yth(i+1));
    df = @(x) 1+2*mu(i+1)*sqrt(3/2)*dA(x);
        
%     [x0,x_store,r_store] = initial(f, nmax);
    x0 = 2;
    [temp,error] = newton(f, df, x0, tol, nmax);
    z = size(error',1);
    error_local(1:z,i) = error';
    norm_sig_dev(i+1) = temp(end);
    
    % updating plastic multiplier and plastic strain
    
    dlambda(i+1) = (norm_sig_tr_dev(i+1)-norm_sig_dev(i+1))/(2*mu(i+1)*sqrt(3/2));     
    if abs(dlambda(i))  ~=  real(dlambda(i))   % RM
       stop dlambda is complex
    end
    
    epsp(i+1,:)  = epsp(i,:) + sqrt(3/2)*dlambda(i+1)*n(i+1,:);
    epsel(i+1,:) = eps(i+1,:) - epsp(i+1,:);
    epsv(i+1)    = epsv(i) + dlambda(i+1);
    

    % calculation of stresses  
    sig(i+1,:) = sig_tr_dev(i+1,:) - 2*mu(i+1)*dlambda(i+1)*sqrt(3/2)*n(i+1,:) + mu(i+1)*id*(epsel(i+1,1)+epsel(i+1,2)+epsel(i+1,3));
    mises(i+1) = sqrt(0.5*((sig(i+1,1) - sig(i+1,2))^2 + (sig(i+1,2)-sig(i+1,3))^2 + (sig(i+1,3)-sig(i+1,1))^2+6*(sig(i+1,4)^2+sig(i+1,5)^2+sig(i+1,6)^2)));
    
    %% calculation of internal hardening variables for each recrystallization cycle

    rho(i+1,:,:)   = 0;
    invLg(i+1,:,:) = 0;
    X(i+1,1,:)     = 1;
    QG(i+1,:)      = 0; 
    
    j = 0;

     if mode(1) == true
        % Bammann model
        j = j+1;
        jp = 1; 
        
        VG(i+1,j) = 0;
        
        for irx = 1:nrx+1
            
            Rrx(i+1,irx+1,j)  = 1/mu(i+1)*exp(-ctheta/T(i+1));
            Qrx1(i+1,irx+1,j) = (1-exp(-B*(zeta(i,irx,j)/mu(i+1))^m));
            Qrx2(i+1,irx+1,j) = (ckappa(jp)*kappa(i,irx,j)^2 + czeta(jp)*zeta(i,irx,j)^2);
            grx(i+1,irx+1,j)  = X(i+1,irx,j)*(X(i,irx+1,j)/X(i+1,irx,j))^(a)*(1-X(i,irx+1,j)/X(i+1,irx,j))^b * (1+c*(1-X(i+1,irx,j)));

            dX(i+1,irx+1,j) = Rrx(i+1,irx+1,j)*Qrx1(i+1,irx+1,j)*Qrx2(i+1,irx+1,j)*grx(i+1,irx+1,j);
            X(i+1,irx+1,j)  = X(i,irx+1,j) + dt*dX(i+1,irx+1,j);
            V(i+1,irx,j)    = X(i+1,irx,j)-X(i+1,irx+1,j);
%           dV(i+1,irx,j)   = dX(i+1,irx,j)-dX(i+1,irx+1,j)/dt;
            VG(i+1,j)       = VG(i+1,j) + V(i+1,irx,j);
              
            if abs(X(i+1,irx+1,j))  ~=  real(X(i+1,irx+1,j))   % RM
                stop mode = 1  X is complex
            end
  
            dXdV = dt*dX(i+1,irx,j)/V(i+1,irx,j); 
            
            Htemp(jp) = c8(jp)*mu(i+1);
            Rd =  c6(jp)*exp(-c7/T(i+1));
                                      
            dzeta(i+1,irx,j) = zeta(i,irx,j)/mu(i+1)*dmu(i+1)*Tdot + hc*mu(i+1)*(zeta(i,irx,j)/mu(i+1))^(1-1/r)*abs(dlambda(i+1)/dt) - zeta(i,irx,j)*dX(i+1,irx,j)/V(i+1,irx,j);
            zeta(i+1,irx,j) = zeta(i,irx,j) + dt*dzeta(i+1,irx,j);
    
            dkappa(i+1,irx,j) = kappa(i,irx,j)/mu(i+1)*dmu(i+1)*Tdot + (Htemp(jp)*(1+zeta(i,irx,j)/kappa(i,irx,j)) - Rd*kappa(i,irx,j))*abs(dlambda(i+1)/dt) - kappa(i,irx,j)*dX(i+1,irx,j)/V(i+1,irx,j);
            kappa(i+1,irx,j) = kappa(i,irx,j) + dt*dkappa(i+1,irx,j);
                  
%            VG(i+1,j) = sum(V(i+1,:,j));
                    
%  postprocessing: compute strain type variables (see Appendix, A2)
            invLg(i+1,irx,j)  = k1(jp)*zeta(i+1,irx,j)/(be*k2(jp)*mu(i+1)*calpha(jp));         
            alphagh           = zeta(i+1,irx,j)/(mu(i+1)*calpha(jp));  
            alphag(i+1,irx,j) = (be*k2(jp))/k1(jp)*invLg(i+1,irx,j);             
             
            rho(i+1,irx,j)    = (kappa(i+1,irx,j)/(mu(i+1)*ceps(jp)*be))^2;         
            epssh             =  kappa(i+1,irx,j)/(mu(i+1)*ceps(jp)) ; 
            epss(i+1,irx,j)   =  be*sqrt(rho(i+1,irx,j));  
            
%          frho(i+1,irx,j) = (k1*sqrt(rho(i,irx,j)) + k2*invLg(i+1,irx,j) - c6*exp(-c7/T(i+1))*rho(i,irx,j));
            
        end
%         czeta(3)  = czeta(1);

        % rule of mixture for averaging the final hardening variable
        
        for irx = 1:nrx
            kappaX(i+1,irx,j) = kappa(i+1,irx,j)*(X(i+1,irx,j)-X(i+1,irx+1,j));
        end
        
        % very important point (not accurate in bammanns paper or just simply wrong)
        
        kappaG(i+1,j) = sum(kappaX(i+1,1:nrx,j)) + kappa(i+1,nrx+1,j)*(X(i+1,nrx+1,j)-X(i+1,nrx+2,j));
  
        QG(i+1,j) = kappaG(i+1,j); 
     
     end
     
    if mode(2) == true
        % own strain based version
        jp = 2; 
        j = j+1;
        VG(i+1,j) = 0; 
        
        for irx = 1: nrx+1
            
            Rrx(i+1,irx+1,j)  = 1/mu(i+1)*exp(-ctheta/T(i+1));
            Qrx1(i+1,irx+1,j) = (1-exp(-B*(zeta(i,irx,j)/mu(i+1))^m));
            Qrx2(i+1,irx+1,j) = (ckappa(jp)*kappa(i,irx,j)^2 + czeta(jp)*zeta(i,irx,j)^2);
            grx(i+1,irx+1,j)  = X(i+1,irx,j)*(X(i,irx+1,j)/X(i+1,irx,j))^(a)*(1-X(i,irx+1,j)/X(i+1,irx,j))^b * (1+c*(1-X(i+1,irx,j)));
            
            dX(i+1,irx+1,j) = Rrx(i+1,irx+1,j)*Qrx1(i+1,irx+1,j)*Qrx2(i+1,irx+1,j)*grx(i+1,irx+1,j);
            X(i+1,irx+1,j)  = X(i,irx+1,j) + dt*dX(i+1,irx+1,j);
            V(i+1,irx,j)    = X(i+1,irx,j)-X(i+1,irx+1,j);
%           dV(i+1,irx,j)   = dX(i+1,irx,j)-dX(i+1,irx+1,j)/dt;
            VG(i+1,j)       = VG(i+1,j) + V(i+1,irx,j);
              
            if abs(X(i+1,irx+1,j))  ~=  real(X(i+1,irx+1,j))   % RM
                stop mode = 2  X is complex
            end
                        
            Rd      = c6(jp)*exp(-c7/T(i+1)); 
            
            fgdl    = clg(jp)*(invLg(i,irx,j))^(1-1/r)*dlambda(i+1); 
%           falpdl  = (be*k2(jp))/k1(jp)*fgdl; 
            frhodl  = (k1(jp)*sqrt(rho(i,irx,j)) + k2(jp)*invLg(i+1,irx,j)- Rd *rho(i,irx,j))*dlambda(i+1);  

            frho(i+1,irx,j) =  frhodl/dlambda(i+1);
%           fepsdl  = b*frhodl/(2*sqrt(rho(i,irx,j)));
        
            dXdV = dt*dX(i+1,irx,j)/V(i+1,irx,j);     
            
            invLg(i+1,irx,j) = invLg(i,irx,j)*(1-dXdV) + fgdl;             
            rho(i+1,irx,j)   = rho  (i,irx,j)*(1-dXdV) + frhodl;
                      
%  postprocessing: compute stress type variables
            
            alphag(i+1,irx,j) = (be*k2(jp))/k1(jp)*invLg(i+1,irx,j);                   
            epss(i+1,irx,j)   = sqrt(rho(i+1,irx,j))*be;       
            
            zeta(i+1,irx,j)  = calpha(jp)*mu(i+1)*alphag(i+1,irx,j);             
            kappa(i+1,irx,j) = ceps(jp)*mu(i+1)*epss(i+1,irx,j);  
            
%            QG(i+1,j) = QG(i+1,j) + 15*mu(i+1)* ceps*frho(i+1,irx,j)* epss(i+1,irx,j)*epss(i+1,irx,j)*V(i+1,irx,j)/...
%                (2*rho(i+1,irx,j)) ;     
%            QG(i+1,j) = QG(i+1,j) + 15* kappa(i+1,irx,j) *frho(i+1,irx,j)*epss(i+1,irx,j)*V(i+1,irx,j)/...
%                (2*rho(i+1,irx,j)) ;     
%           QG(i+1,j) = QG(i+1,j) + mu(i+1)* ceps*frho(i+1,irx,j)*be*be*V(i+1,irx,j)/2;   
%            QG(i+1,j) = QG(i+1,j) + 3*mu(i+1)* ceps*frho(i+1,irx,j)*be*be/2;   
%           QG(i+1,j) = QG(i+1,j) + be*(k1*kappa(i+1,irx,j) +  zeta(i+1,irx,j))/2;   

            fk  = be*k1(jp)*k3(jp)/2; 
            fz1 = be*k1(jp)*ceps(jp)/(calpha(jp)*2);
            fz2 = hc*(zeta(i+1,irx,j)/mu(i+1))^(1-1/r)/calpha(jp);
            fz  = fz1 + fz2; 
            
%             fk = 1.0; fz = 0;  % simple yield stress 
            QG(i+1,j) = QG(i+1,j) + (fk*kappa(i+1,irx,j) +  fz*zeta(i+1,irx,j))*V(i+1,irx,j);  
                                   
        end 
        
        for irx = 1:nrx
            kappaX(i+1,irx,j) = kappa(i+1,irx,j)*(V(i+1,irx,j));
        end
        
        kappaG(i+1,j) = sum(kappaX(i+1,1:nrx,j));
%        QG(i+1,j) = sum(kappa(i+1,irx)*epss(i+1,irx)/rho(i+1,irx)*frho(i+1,irx))/2;
%        QG(i+1,j) = mu(i+1)/2*ceps*be^2*sum(frho(i+1,irx));
        
         kappaG(i+1,j) = QG(i+1,j);

 % Gleichung von mode 1         
%          kappaG(i+1,j) = sum(kappaX(i+1,1:nrx,j)) + kappa(i+1,nrx+1,j)*(X(i+1,nrx+1,j)-X(i+1,nrx+2,j));  
%           QG(i+1,j) = kappaG(i+1,j); 
         
    end

    
     if mode(3) == true
        % Own stress based formulation: difference is factor2
        jp = 3; 
        j = j+1;
        
        VG(i+1,j) = 0;
        
        for irx = 1:nrx+1
            
            Rrx(i+1,irx+1,j)  = 1/mu(i+1)*exp(-ctheta/T(i+1));
            Qrx1(i+1,irx+1,j) = (1-exp(-B*(zeta(i,irx,j)/mu(i+1))^m));
            Qrx2(i+1,irx+1,j) = (ckappa(jp)*kappa(i,irx,j)^2 + czeta(jp)*zeta(i,irx,j)^2);
            grx(i+1,irx+1,j)  = X(i+1,irx,j)*(X(i,irx+1,j)/X(i+1,irx,j))^a*...
                (1-X(i,irx+1,j)/X(i+1,irx,j))^b * (1+c*(1-X(i+1,irx,j)));
            
            dX(i+1,irx+1,j) = Rrx(i+1,irx+1,j)*Qrx1(i+1,irx+1,j)*Qrx2(i+1,irx+1,j)*grx(i+1,irx+1,j);
            X(i+1,irx+1,j)  = X(i,irx+1,j) + dt*dX(i+1,irx+1,j);
            V(i+1,irx,j)    = X(i+1,irx,j)-X(i+1,irx+1,j);
%           dV(i+1,irx,j)   = dX(i+1,irx,j)-dX(i+1,irx+1,j)/dt;
            VG(i+1,j)       = VG(i+1,j) + V(i+1,irx,j);
              
            if abs(X(i+1,irx+1,j))  ~=  real(X(i+1,irx+1,j))   % RM
                stop mode = 3  X is complex
            end
  
            dXdV = dt*dX(i+1,irx,j)/V(i+1,irx,j); 
            
            fzet     = hc*mu(i+1)*(zeta(i,irx,j)/mu(i+1))^(1-1/r);
            fzetdl   = fzet * dlambda(i+1);        
           
            Htemp = c8(jp)*mu(i+1); Rd =  c6(jp)*exp(-c7/T(i+1)); 

            fkap     = (Htemp*(k1(jp)*k3(jp) + ...
                  k1(jp)*ceps(jp)*zeta(i,irx,j)/(calpha(jp)*kappa(i,irx,j))) - ...
                  kappa(i,irx,j)*Rd)/factor2;
            fkapdl   = fkap * dlambda(i+1); 
            
            zeta(i+1,irx,j)  = zeta(i,irx,j) *(1 + dt*dmu(i+1)*Tdot/mu(i+1)-dXdV) + fzetdl;                                             
            kappa(i+1,irx,j) = kappa(i,irx,j)*(1 + dt*dmu(i+1)*Tdot/mu(i+1)-dXdV) + fkapdl;                                  
                                 
%            VG(i+1,j) = sum(V(i+1,:,j));
                       
%  postprocessing: compute strain type variables (see Appendix, A2)                            

            invLg(i+1,irx,j)  = k1(jp)*zeta(i+1,irx,j)/(be*k2(jp)*mu(i+1)*calpha(jp));         
            alphagh           = zeta(i+1,irx,j)/(mu(i+1)*calpha(jp));  
            alphag(i+1,irx,j) = (be*k2(jp))/k1(jp)*invLg(i+1,irx,j);             
             
            rho(i+1,irx,j)    = (kappa(i+1,irx,j)/(mu(i+1)*ceps(jp)*be))^2;         
            epssh             =  kappa(i+1,irx,j)/(mu(i+1)*ceps(jp)) ; 
            epss(i+1,irx,j)   =  be*sqrt(rho(i+1,irx,j));  
            
            xx = 1; 
            
%           epss(i+1,irx,j)   = sqrt(rho(i+1,irx,j))*be;      
        end
        
        % rule of mixture for averaging the final hardening variable
        
        for irx = 1:nrx
            kappaX(i+1,irx,j) = kappa(i+1,irx,j)*(X(i+1,irx,j)-X(i+1,irx+1,j));
        end
        
       % kappaX(i+1,irx,j) 
        
     % very important point (not accurate in bammanns paper or just simply wrong)
        
        kappaG(i+1,j) = sum(kappaX(i+1,1:nrx,j)) + kappa(i+1,nrx+1,j)*(X(i+1,nrx+1,j)-X(i+1,nrx+2,j));
  
        QG(i+1,j) = kappaG(i+1,j); 
     end
     
  

    
    
    
end

%% command line output
disp('------------------------')
disp('Completed with no errors.')
disp(['CPU Time    ',num2str(cputime),'s'])
disp(['Wall Time   ',num2str(toc),'s'])
disp('------------------------')
disp(['Total Time  ',num2str(tmax),'s'])
disp(['Timesteps   ',num2str(nt)])
disp(['Cycles      ',num2str(nrx)])
disp(['Modes       ',num2str(find(mode))])
disp('------------------------')

%% graphical output

% input/output
figure('Name','Input and Stresses','NumberTitle','off','units','normalized','outerposition',[0 0 1 1])
%subplot(2,2,1);
subplot(4,4,1);
%yyaxis left
plot(t,eps(:,1),'black')
hold on
title('Strain and Temperature')
ylabel(['Strain ',' [-] '])
%yyaxis right
plot(t,T-273.15,'r')
axis([0 tmax 0 T(1)*1.2]);
ylabel('Temperature T [°C] ')
xlabel('Time [s]')
% legend(char(949),'T','Location','NorthWest')
subplot(4,4,2); plot_function(eps(:,1),mises,'True Strain',char(949),'[-]','Equivalent Stress','\sigma_v','[MPa]','\sigma_v',mode)
subplot(4,4,3); plot_function(t,epsv,'Time','t','[s]','Equivalent Plastic Strain','epsv','[-]','EV',mode)

%internal hardening variables
%figure('Name','Internal Variables','NumberTitle','off','units','normalized','outerposition',[0 0 1 1])
subplot(4,4,4); 
for k = 1:j
  plot(eps(:,1),kappaG(:,1,k),'black')
  hold on
end
plot(eps(:,1),kappaX(:,:))
title('Averaged Hardening Variable \kappa_G '); ylabel('Averaged Hardening Variable \kappa_G [MPa] '); xlabel('True Strain [-]')
% legend({'KG1','KG2','KV1','KV2','KV3','KV4','KV5','KV6'},'Location','SouthWest')

% internal hardening variables
subplot(4,4,5); plot_function(eps(:,1),kappa,'True Strain',char(949),'[-]','Hardening Variable ','\kappa','[MPa]','K',mode)
subplot(4,4,6); plot_function(eps(:,1),X,'True Strain',char(949),'[-]','Recrystallized Fraction','X','[-]','X',mode)
subplot(4,4,7); plot_function(eps(:,1),zeta,'True Strain',char(949),'[-]','Misorientation Variable','\zeta','[MPa]','Z',mode)

% dislocation density and average spacing
%figure('Name','Dislocation Density','NumberTitle','off','units','normalized','outerposition',[0 0 1 1])
subplot(4,4,8); plot_function(eps(:,1),rho,'True Strain',char(949),'[-]','Dislocatin Density','\rho_s','[m/m^3]','R',mode)
subplot(4,4,9); plot_function(eps(:,1),invLg,'True Strain',char(949),'[-]','Average Spacing','1/Lg','[m/m^3]','1/Lg',mode)
subplot(4,4,10); plot_function(eps(:,1),V,'True Strain',char(949),'[-]','Volume fraction','V','[-]','V',mode)
subplot(4,4,11); plot_function(eps(:,1),QG,'True Strain',char(949)','[-]','Hardening Variable','Q','[MPa]','Q',mode)

% evolution of x
%figure('Name','Recrystallized Fraction X - Components','NumberTitle','off','units','normalized','outerposition',[0 0 1 1])
subplot(4,4,12); plot_function(eps(:,1),Qrx1,'True Strain',char(949)','[-]','Part 1 -',['Q(',char(949),',\alpha)'],'[?]','Q',mode)
subplot(4,4,13); plot_function(eps(:,1),Qrx2,'True Strain',char(949)','[-]','Part 2 -',['Q(',char(949),',\alpha)'],'[?]','Q',mode)
subplot(4,4,14); plot_function(eps(:,1),grx,'True Strain',char(949)','[-]','Correction Function','g(X)','[?]','g',mode)
subplot(4,4,15); plot_function(eps(:,1),Rrx,'True Strain',char(949)','[-]','Temperature Function','R(\theta)','[?]','R',mode)

figure(3)
plot_function(eps(:,1),QG,'True Strain',char(949)','[-]','Hardening Variable','Q','[MPa]','Q',mode)