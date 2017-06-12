%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann Method        %
% of a channel flow using a D2Q9  %
% speed model and single-         %
% relaxation time collision model %
%                                 %
% Tom-Robin Teschner 2017         %
% Cranfield University            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear screen and cache
clc;
clear;
close all;

%% initialise all 

% LB fluid information
% variables===========================description==========================
nu     = 0.1;                       % lattice viscosity
uref   = 0.1;                       % reference velocity
omega  = 1/(3*nu+0.5);              % collision frequency

% Geometrical information
% variables===========================description==========================
nx     = 100;                       % number of nodes in x
ny     = 100;                       % number of nodes in y
lx     = 1;                         % domain length in x
ly     = 1;                         % domain length in y
dx     = lx/(nx-1);                 % mesh spacing in x
dy     = ly/(ny-1);                 % mesh spacing in y

% calculation specifics
% variables===========================description==========================
tmax   = 10000;                     % total amount of timesteps
t      = 1;                         % first timestep, default=1
eps    = 10e-5;                     % convergence criterion (L2 norm)
ures   = 1;                         % initial residual of u velocity
vres   = 1;                         % initial residual of v velocity
rres   = 1;                         % initial residual of density
rf     = 1;                         % report frequency of iterations to screen
pa     = 0;                         % plot animation

% initialise arrays
% variables===========================description==========================
f      = zeros(nx,ny,9);            % distribution function at time level n
fn05   = zeros(nx,ny,9);            % intermediate distribution function
fn1    = zeros(nx,ny,9);            % distribution function at time level n+1
feq    = zeros(nx,ny,9);            % distribution function at equilibrium
rho    = zeros(nx,ny);              % lattice density
u      = zeros(nx,ny);              % lattice velocity in x
v      = zeros(nx,ny);              % lattice velocity in v
x      = zeros(nx,ny);              % x spacing vector
y      = zeros(nx,ny);              % y spacing vector

% ===== set speedmodel's directions =======================================
cx = [0 1 0 -1 0 1 -1 -1 1];
cy = [0 0 1 0 -1 1 1 -1 -1];

% ===== set speedmodel's weights ==========================================
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];

% ===== lattice (mesh) ====================================================
for i=1:nx
    for j=1:ny
        x(i,j) = dx*(i-1);
        y(i,j) = dy*(j-1);
    end
end

% ===== initialise arrays =================================================
for i=1:nx
    for j=1:ny
        for k=1:9
            % initialise distribution function with small, non-zero values
            f(i,j,k)   = 1.0*w(k);
            fn1(i,j,k) = 1.0*w(k);
            feq(i,j,k) = 1.0*w(k);
        end
        % initialise velocity and density
        rho(i,j) = 1.0;
        u(i,j)   = 0.0;
        v(i,j)   = 0.0;
    end
end

%% Run Simulation

% ===== main timeloop =====================================================
while(t<=tmax)

% ===== write updated values into old variable array ======================
    f=fn1;
    
% ===== variables needed for calculation of residuals =====================
    usum0 = 0;
    vsum0 = 0;
    rsum0 = 0;

    for i=1:nx
        for j=1:ny
            usum0 = usum0 + u(i,j);
            vsum0 = vsum0 + v(i,j);
            rsum0 = rsum0 + rho(i,j);
        end
    end
   
% ===== perform collision step ============================================
    for i=1:nx
        for j=1:ny
            % temporary variable 1
            cu1 = u(i,j)*u(i,j) + v(i,j)*v(i,j);
            
            for k=1:9
                % temporary variable 2
                cu2 = u(i,j)*cx(k) + v(i,j)*cy(k);
                
                % calculate equilibrium distribution function
                feq(i,j,k) = rho(i,j)*w(k)*(1+3*cu2+4.5*cu2*cu2-1.5*cu1);
                
                % calculate updated distribution function
                fn05(i,j,k) = omega*feq(i,j,k)+(1-omega)*f(i,j,k);
            end
        end
    end
    
% ===== perform streaming of distribution function ========================
    
    % stream on interior first
    for i=2:nx-1
        for j=2:ny-1
            fn1(i,j,1) = fn05(i,j,1);       
            fn1(i,j,2) = fn05(i-1,j,2);   
            fn1(i,j,3) = fn05(i,j-1,3);   
            fn1(i,j,4) = fn05(i+1,j,4);   
            fn1(i,j,5) = fn05(i,j+1,5);   
            fn1(i,j,6) = fn05(i-1,j-1,6);   
            fn1(i,j,7) = fn05(i+1,j-1,7);   
            fn1(i,j,8) = fn05(i+1,j+1,8);   
            fn1(i,j,9) = fn05(i-1,j+1,9);   
        end
    end
    
    % stream on north boundary
    for i=2:nx-1
        j = ny;
        fn1(i,j,1) = fn05(i,j,1);       
        fn1(i,j,2) = fn05(i-1,j,2);   
        fn1(i,j,3) = fn05(i,j-1,3);   
        fn1(i,j,4) = fn05(i+1,j,4);    
        fn1(i,j,6) = fn05(i-1,j-1,6);   
        fn1(i,j,7) = fn05(i+1,j-1,7);    
    end
    
    % stream on south boundary
    for i=2:nx-1
        j=1;
        fn1(i,j,1) = fn05(i,j,1);       
        fn1(i,j,2) = fn05(i-1,j,2);    
        fn1(i,j,4) = fn05(i+1,j,4);   
        fn1(i,j,5) = fn05(i,j+1,5);     
        fn1(i,j,8) = fn05(i+1,j+1,8);   
        fn1(i,j,9) = fn05(i-1,j+1,9);   
    end
    
    % stream on east boundary
    for j=2:ny-1
        i=nx;
        fn1(i,j,1) = fn05(i,j,1);       
        fn1(i,j,2) = fn05(i-1,j,2);   
        fn1(i,j,3) = fn05(i,j-1,3);    
        fn1(i,j,5) = fn05(i,j+1,5);   
        fn1(i,j,6) = fn05(i-1,j-1,6);      
        fn1(i,j,9) = fn05(i-1,j+1,9);   
    end
    
    % stream on west boundary
    for j=2:ny-1
        i=1;
        fn1(i,j,1) = fn05(i,j,1);        
        fn1(i,j,3) = fn05(i,j-1,3);   
        fn1(i,j,4) = fn05(i+1,j,4);   
        fn1(i,j,5) = fn05(i,j+1,5);    
        fn1(i,j,7) = fn05(i+1,j-1,7);   
        fn1(i,j,8) = fn05(i+1,j+1,8);      
    end
    
    % north-east corner
    i=nx; j=ny;
    fn1(i,j,1) = fn05(i,j,1);       
    fn1(i,j,2) = fn05(i-1,j,2);   
    fn1(i,j,3) = fn05(i,j-1,3);     
    fn1(i,j,6) = fn05(i-1,j-1,6);   
            
    % north-west corner
    i=1; j=ny;
    fn1(i,j,1) = fn05(i,j,1);          
    fn1(i,j,3) = fn05(i,j-1,3);   
    fn1(i,j,4) = fn05(i+1,j,4);     
    fn1(i,j,7) = fn05(i+1,j-1,7);   
            
    % south-east corner
    i=nx; j=1;
    fn1(i,j,1) = fn05(i,j,1);       
    fn1(i,j,2) = fn05(i-1,j,2);     
    fn1(i,j,5) = fn05(i,j+1,5);      
    fn1(i,j,9) = fn05(i-1,j+1,9); 
            
    % south-west corner
    i=1; j=1;
    fn1(i,j,1) = fn05(i,j,1);         
    fn1(i,j,4) = fn05(i+1,j,4);   
    fn1(i,j,5) = fn05(i,j+1,5);     
    fn1(i,j,8) = fn05(i+1,j+1,8);
    
% ===== apply Zou-He BC on boundary edges =================================

    % north boundary
    for i=2:nx-1
        j = ny;
        utemp = uref;
        vtemp = 0;       
        rtemp = (1/(1+vtemp))*( fn1(i,j,1)+fn1(i,j,2)+fn1(i,j,4)+2*(fn1(i,j,3)+fn1(i,j,7)+fn1(i,j,6)) );
        
        % microscopic boundary conditions
        fn1(i,j,5) = fn1(i,j,3)-(2/3)*rtemp*vtemp;
        fn1(i,j,8) = fn1(i,j,6)+0.5*(fn1(i,j,2)-fn1(i,j,4))-(1/6)*(rtemp*vtemp)-0.5*(rtemp*utemp);
        fn1(i,j,9) = fn1(i,j,7)+0.5*(fn1(i,j,4)-fn1(i,j,2))-(1/6)*(rtemp*vtemp)+0.5*(rtemp*utemp);
        
        % macroscopic boundary conditions
        rho(i,j) = rtemp;
        u(i,j)   = utemp;
        v(i,j)   = vtemp;
    end
    
    % south boundary
    for i=2:nx-1
        j = 1;
        utemp = 0;
        vtemp = 0;       
        rtemp = (1/(1-vtemp))*( fn1(i,j,1)+fn1(i,j,2)+fn1(i,j,4)+2*(fn1(i,j,5)+fn1(i,j,8)+fn1(i,j,9)) );
        
        % microscopic boundary conditions
        fn1(i,j,3) = fn1(i,j,5)+(2/3)*rtemp*vtemp;
        fn1(i,j,6) = fn1(i,j,8)-0.5*(fn1(i,j,2)-fn1(i,j,4))+(1/6)*(rtemp*vtemp)+0.5*(rtemp*utemp);
        fn1(i,j,7) = fn1(i,j,9)+0.5*(fn1(i,j,2)-fn1(i,j,4))+(1/6)*(rtemp*vtemp)-0.5*(rtemp*utemp); 
        
        % macroscopic boundary conditions
        rho(i,j) = rtemp;
        u(i,j)   = utemp;
        v(i,j)   = vtemp;
    end
    
    % east boundary
    for j=2:ny-1
        i = nx;
        utemp = 0;
        vtemp = 0;       
        rtemp = (1/(1+utemp))*( fn1(i,j,1)+fn1(i,j,3)+fn1(i,j,5)+2*(fn1(i,j,2)+fn1(i,j,6)+fn1(i,j,9)) );
        
        % microscopic boundary conditions
        fn1(i,j,4) = fn1(i,j,2)-(2/3)*rtemp*utemp;
        fn1(i,j,8) = fn1(i,j,6)+0.5*(fn1(i,j,3)-fn1(i,j,5))-(1/6)*(rtemp*utemp)-0.5*(rtemp*vtemp);
        fn1(i,j,7) = fn1(i,j,9)-0.5*(fn1(i,j,3)-fn1(i,j,5))-(1/6)*(rtemp*utemp)+0.5*(rtemp*vtemp);
        
        % macroscopic boundary conditions
        rho(i,j) = rtemp;
        u(i,j)   = utemp;
        v(i,j)   = vtemp;
    end
    
    % west boundary
    for j=2:ny-1
        i = 1;
        utemp = 0;
        vtemp = 0;       
        rtemp = (1/(1-utemp))*( fn1(i,j,1)+fn1(i,j,3)+fn1(i,j,5)+2*(fn1(i,j,4)+fn1(i,j,7)+fn1(i,j,8)) );
        
        % microscopic boundary conditions
        fn1(i,j,2) = fn1(i,j,4)+(2/3)*rtemp*utemp;
        fn1(i,j,6) = fn1(i,j,8)-0.5*(fn1(i,j,3)-fn1(i,j,5))+(1/6)*(rtemp*utemp)+0.5*(rtemp*vtemp);
        fn1(i,j,9) = fn1(i,j,7)+0.5*(fn1(i,j,3)-fn1(i,j,5))+(1/6)*(rtemp*utemp)-0.5*(rtemp*vtemp);
        
        % macroscopic boundary conditions
        rho(i,j) = rtemp;
        u(i,j)   = utemp;
        v(i,j)   = vtemp;
    end
    
% ===== apply bounce-back BC at corners ===================================  
 
    % north-east corner
    i=nx; j=ny;
    
        % microscopic boundary conditions
        fn1(i,j,4) = fn1(i,j,2);
        fn1(i,j,5) = fn1(i,j,3);
        fn1(i,j,8) = fn1(i,j,6);

        % macroscopic boundary conditions
        rho(i,j) = fn1(i,j,1)+fn1(i,j,2)+fn1(i,j,3)+fn1(i,j,4)+fn1(i,j,5)+fn1(i,j,6)+fn1(i,j,7)+fn1(i,j,8)+fn1(i,j,9);
        u(i,j)   = 0;
        v(i,j)   = 0;
    
    % north-west corner
    i=1; j=ny;
    
        % microscopic boundary conditions        
        fn1(i,j,5) = fn1(i,j,3);
        fn1(i,j,2) = fn1(i,j,4);
        fn1(i,j,9) = fn1(i,j,7);

        % macroscopic boundary conditions
        rho(i,j) = fn1(i,j,1)+fn1(i,j,2)+fn1(i,j,3)+fn1(i,j,4)+fn1(i,j,5)+fn1(i,j,6)+fn1(i,j,7)+fn1(i,j,8)+fn1(i,j,9);
        u(i,j)   = 0;
        v(i,j)   = 0;
    
    % south-east corner
    i=nx; j=1;
        
        % microscopic boundary conditions     
        fn1(i,j,4) = fn1(i,j,2);
        fn1(i,j,3) = fn1(i,j,5);
        fn1(i,j,7) = fn1(i,j,9);

        % macroscopic boundary conditions
        rho(i,j) = fn1(i,j,1)+fn1(i,j,2)+fn1(i,j,3)+fn1(i,j,4)+fn1(i,j,5)+fn1(i,j,6)+fn1(i,j,7)+fn1(i,j,8)+fn1(i,j,9);
        u(i,j)   = 0;
        v(i,j)   = 0;
    
    % south-west corner
    i=1; j=1;
        
        % microscopic boundary conditions     
        fn1(i,j,2) = fn1(i,j,4);
        fn1(i,j,3) = fn1(i,j,5);
        fn1(i,j,6) = fn1(i,j,8);

        % macroscopic boundary conditions
        rho(i,j) = fn1(i,j,1)+fn1(i,j,2)+fn1(i,j,3)+fn1(i,j,4)+fn1(i,j,5)+fn1(i,j,6)+fn1(i,j,7)+fn1(i,j,8)+fn1(i,j,9);
        u(i,j)   = 0;
        v(i,j)   = 0;
    
% ===== calculate macroscopic properties ==================================

    % Density and Velocity
    for i=2:nx-1
        for j=2:ny-1
            % temporary variables
            fsum  = 0;
            utemp = 0;
            vtemp = 0;
            
            for k=1:9
                fsum  = fsum  + fn1(i,j,k);
                utemp = utemp + fn1(i,j,k)*cx(k);
                vtemp = vtemp + fn1(i,j,k)*cy(k);
            end
            rho(i,j) = fsum;
            u(i,j)   = utemp/rho(i,j);
            v(i,j)   = vtemp/rho(i,j);
        end
    end    
    
% ===== plot solution (contour plot of velocity) ==========================
    if(pa==1)
        contourf(u'./uref)
        colorbar
        drawnow;
    end
    
% ===== calculate residuals ===============================================
    usum = 0;
    vsum = 0;
    rsum = 0;
    
    for i=1:nx
        for j=1:ny
            usum = usum + u(i,j);
            vsum = vsum + v(i,j);
            rsum = rsum + rho(i,j);
        end
    end
    
    % values for normalising residuals
    if(t==1)
        if(abs(usum-usum0)==0)
            unorm = 1;
        else
            unorm = abs(usum-usum0);
        end
        
        if(abs(vsum-vsum0)==0)
            vnorm = 1;
        else
            vnorm = abs(vsum-vsum0);
        end
        
        if(abs(rsum-rsum0)==0) 
            rnorm = 1;
        else
            rnorm = abs(rsum-rsum0);
        end
    end
    
    % calculate residuals
    ures = (abs(usum-usum0))/unorm;
    vres = (abs(vsum-vsum0))/vnorm;
    rres = (abs(rsum-rsum0))/rnorm; 
       
% ===== output iterations to screen =======================================
    if(mod(t,rf)==0)
        if(mod(t,100)==0 || t==1)
        fprintf('iteration---u residual-----v residual-----rho residual\n'); 
        end
        fprintf('%8d   %12.5e   %12.5e   %12.5e\n',t,ures,vres,rres);
    end
    
    % check if solution has converged
    if(ures<eps && vres<eps && rres<eps)
        break;
    end
    
% ===== increase counter ==================================================
    t = t+1;

end

%% Output data for postprocessing

fid=fopen('output.dat','w');
fprintf(fid,'TITLE="LBM SIMULATION"\n');
fprintf(fid,'VARIABLES = "x", "y", "u", "v", "p", "rho"\n');
fprintf(fid,'ZONE T="FLUID", I=%d, J=%d, F=POINT\n', nx, ny);
for j=1:ny
    for i=1:nx
        fprintf(fid,'%12.5e   %12.5e   %12.5e   %12.5e   %12.5e   %12.5e\n', x(i,j), y(i,j), u(i,j), v(i,j), rho(i,j)/3, rho(i,j));
    end
end
fclose(fid);

%  ===== print contour plot to screen =====================================
contourf(u'./uref)
colorbar












