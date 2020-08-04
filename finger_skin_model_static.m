close all

global km kt Kn kn B Ndx posi nu

disp('Initializing')
%% --------------------Parametrization--------------------
% % Skin parameters
E_ext = 1.54e6;             %% young's modulus (Pa) of the epidermis
d = 1e-2;                   %% diameter of the contact area
e_ext = 20e-6;              %% thickness of the external layer
e_epi = 1e-6;               %% thickness of the epidermis
E_tis = 0.025e6;            %% young's modulus (Pa) of the skin tissue
nu = 0.48;                  %% poisson's ratio (Fung)

% % Skin structure
r = 8e-3;                   %% finger radius in m
Ndx = 501;                  %% number of elements

%% ----------------Shape function of the finger----------------
d0 = 4e-4;                  %% initial distance from the surface
S = @(x) r-sqrt(r^2-x.^2);  %% shape function of the finger surface (half-circle)
dx = 1e-6;
x = -r+1e-6:dx:r-1e-6;
p = [x',S(x)'];
[q,dL] = curvspace(p,Ndx);  %% allows to have an equally spaced discretization along the curve
xi = q(:,1);                %% abscisses of the equally spaced points

bone_i = [r;0];             %% initial position of the bone

posi = [r;S(xi);0;xi];

%% ------------------Matrix of dependency------------------
    %%Stiffness matrix (spring)
km = E_ext*d*e_ext/dL/2;    %% stiffness of the membrane
kt = E_tis*pi*r^2/2/r;      %% stiffness of the underlying tissue
kn = 5e3;

Kn = sparse([ 1,2,Ndx+1,1, 2,1,     Ndx+1],...
            [ 1,1,1,    2, 2,Ndx+1, Ndx+1],...
            [-2,1,1,    1,-1,1,    -1],...
            2*(Ndx+1),2*(Ndx+1));

%% ---------------------------Contact---------------------------
kc = 100*kt;                 %% stiffness of the contact spring for penalty method

% % Friction model (Dahl)              
fclb = 0.6;                  %% coulomb friction coefficient
sigma0 = 5e2;                %% rest stiffness of the bristle (N/m)
n = 0.7;                     %% exponent: material dependent parameter
                                % 0<=n<=1 for brittle materials
                                
% % External load
P = 3;                       %% initial pression in N

    %%Damping matrix (dashpot)
zeta = 0.1;
B = -zeta*speye(2*(Ndx+1));

% % Temporal parameters
Fsmin = kc/(2*zeta);
Fs = 5e5;                    %% sampling frequency
dt = 1/Fs;
time = 0:dt:0.05;            %% time vector
Ndt = length(time);

%% -------------------------RK4 function:-------------------------
    % Uc is the current displacement vector
    % Fc the current force vector
    % tc the current instant
% % see the function RK4_f_U

%% --------------------------Resolution--------------------------
% % Initial conditions
U = zeros(2*(Ndx+1),length(time));  %% displacement U = [UN UT]
pos(1:2*(Ndx+1),1) = posi;
Upt = zeros(Ndx+1,length(time));
Upn = zeros(Ndx+1,length(time));
F = zeros(2*(Ndx+1),length(time));	%% force F = [FN FT]
F(:,1:Ndt) = repmat([-P;zeros(Ndx,1);0;zeros(Ndx,1)],1,Ndt);

contact_area = NaN(Ndx,Ndt);
mu = NaN(Ndx,Ndt);
% f = waitbar(0,'iterations','Name','Computation...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
t=2;
% % Impedance loading
disp('Computing')
while(t<=2 || abs(sum(F(1:Ndx+1,t-1),1))>1e-3)      %% Stop condition: global net force <10%
    % % Force resolution
    for i=1:Ndx+1
        if(t>2)
            Upt(i,t-1) = (U(i+Ndx+1,t-1)-U(i+Ndx+1,t-2))/(dt);    %% speed along x
        end
        if(pos(i,t-1) <0)
            contact_area(i-1,t) = pos(i+Ndx+1,t-1);
            % Normal force by penalty method
            Fr = kc*(abs(pos(i,t-1)));
            F(i,t) = F(i,t) + Fr;
            % Tangential force by Dahl model of friction
            F(i+Ndx+1,t) = F(i+Ndx+1,t) + F(i+Ndx+1,t-1) - dt*...
                    sigma0*Upt(i,t-1)*(abs(1-F(i+Ndx+1,t-1)/(fclb*Fr)*sign(Upt(i,t-1)))).^n*...
                    sign(1-F(i+Ndx+1,t-1)/(fclb*Fr)*sign(Upt(i,t-1)));
            if(i == floor(Ndx/2)+2)
                F(i+Ndx+1,t) = 0;
            end
            mu(i-1,t) = abs(F(i+Ndx+1,t)./F(i,t));
        else
            F(:,t) = F(:,t);
        end
    end
    
    % % Runge-Kutta4 : displacement resolution
    k1 = RK4_f_U(t*dt,U(:,t-1),F(:,t));
    k2 = RK4_f_U((t+0.5)*dt,U(:,t-1)+0.5*dt*k1,F(:,t));
    k3 = RK4_f_U((t+0.5)*dt,U(:,t-1)+0.5*dt*k2,F(:,t));
    k4 = RK4_f_U((t+1)*dt,U(:,t-1)+dt*k3,F(:,t));
    U(:,t) = U(:,t-1)+ dt/6*(k1+2*k2+2*k3+k4);
    
    % % Adding shape function -> final position of each element
    pos(:,t) = U(:,t) + posi;
    
    t = t+1;
end
tf=t-1;

%% --------------------------Plotting--------------------------
h = figure;
hold on
axis off
ndraw = 50;
vec_plot = floor((Ndx+1)/2)+1-50:floor((Ndx+1)/2)+1+50;
sz = 20;
for k=1:tf
    if (~mod(k-1,ndraw))
        clf;
        subplot 211
        axis off
        hold on
        area([xi(1) xi(end)],[0 0],-1e-4,'Facecolor','k')
        cd = [NaN;mu(vec_plot(1:end-1),k)];
        cd(isnan(cd)) = zeros(length(cd(isnan(cd))),1);
        scatter(pos(Ndx+1+vec_plot,k),pos(vec_plot,k),sz,cd,'filled');
        plot(pos(Ndx+1+vec_plot,k),pos(vec_plot,k),'-b');
        caxis([0 1])
        
        quiver(pos(Ndx+1+vec_plot,k),pos(vec_plot,k),...
            F(Ndx+1+vec_plot,k),F(vec_plot,k),0.5,'Color','b')
        
        xlim([-4e-3 4e-3])
        xlabel('x (m)')
        ylabel('y (m)')
        axis equal
        subplot 212
        hold on
        muP = area(pos(Ndx+1+vec_plot,k),-fclb.*F(vec_plot,k),-1e-4,'Facecolor',[0 0.5 0.5]);
        alpha(muP,0.2);
        plot(pos(Ndx+1+vec_plot,k),F(Ndx+1+vec_plot,k));
        
        xlabel('Position')
        ylabel('Tangential traction (N)')
        xlim([-4e-3 4e-3])
        ylim([-0.05 0.1])
        
        drawnow
    pause(1e-3);
    end
end

%% ----------------------Strains----------------------
vec_plot = floor((Ndx+1)/2)+1-50:2:floor((Ndx+1)/2)+1+50;
sz = 10;
% % Plot
figure; hold on
subplot(3,1,1)
axis off
hold on
area([xi(1) xi(end)],[0 0],-1e-4,'Facecolor','k')
plot(pos(Ndx+1+vec_plot,tf),pos(vec_plot,tf),'-b');
cd = [NaN;mu(vec_plot(1:end-1),tf)];
cd(isnan(cd)) = zeros(length(cd(isnan(cd))),1);
scatter(pos(Ndx+1+vec_plot,tf),pos(vec_plot,tf),sz,cd,'filled');

caxis([0 1])
quiver(pos(Ndx+1+vec_plot,tf),pos(vec_plot,tf),...
    F(Ndx+1+vec_plot,tf),F(vec_plot,tf),0.5,'Color','b')

subplot(3,1,2)
hold on
muP = area(pos(Ndx+1+vec_plot,tf),-fclb.*F(vec_plot,tf),-1e-4,'Facecolor',[0 0.5 0.5]);
alpha(muP,0.3);
plot(pos(Ndx+1+vec_plot(1:end),tf),F(Ndx+1+vec_plot,tf));
xlabel('Position (N)')
ylabel('Tangential traction')

subplot(3,1,3)
plot(pos(Ndx+1+vec_plot(1:end-1),tf).*1e3,diff(U(Ndx+1+vec_plot,tf))./diff(pos(Ndx+1+vec_plot,tf)));
xlabel('Position (mm)')
ylabel('Tangential strain')

%% Divergence computation
ac = find(isnan(contact_area(:,tf))==0);
for kk=1:tf
    div(:,kk) = 2*gradient(U(Ndx+3:end,kk),pos(Ndx+3:end,kk));
    div_median(:,kk) = nanmedian(div(ac,kk),1);
end

figure;
plot(abs(sum(F(2:Ndx+1,1:tf),1)),div_median);
xlabel('Normal force (N)')
ylabel('Divergence')
xlim([0 P])
ylim([-2e-3 0.05])
