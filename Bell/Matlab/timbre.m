clear all;
%---------------------

%  Paramètres

%---------------------

N=90; % Nombre de points en theta

K=180; % Nombre de points en phi

n0=40; % numéro du point max en theta <=N

k0=100; % numéro du point max en phi <=K

u0_i=0.0001; % Amplitude au point d'impact à t0

Nt=60000; %nombre de points en temps

Tau=7.; %durée du signal en secondes

R=0.05; % rayon du timbre (m)

h=0.001; % épaisseur du timbre (m)

rho=8500. ; % masse volumique kg/m3

sigma=0;%0.001  % Amortissement en kg/m2/s

E=0; % 110*10^9; % module d'Young

nu=0.37; % coefficient de poisson

Gain=100. ; % Gain d'amplification sur le niveau émis

% Noms du fichier wav

Lenom_Sce='Son_Timbre_Test.wav';

%---------------------

%  Calcul des valeurs initiales et aux bornes

%---------------------

D=E*h^3/12/(1-nu^2); % Rigidité de flexion

delta_theta=pi/2/N; %pas en theta (radian)

delta_phi=2*pi/K; %pas en phi (radian)

delta_t=Tau/Nt; % pas en temps (secondes)

%

alpha=rho*h/delta_t^2+sigma/2/delta_t;

gamma=2*rho*h/delta_t^2;

zeta=-rho*h/delta_t^2+sigma/2/delta_t;

%

theta=[1:1:N]*delta_theta;% angle antilatitude

phi=[1:1:K]*delta_phi;% angle longitude

u1=zeros(N,K); % Elongations à l'instant t+1

u0=zeros(N,K);% Elongations à l'instant t

u_=zeros(N,K);% Elongations à l'instant t-1

Fs=60000000; % fréquence d'échantillonnage

d_dtheta=zeros(N,K);% Dérivée par rapport à theta

d2_dtheta2=zeros(N,K);% Dérivée seconde par rapport à theta

d2_dphi2=zeros(N,K);% Dérivée seconde par rapport à phi

d3_dtheta3=zeros(N,K);% Dérivée troisième par rapport à theta

d4_dtheta4=zeros(N,K);% Dérivée quatrième par rapport à theta

d4_dphi4=zeros(N,K);% Dérivée quatrième par rapport à phi

d3_dphi_dtheta2=zeros(N,K);% Dérivée troisième par rapport à phi et theta carré

d3_dtheta_dphi2=zeros(N,K);% Dérivée troisième par rapport à theta et phi carré

d4_dtheta2_dphi2=zeros(N,K);% Dérivée quatrième par rapport à theta carré et phi carré

nabla4=zeros(N,K);% Double laplacien

% Mise à zéro du niveau moyen

y_moy(1:Nt)=zeros(Nt,1);

 

% Forme de la corde et vitesse initiales

u0(n0,k0)=u0_i;

% Valeurs nulles juste avant et partout ailleurs

% Déjà mises à zéro dans l'initialisation

%---------------------

%---------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRANDE BOUCLE SUR LE TEMPS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:1:Nt

    t_courant=t*delta_t;

    disp(t_courant);

    y_moy(t)=Gain*mean(u0(:,:),'all');

%---------------------

% Calcul des dérivées successives

%---------------------

%

for k=1:1:K

    for n=2:1:N-1

d_dtheta(n,k)=1/2/delta_theta*(u0(n+1,k)-u0(n-1,k));

d2_dtheta2(n,k)=1/delta_theta^2*(u0(n+1,k)-2*u0(n,k)+u0(n-1,k));

    end

end

%

for k=2:1:K-1

    for n=1:1:N

d2_dphi2(n,k)=1/delta_phi^2*(u0(n,k+1)-2*u0(n,k)+u0(n,k-1));

    end

end

for k=K

    for n=1:1:N

d2_dphi2(n,k)=1/delta_phi^2*(u0(n,1)-2*u0(n,K)+u0(n,K-1));

    end

end

%

for k=1:1:K

    for n=3:1:N-2

d3_dtheta3(n,k)=3./8/delta_theta^3*(u0(n+2,k)-u0(n-2,k)-2*(u0(n+1,k)-u0(n-1,k)));

d4_dtheta4(n,k)=1/delta_theta^4*(u0(n+2,k)+u0(n-2,k)-4*u0(n+1,k)-4*u0(n-1,k) +6*u0(n,k));

    end

end

%

for k=3:1:K-2

    for n=1:1:N

d4_dphi4(n,k)=1/delta_phi^4*(u0(n,k+2)+u0(n,k-2)-4*u0(n,k+1)-4*u0(n,k-1)+6*u0(n,k));

    end

end

for k=1

    for n=1:1:N

d4_dphi4(n,k)=1/delta_phi^4*(u0(n,k+2)+u0(n,K-1)-4*u0(n,k+1)-4*u0(n,K) +6*u0(n,k));

    end

end

for k=2

    for n=1:1:N

d4_dphi4(n,k)=1/delta_phi^4*(u0(n,k+2)+u0(n,K)-4*u0(n,k+1)-4*u0(n,k-1)+6*u0(n,k));

    end

end

for k=K-1

    for n=1:1:N

d4_dphi4(n,k)=1/delta_phi^4*(u0(n,1)+u0(n,k-2)-4*u0(n,k+1)-4*u0(n,k-1) +6*u0(n,k));

    end

end

for k=K

    for n=1:1:N

d4_dphi4(n,k)=1/delta_phi^4*(u0(n,2)+u0(n,k-2)-4*u0(n,1)-4*u0(n,k-1) +6*u0(n,k));

    end

end

%

for k=2:1:K-1

    for n=2:1:N-1

d3_dphi_dtheta2(n,k)=1/2/delta_phi/delta_theta^2*(u0(n+1,k+1)-2*u0(n,k+1)   +u0(n-1,k+1)-u0(n+1,k-1)+2*u0(n,k-1)-u0(n-1,k-1));

d3_dtheta_dphi2(n,k)=1/2/delta_theta/delta_phi^2*(u0(n+1,k+1)-2*u0(n+1,k)  +u0(n+1,k-1)-u0(n-1,k+1)+2*u0(n-1,k)-u0(n-1,k-1));

d4_dtheta2_dphi2(n,k)=1/delta_theta^2/delta_phi^2*(u0(n+1,k+1)-2*u0(n+1,k)+u0(n+1,k-1)) -2/delta_theta^2/delta_phi^2*(u0(n,k+1)-2*u0(n,k)+u0(n,k-1)) +1/delta_theta^2/delta_phi^2*(u0(n-1,k+1)-2*u0(n-1,k)+u0(n-1,k-1));

    end

end

for k=1

    for n=2:1:N-1

d3_dphi_dtheta2(n,k)=1/2/delta_phi/delta_theta^2*(u0(n+1,k+1)-2*u0(n,k+1)+u0(n-1,k+1)-u0(n+1,K)+2*u0(n,K)-u0(n-1,K));

d3_dtheta_dphi2(n,k)=1/2/delta_theta/delta_phi^2*(u0(n+1,k+1)-2*u0(n+1,k)   +u0(n+1,K)-u0(n-1,k+1)+2*u0(n-1,k)-u0(n-1,K));

d4_dtheta2_dphi2(n,k)=1/delta_theta^2/delta_phi^2*(u0(n+1,k+1)-2*u0(n+1,k)+u0(n+1,K))  -2/delta_theta^2/delta_phi^2*(u0(n,k+1)-2*u0(n,k)+u0(n,K)) +1/delta_theta^2/delta_phi^2*(u0(n-1,k+1)-2*u0(n-1,k)+u0(n-1,K));

    end

end

for k=K

    for n=2:1:N-1

d3_dphi_dtheta2(n,k)=1/2/delta_phi/delta_theta^2*(u0(n+1,1)-2*u0(n,1) +u0(n-1,1)-u0(n+1,k-1)+2*u0(n,k-1)-u0(n-1,k-1));

d3_dtheta_dphi2(n,k)=1/2/delta_theta/delta_phi^2*(u0(n+1,1)-2*u0(n+1,k) +u0(n+1,k-1)-u0(n-1,1)+2*u0(n-1,k)-u0(n-1,k-1));

d4_dtheta2_dphi2(n,k)=1/delta_theta^2/delta_phi^2*(u0(n+1,1)-2*u0(n+1,k)+u0(n+1,k-1))  -2/delta_theta^2/delta_phi^2*(u0(n,1)-2*u0(n,k)+u0(n,k-1))  +1/delta_theta^2/delta_phi^2*(u0(n-1,1)-2*u0(n-1,k)+u0(n-1,k-1));

    end

end

%---------------------

% Calcul du double laplacien

%---------------------

for k=1:1:K

    for n=1:1:N

        nabla4(n,k)=d4_dtheta4(n,k)+2*cot(theta(n))*d3_dtheta3(n,k)+( (cos(theta(n)))^2-2 )/(sin(theta(n)))^2*d2_dtheta2(n,k)+ cot(theta(n))/(sin(theta(n)))^2*d_dtheta(n,k) +1/(sin(theta(n)))^4*d4_dphi4(n,k)+2*( (cos(theta(n)))^2+1 )/(sin(theta(n)))^4*d2_dphi2(n,k)+2/(sin(theta(n)))^2*d4_dtheta2_dphi2(n,k)-2*cot(theta(n))/(sin(theta(n)))^2*d3_dtheta_dphi2(n,k);

    end

end

%

% Calcul du vecteur à t+1

for k=1:1:K

    for n=1:1:N

        u1(n,k)=1/alpha*(-D/R^4*nabla4(n,k)+gamma*u0(n,k)+zeta*u_(n,k)); %vecteur à t+1

    end

end

%

% Passage à l'instant suivant

u_=u0;

u0=u1;

 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIN DE GRANDE BOUCLE SUR LE TEMPS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%plot(y_moy);

%mesh(u0);

% Ecriture du fichier de wav

audiowrite(Lenom_Sce,y_moy,Fs); %

 