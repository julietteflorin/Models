clear all;

%---------------------

%  Paramètres

%---------------------

K=10 ; % nombre de fréquences

f=zeros(K,1); % Tableau des fréquences mis à zéro

R=0.015; % rayon du timbre (m)

h=0.0008; % épaisseur du timbre (m)

rho=8470. ; % masse volumique kg/m3

sigma=10;%0.001  % Amortissement en kg/m2/s

E= 90*10^9; % module d'Young (kg m?1 s?2)(N/m2)

nu=0.37; % coefficient de poisson

%----------------------------------------------

%   Choix des paramètres ajustables

%----------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m0=2;% M0 C EST UN PARAMETRE
% Mise à zéro du niveau moyen
N=90; % Nombre de points en theta
    
K=360; % Nombre de points en phi

h0=1;

Tau=1.*4;   % durée en secondes

Nt=441000*4; %nombre de points en temps

Gain=0.1; % Gain d'amplification ou atténuation sur le niveau émis

% Noms du fichier wav

Lenom_Sce='Son_Timbre_Test2.wav';

%---------------------

%  Calcul des valeurs initiales et aux bornes

%---------------------

D=E*h^3/12/(1-nu^2) % Rigidité de flexion (kg m2 s?2)
%

gamma=sigma/2/rho/h  % atténuation

alpha=D/rho/h/R^4 

%

Fs=Nt/Tau; % fréquence d'échantillonnage



y_moy(1:Nt)=zeros(Nt,1);
for n0=0:1:6
    %n0=0;% N0 C EST UN PARAMETRE
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %-----------------------
    
    % Calculs des grandeurs
    
    %-----------------------
    
    delta_theta=pi/2/N; %pas en theta (radian)
    
    delta_phi=2*pi/K; %pas en phi (radian)
    
    delta_t=Tau/Nt; % pas en temps (s)
    
    theta=[1:1:N]*delta_theta;% angle antilatitude
    
    phi=[1:1:K]*delta_phi;% angle longitude
    
    P_theta=zeros(N,1); %Initialisation des valeurs à zéro
    
    W_theta=zeros(N,1); %Initialisation des valeurs à zéro
    
    
    h_m0=zeros(K,1); %Initialisation des valeurs à zéro
    
    u0=zeros(K,N); %Initialisation des valeurs à zéro
    
    k0=(n0+m0);
    
    Omega0_2=k0*(k0+1);
    
    f_k0=1/2/pi*sqrt(alpha*k0^2*(k0+1)^2-gamma^2)
    
    a=zeros(n0+3,1); %Initialisation des coefficients à zéro
    
    % Test de parité de n0
    
    if mod(n0,2)==0
    
        a(0+1)=1;
    
        for n=0:2:n0
    
            a(n+2+1)=a(n+1)*(n^2+n*(2*m0+1)+m0*(m0+1)-Omega0_2)/(n+2)/(n+1);
    
        end
    
    else
    

    
        a(1+1)=1;
    
        for n=1:2:n0
    
            a(n+2+1)=a(n+1)*(n^2+n*(2*m0+1)+m0*(m0+1)-Omega0_2)/(n+2)/(n+1);
    
        end
    
    end
    

    
    %-------------------------------
    
    % Calcul de p(theta)
    
    %-------------------------------
    
    for l=1:1:N
    
        P_theta(l)=0;
    
        for n=0:1:n0
    
            P_theta(l)=P_theta(l)+a(n+1)*(cos(theta(l)))^n;
    
        end
    
        W_theta(l)=(sin(theta(l)))^(m0/2)*P_theta(l);
    
    end
    
    
    for k=1:1:K
    
        h_m0(k)=h0*cos(m0*phi(k));
    
    end
    
    plot(h_m0);
    
    Energie=0;
    
    for l=1:1:N
    
        for k=1:1:K
    
            u0(k,l)=h_m0(k)*W_theta(l);
    
            Energie=Energie+(u0(k,l))^2/K/N;
    
        end
    
    end
    
    mesh(u0);
    
    %
    
    % Calcul du niveau moyen
    
    for t=1:1:Nt
    
        y_moy(t)= y_moy(t) + Gain*sqrt(Energie)*exp(-gamma*t*delta_t)*cos(2*pi*f_k0*t*delta_t);
    
    end
    
end
% Ecriture du fichier de wav

audiowrite(Lenom_Sce,y_moy,Fs); %