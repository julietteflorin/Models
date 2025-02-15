clear;
% MODE impulsion à n0
%---------------------
%  Initialisation
%---------------------
N=100; % Nombre de points en espace
K=1000; %nombre de points en temps
c=1; %vitesse du son en m/s
delta_t=0.01; %pas temporel en s
delta_x=0.01; %pas physique en m


T=K*delta_t;   %durée du changement de fréquence en secondes

gamma=(c*delta_t/delta_x)^2 % idéalement <=1.00;

u0=1.0; %valeur de l'impulsion
u=zeros(N,K);% nombre de points situées en 0
Freq=zeros(1,K);
freqpul=1/2*5


for n=1:1:N %initVal:step:endVal — Incremente index = n par la valeur = 1 sur chaque iterations jusqu'à la valeur suivante
    u(n,1)=0;
    u(n,2)=0;
end
%---------------------
Longueur=N*delta_x %m
Duree=K*delta_t %s
for k=1:1:K % conditions aux limites
        u(1,k)=u0*sin(k*delta_t*2*pi*freqpul);
        u(N,k)=0;
end
for k=2:1:K-1 
    for n=2:1:N-1
        u(n,k+1)=gamma*(u(n-1,k)+u(n+1,k))+2*(1-gamma)*u(n,k)-u(n,k-1);
 
    end
end
%plot(u(1:N,1));
mesh(u');
%waterfall(u(1:N,1:500)')
