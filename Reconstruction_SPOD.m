%% %%%%%% Reconstruction SPOD -- Tiré de Towne et al. (Jfm 2018)%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nx: dimension grille PIV en x;
% Ny: dimension grille PIV en y;
% Nt: dimension temporelle avec Nt=3580;
% Fs: fréquence d'échantillonnage avec Fs=400Hz;
% U_Runi et V_Runi: champs de vitesse de dimension: [Nx Ny Nt] et i=1,2,3;
% Ndim: dimension du problème: ici en 2D (U et V) donc Ndim=2;
% X: data concaténées avec dimension temporelle d'abord: [3*Nt Nx Ny Ndim];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A. Définition des paramètres et chargement des data à effectuer: %%%%%%%
Nx=193;
Ny=67;
Nt=3580;
Ndim=2;
Fs=400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B. Réarrangement des data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tic
U_tild=cat(3,U_Run1,U_Run2,U_Run3); 
V_tild=cat(3,V_Run1,V_Run2,V_Run3);
X=cat(4,U_tild,V_tild);
X=permute(X,[3 1 2 4]); 
toc
clear U_tild V_tild  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C. Calcul de la SPOD (Towne et al. 2018): %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[L,P,F,Lc,Q,Theta,Lreal]=spod(X,ones(Nt,1),[],0,1/Fs);% ici 3 paquets de Nt 
toc
clear X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% D. Reconstruction des champs SPOD de la vitesse fluctuante: %%%%%%%%%%%%
n_block=1;% on choisit de reconstruire par exemple le premier bloc (le Run1)
n_modes=3;% on a choisi 3 paquets de Nt donc on a 3 modes
Nf=length(F);% longueur du vecteur des fréquences jusqu'à f_Nyquist=fs/2;
Q_rec=zeros(2*Nx*Ny, Nf, n_modes);% Q_reconstruction initialisé  

for n_freq=1:Nf     % on prend toutes les fréquences
    for n_mode=1:n_modes  % on prend tous les modes
        Psi_reshaped=reshape(squeeze(P(n_freq,:,:,:,n_mode)),[],1);
        Q_rec(:,n_freq,n_mode)=sqrt(n_modes*Lreal(n_freq,n_mode))*...
            Psi_reshaped*conj(Theta(n_block,n_mode,n_freq));
    end
    iteration_FREQ=n_freq
end

Q_rec = reshape(Q_rec(:), [Nx Ny Ndim Nf n_modes]);

Q_SUM = sum(Q_rec,5); % soit on vient sommer la contribution des 3 modes
% Q_SUM =squeeze(Q_rec(:,:,:,:,1));% soit on prend le mode 1 seulement! :) 

time=([0:(Nt-1)]/Nt);       % temps unitaire des Nt=3580 snapshots
time=permute(time,[3 1 2]); % le temps doit passer en dimension: 1x1x3
Flow_reco=zeros(Nx,Ny,Ndim,Nt);% on initialise la matrice des fluctuations
for t=1:Nt
    for m=1:Nf % on peut choisir ici l'intervalle de fréquence qu'on veut! 
        Flow_reco(:,:,:,t)=Flow_reco(:,:,:,t)+...
            real(Q_SUM(:,:,:,m)).*cos(2*pi*(m-1)*time(1,1,t))-...
            imag(Q_SUM(:,:,:,m)).*sin(2*pi*(m-1)*time(1,1,t));
    end
    iteration=t
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% E. On rajoute la composante moyenne des data concaténées: %%%%%%%%%%%%%%
U_Run1MOY=mean(U_Run1,3);
U_Run2MOY=mean(U_Run2,3);
U_Run3MOY=mean(U_Run3,3);
U_moy=1/3*(U_Run1MOY+U_Run2MOY+U_Run3MOY);%Moyenne U des data concaténées

V_Run1MOY=mean(V_Run1,3);
V_Run2MOY=mean(V_Run2,3);
V_Run3MOY=mean(V_Run3,3);
V_moy=1/3*(V_Run1MOY+V_Run2MOY+V_Run3MOY);%Moyenne V des data concaténées

% Puis on ajoute la partie moyenne aux fluctuations obtenues du bloc 1:
U_SPOD=U_moy+(squeeze(Flow_reco(:,:,1,:))); % matrice 3D [Nx Ny Nt]
V_SPOD=V_moy+(squeeze(Flow_reco(:,:,2,:))); % matrice 3D [Nx Ny Nt] 
%% F. Vérification numérique qu'on retombe bien sur la valeur brute de la 
% vitesse originale du Run 1 en sommant sur toutes les fréquences et les 3 
% modes SPOD obtenus: 
t=5; % je prends le 5ème snapshot au hasard

verif_U=U_Run1(:,:,t)-U_SPOD(:,:,t);
verif_V=V_Run1(:,:,t)-V_SPOD(:,:,t);

figure()
plot(reshape(verif_U,[],1),'k')
hold on
plot(reshape(verif_V,[],1),'r')
legend('écart \deltaU=orignal-SPOD','écart \deltaV=orignal-SPOD',...
       'location','northwest')
xlabel('Nombre N_x*N_y de points dans le domaine PIV')
ylabel('\deltaVitesse: (brute - SPOD_{reconstruite}) \sim 10^{-13}')
grid on
title('Écart \deltaV de l''ordre de 10^{-13}: la reconstruction SPOD marche! :)')
set(gca,'FontName','Times New Roman','FontSize',11)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
