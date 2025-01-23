import matplotlib.pyplot as plt
def somme_1(mu,a0,n_an):
    """Calcule p(1) en utilisant a0=1 (tout est proportionnel à a0)
    et a1=0 (pour avoir la condition p'(0)=0).
    Les paramètres sont m et mu.
    On choisira m et mu (ou l'un en fonction de l'autre)
    pour avoir p(1)=0..."""
    an = a0 #a0
    somme_p = an
    #an+2
    for n in range(0, n_an):
        n =n*2
        an = an *(n*(n + 1) - mu) / ((n + 2) * (n + 1))
        somme_p += an
    return somme_p
# liste : valeur de la somme pour des mu donnés
def list_somme(MU):
    somme_list = []
    for mu in MU :
        somme_list.append(somme_1(mu,a0,n_an))
    return somme_list
# Recherche (par parcours) des points où p(1) change de signe (donc s'annule)
def zeros_found(MU,somme_list) :
    # On les rajoutera en rouge, en plus de les afficher
    ListeMus = [] # Listes des endroit où on a convergence
    for i in range(len(MU) - 1):
        
        if somme_list[i] == 0:
            ListeMus.append(somme_list[i])
        if somme_list[i + 1] == 0:
            
            ListeMus.append(somme_list[i + 1])
        if somme_list[i] * somme_list[i + 1] < 0:
            mu = (MU[i]+MU[i + 1]) / 2
            ListeMus.append(mu)
            print("Valeur possible de mu :", mu, 'pour m =', 0)
    return ListeMus
#%% Initialisation
n_mu = 100000
MU = [i / 1000 for i in range(0, n_mu)]
n_an = 32000
a0 = 1
#%% Calcul de la somme en fonction des Mu et zéros
somme_p_list = list_somme(MU)
list_zeros = zeros_found(MU,somme_p_list)
#%% Affichage sommes en fonciton de mu
plt.figure(1, (15, 15))
plt.title(r'Graphe de p(1) en fonction de $\mu$ pour m='+str(0))
plt.plot(MU, somme_p_list, 'b')
plt.plot([0, MU[-1]], [0,0], 'k')
plt.plot(list_zeros, [0]*len(list_zeros), 'or')
plt.xlabel(r'$\mu$')
plt.ylabel('p(1)')
plt.show()