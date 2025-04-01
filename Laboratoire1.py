import numpy as np

def MDF(m, n, eps_r1, eps_r2, d, w, tol):
    V = np.zeros((n+1, m+1))
    V[:, 0] = 0  # Bord à 0V
    V[:, -1] = 0
    V[0, :] = 0
    V[-1, :] = 0
    V[4, 2]=V[4, 3]=V[4, 4]=V[4, 5]=1

    # Convergence par relaxation
    diff = tol + 1
    while diff > tol:
        V_old = V.copy()
        for i in range(1, n):
            for j in range(1, m):
                if (((i==4) & (j==1)) | ((i==4) & (j==6))):
                    V[i, j] = 0.25*(V[i, j+1] + V[i, j-1])+ (0.5/(eps_r1*eps_r2))*(eps_r1*V[i+1, j] + eps_r2*V[i-1, j])
                else:
                    V[i, j] = 0.25 * (V[i+1, j] + V[i-1, j] + V[i, j+1] + V[i, j-1])
                V[4, 2]=V[4, 3]=V[4, 4]=V[4, 5]=1
        diff = np.max(np.abs(V - V_old))

    return V

def IGauss(V, eps_r1, eps_r2, d, w):
    """
    Calcule la capacité en utilisant l'intégrale de Gauss.
    
    Paramètres :
    - V : matrice des potentiels
    - eps_r1, eps_r2 : constantes diélectriques des matériaux
    - d, w : dimensions physiques (non utilisées ici mais laissées pour cohérence)
    
    Retourne :
    - C : capacité calculée
    """
    n, m = V.shape
   
  
    Q = eps_r1*(8.85e-12)*(np.sum(V[1, :])+np.sum( V[:, 1])+np.sum( V[-2, :])+np.sum( V[:, -2]))
    V_max = np.max(V)
    if V_max == 0 or np.isnan(V_max):
        raise ValueError("V_max est nul ou invalide, vérifiez les calculs de V.")

    C = Q / V_max
    return C

def MicroPar(m, n, eps_r1, eps_r2, d, w, tol):
    """
    Calcule les paramètres Zo et vp d'une ligne de transmission.
    
    Paramètres :
    - m, n : dimensions de la grille
    - eps_r1, eps_r2 : constantes diélectriques des matériaux
    - d, w : dimensions physiques
    - tol : tolérance pour la convergence
    
    Retourne :
    - Zo : impédance caractéristique
    - vp : vitesse de propagation
    """
    # Calculer les potentiels en utilisant la méthode des différences finies
    V = MDF(m, n, eps_r1, eps_r2, d, w, tol)
    
    # Calculer la capacitance avec les constantes diélectriques
    C = IGauss(V, eps_r1, eps_r2, d, w)
    
    # Calculer la capacitance dans le vide (eps_r1 = eps_r2 = 1)
    C0 = IGauss(V, 1, 1, d, w)
    
    # Calculer Zo et vp
    epsilon=C/C0
    vp = 3e8 / epsilon**0.5
    Zo = 1/(vp*C)
    
    return Zo, vp

v0=MDF(7, 6, 1, 1, 2, 3, 0.002)
print(f"V0={v0}")
c0=IGauss(v0, 1, 1, 2, 3)
print(f"c0={c0}")
VZ0=MicroPar(7, 6, 1, 1, 2, 3, 0.002)
print(f"VZ0={VZ0}")
v1=MDF(7, 6, 10, 1, 2, 3, 0.002)
print(f"v1={v1}")
c1=IGauss(v1, 10, 1, 2, 3)
print(f"c1={c1}")
VZ1=MicroPar(7, 6, 10, 1, 2, 3, 0.002)
print(f"VZ1={VZ1}")

# Définition des paramètres fixes
m, n = 100, 80  # Domaine de définition
eps_r1, eps_r2 = 12, 1  # Constantes diélectriques
tol = 1e-4  # Tolérance 0.01% (changer pour tester 0.002%)

# Listes des valeurs à tester
d_values = [2, 3, 4, 6, 10, 20]
w_values = [3, 4, 6, 10, 16, 20, 30, 40]

w_d_ratios = []
Z0_values = []
vp_values = []
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
# Étude pour w = 2 et d variable
w_fixed = 2
for d in d_values:
    w_d_ratios.append(w_fixed / d)
    Zo, vp = MicroPar(m, n, eps_r1, eps_r2, d, w_fixed, tol)
    Z0_values.append(Zo)
    vp_values.append(vp)

# Étude pour d = 2 et w variable
d_fixed = 2
for w in w_values:
    w_d_ratios.append(w / d_fixed)
    Zo, vp = MicroPar(m, n, eps_r1, eps_r2, d_fixed, w, tol)
    Z0_values.append(Zo)
    vp_values.append(vp)

# Tracé des courbes
plt.figure(figsize=(10, 5))

# Graphique de Z0
plt.subplot(1, 2, 1)
plt.semilogx(w_d_ratios, Z0_values, 'o-', label="Z0")
plt.xlabel("w/d (échelle logarithmique)")
plt.ylabel("Impédance caractéristique Zo")
plt.title("Zo en fonction de w/d")
plt.grid(True, which="both", linestyle="--")
plt.legend()

# Graphique de vp
plt.subplot(1, 2, 2)
plt.semilogx(w_d_ratios, vp_values, 's-', label="vp", color="red")
plt.xlabel("w/d (échelle logarithmique)")
plt.ylabel("Vitesse de propagation vp (m/s)")
plt.title("vp en fonction de w/d")
plt.grid(True, which="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.show()
