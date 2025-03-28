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
                V[i, j] = 0.25 * (V[i+1, j] + V[i-1, j] + V[i, j+1] + V[i, j-1])
                V[4, 2]=V[4, 3]=V[4, 4]=V[4, 5]=1
        diff = np.max(np.abs(V - V_old))

    return V

def IGauss(V, eps_r1, eps_r2, d, w):
    n, m = V.shape
   
  
    Q = eps_r1*(8.85e-12)*(np.sum(V[1, :])+np.sum( V[:, 1])+np.sum( V[-2, :])+np.sum( V[:, -2]))
    V_max = np.max(V)
    if V_max == 0 or np.isnan(V_max):
        raise ValueError("V_max est nul ou invalide, vérifiez les calculs de V.")

    C = Q / V_max
    return C

def MicroPar(m, n, eps_r1, eps_r2, d, w, tol):
    # Calculer les potentiels en utilisant la méthode des différences finies
    V = MDF(m, n, eps_r1, eps_r2, d, w, tol)
    
    # Calculer la capacitance avec les constantes diélectriques
    C = IGauss(V, eps_r1, eps_r2, d, w)
    
    # Calculer la capacitance dans le vide (eps_r1 = eps_r2 = 1)
    C0 = IGauss(V, 1, 1, d, w)
    
    # Calculer Zo et vp
    Zo = 1 / (3e8 * (C / C0)**0.5)
    vp = 3e8 / (eps_r1**0.5 if eps_r1 == eps_r2 else ((eps_r1 + eps_r2) / 2)**0.5)
    
    return Zo, vp

v0=MDF(7, 6, 1, 1, 2, 3, 0.002)
print(v0)
c0=IGauss(v0, 1, 1, 2, 3)
print(c0)
v1=MDF(7, 6, 10, 1, 2, 3, 0.002)
print(v0)
c1=IGauss(v0, 10, 1, 2, 3)
print(c0)
