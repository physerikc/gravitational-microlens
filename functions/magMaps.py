import numpy as np
import matplotlib.pyplot as plt

num_pontos = 50000
cart_pontos = 200

def make_complex(x, y):
    return x + 1j*y
    
def pos_final(zi, n, r):
    
    if n == 1:
        e1 = 1
        e2 = 0
        e3 = 1 - e1 - e2
    elif n == 2:
        e1 = 0.5
        e2 = 0.5
        e3 = 1 - e1 - e2
    elif n == 3:
        e1 = 0.333
        e2 = 0.333
        e3 = 1 - e1 - e2
        
    soma = 0
    e = np.array([e1, e2, e3])
    
    for j in range(1, n+1):
        z_r = zi - r[j-1]
        soma += (e[j-1]*(z_r))/abs(z_r)**2
        z = zi - soma
    return z
        
def zs_final_position(Zp, n, r):
    zs = []
    
    for zi in Zp:
        z = pos_final(zi, n, r)
        zs.append(z)
    x = np.real(zs)
    y = np.imag(zs)
    
    return x, y

def arquimedes(num_points):
    
    # Definindo o ângulo de ouro
    golden_angle = np.pi * (3 - np.sqrt(5))
    
    # Cálculo dos ângulos usando o ângulo de ouro
    theta = golden_angle * np.arange(num_points)
    
    # Raio
    r = np.sqrt(theta)
    
    # Coordenadas
    x = r * np.cos(theta)/100
    y = r * np.sin(theta)/100

    return x, y

def plot_func(xp, yp, title):
    plt.figure(figsize=(6, 6))
    
    plt.scatter(xp, yp, s = 1, color='k', alpha=0.5)  # s é o tamanho dos pontos, alpha é a transparência

    plt.title(title)
    plt.xlabel(r'$x/R_{E}$')
    plt.ylabel(r'$y/R_{E}$')

    plt.axis('equal')
    plt.xlim(-1.5, 1.5)
    plt.ylim(-1.5, 1.5)

    plt.show()

# def make_complex(xp, yp):
#     # Vectorize a função para aplicá-la a todos os elementos de X e Y
#     make_complex_vec = np.vectorize(make_complex)  
    
#     # Crie o array Z de números complexos
#     zp = make_complex_vec(xp, yp)
#     return zp

def sub_plots(xr, yr, xc, yc, xa, ya):
    # Criar a figura e os subplots
    fig, axs = plt.subplots(1, 3, figsize=(11, 4))
    
    # Plotar cada gráfico nos subplots correspondentes
    # axs[0].plot_func(Xp_rand, Yp_rand, title)
    axs[0].scatter(xr, yr, s = 1, color='k', alpha=0.5)  # s é o tamanho dos pontos, alpha é a transparência
    axs[0].set_title('Distribuição randômica')
    axs[0].set_xlabel(r'$x/R_{E}$')
    axs[0].set_ylabel(r'$y/R_{E}$')
    axs[0].set_aspect('equal')
    axs[0].set_xlim(-1.5, 1.5)
    axs[0].set_ylim(-1.5, 1.5)

    # pts = np.column_stack([xc.flatten(), yc.flatten()])
    # plot_func(points[:, 0], points[:, 1], title)

    axs[1].scatter(xc, yc, s = 1, color='k', alpha=0.5)  # s é o tamanho dos pontos, alpha é a transparência
    axs[1].set_title('Distribuição cartesiana')
    axs[1].set_xlabel(r'$x/R_{E}$')
    axs[1].set_aspect('equal')
    axs[1].set_xlim(-1.5, 1.5)
    axs[1].set_ylim(-1.5, 1.5)

    axs[2].scatter(xa, ya, s = 1, color='k', alpha=0.5)  # s é o tamanho dos pontos, alpha é a transparência
    axs[2].set_title('Distribuição quase regular')
    axs[2].set_xlabel(r'$x/R_{E}$')
    axs[2].set_aspect('equal')
    axs[2].set_xlim(-1.5, 1.5)
    axs[2].set_ylim(-1.5, 1.5)
    # Ajustar o layout para evitar sobreposição
    plt.tight_layout()
    
    # Mostrar os gráficos
    plt.show()
