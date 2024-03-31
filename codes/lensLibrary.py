import pandas as pd
import numpy as np
import sympy as sp
from numpy import *
import math
import warnings
from numpy.polynomial import Polynomial
from sympy import *
import matplotlib.pyplot as plt

def plot_ticks():
    plt.figure(figsize = (8,5))     
    plt.minorticks_on()     
    plt.tick_params(axis='x', which='both', top = True, right = True, direction = 'in', length=5, width=1)     
    plt.tick_params(axis='y', which='both', top = True, right = True, direction = 'in', length=5, width=1)

def path(mu, alpha, tE, t0, pontosLC, area, tEx):
  raz = (tEx - t0)/tE

  Yp = math.tan(alpha)*math.cos(alpha)*raz + mu/math.cos(alpha)
  Xp = math.cos(alpha)*raz

  return np.array(Xp*area), np.array(Yp*area)


def produtorioG(zs, j, n):
  prdg = 1.0
  for i in range(0, n):
    if i != j:
      prdg *= zs[i]
  return prdg


def produtorioH(zs, j):
  return zs[j]


def GeH(lentes, n, e, zs, z):
  # secao para o array G
  G = Integer(0)

  for j in range(0, n):
    G += e[j]*produtorioG(zs, j, n)
  G_exp = expand(G)
  coffG = np.array(G_exp.as_poly(z).all_coeffs())

    # secao para o array H
  H = Integer(1)

  for j in range(0, n):
    H *= produtorioH(zs, j)

  H_exp = expand(H)
  coffH = np.array(H_exp.as_poly(z).all_coeffs())

    # consertando tamanho dos arrays
  if len(coffG) < n+1:
    num_zeros = n+1-len(coffG)
    coffG = np.pad(coffG, (num_zeros, 0), 'constant')
  coffG = coffG[::-1]

  if len(coffH) < n+1:
    num_zeros = n+1-len(coffH)
    coffH = np.pad(coffH, (num_zeros, 0), 'constant')
  coffH = coffH[::-1]
  return coffG, coffH


def ETA(G, H):

    n = 2
  # funcao para rearranjar os polinomios pelo grau
    def grau(expr):
        return max([0] + [p.as_poly(z).degree() for p in expr.atoms(Pow)])

    # criacao dos arrays de simbolos

    var_G = symbols('G:%d' % (n+1))
    G_pol = Polynomial(var_G)
    
    # var_G = np.array(symbols('G:%d' % (n+1))).tolist()
    var_H = symbols('H:%d' % (n+1))
    H_pol = Polynomial(var_H)

    # para cada valor de i, quero criar um array de tamanho 5 onde armazenarei G^i * H^n-i
    all_list = []

    for i in range(0, n+1):
        polinomio = (G_pol**i)*(H_pol**(n-i))
        
        lst = sorted(polinomio, key=grau)
        all_list.append(lst)
    
        # mapeamentos
    mapG = {}
    for simbolo, valor in zip(var_G, G):
        mapG[simbolo] = valor

    mapH = {}

    for simbolo, valor in zip(var_H, H):
        mapH[simbolo] = valor


    # atribuicao dos mapeamentos aos seus simbolos
    lst_with_data = []

    for lista in all_list:
        nova_lista = [elemento.subs(mapG) for elemento in lista]
        nova_lista = [elemento.subs(mapH) for elemento in nova_lista]
        lst_with_data.append(nova_lista)

    etas = np.array(lst_with_data)

    return etas


def somatoria(eta, dataX, W, l):

  n = 2
  cl = []
  soma = 0

  for i in range(0, n+1):
    if l == 0:
      soma += - eta[i][l]*W[i]
    elif l <= n**2:
      soma += eta[i][l-1]*dataX[i] - eta[i][l]*W[i]
    else:
      soma += eta[i][l-1]*dataX[i]

  cl.append(re(expand(soma)))
  # print(cl)
  return cl

z = symbols('z')

def grau(expr):
  return max([0] + [p.as_poly(z).degree() for p in expr.atoms(Pow)])


def jacobian(Xp, Yp, e, z1, z2, z3, eta, G, H, precisao, i, rs):

  n = 2
  w = complex(Xp[i], Yp[i])
  # w_data = np.array(w)


  w1, w2, w3 = rs[0] - w, rs[1] - w, rs[2] - w
  w1_conj, w2_conj, w3_conj = np.conj(w1), np.conj(w2), np.conj(w3)

  ws = ([w1, w2, w3])
  ws_conj = ([w1_conj, w2_conj, w3_conj])

  # criacao de arrays necessarios
  def produtorioV(ws_conj, j):
    prdg = 1.0
    for i in range(1, n+1):
      if i != j:
        prdg *= z - ws_conj[i-1]
    return prdg


  V = Integer(0)
  for j in range(1, n+1):
    V += e[j-1]*produtorioV(ws_conj, j)
  # organizando os polinomios das linhas de V
  valueV = collect(expand(V), z)
  coeffV = Poly(valueV, z).coeffs()[::-1]

  while len(coeffV) < n+1:
    coeffV = np.append(coeffV, 0)

  dataV = np.array(coeffV)
  X = Integer(1)
  for j in range(1, n+1):
    X *= z - ws_conj[j-1]
  # organizando os polinomios das linhas de X
  valueX = collect(expand(X), z)
  coeffX = Poly(valueX, z).coeffs()[::-1]
  # coeffX = Poly(X, z).coeffs()[::-1]

  while len(coeffX) < n+1:
    coeffX = np.append(coeffX, 0)

  dataX = np.array(coeffX)

  # calculo de W

  W0 = np.array(w*dataX[0] - dataV[0])
  W1 = np.array(w*dataX[1] - dataV[1])
  W2 = np.array(w*dataX[2] - dataV[2])
  # W3 = np.array(w*dataX[3] - dataV[3])

  W = np.array([W0, W1, W2])

  cl = []
  for l in range(0, n**2+2):
    cl.append(somatoria(eta, dataX, W, l))


  cl = np.array(cl)
  coeffs_cl = [elemento for subarray in cl for elemento in subarray]

  raizes_zs = np.roots(coeffs_cl)

  def soma_frac(e, raizes_zs, rs, i):

    soma = 0
    for j in range(0, n):
      soma += e[j]/raizes_zs[i]-rs[j]

    fracao = 1/(1-(abs(soma)**2))
    # print(fracao)
    # fracao = np.nan_to_num(fracao)

    return fracao

  MagTotal = 0
  for i in range(0, n**2+1):
    MagTotal += soma_frac(e, raizes_zs, rs, i)

  return MagTotal


def caustic(G, H, lincaustics, lentes, e):
  q = e[1]/e[1]
  s = np.arange(0.4, 2.4, 0.2)
  Xpc = (1/(1+q))*(s - ((1-q)/s))
  return Xpc

