print('TRABALHO FINAL, THIAGO HEIDI HORITA')
print('AGRADEÇO A DEUS, PELO MILAGRE')
print('AO PROFESSOR FABIO, QUE ADMIRO MUITO E A DAVID, UM IRMÃO')

                            #EULER EXPLICITO

#bibiliotecas
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt #graphics library

#criando a malha espacial
x_inicial = 0.0
x_final   = 100.0
N         = 100 # número de pontos da malha
Delta_x   = (x_final - x_inicial)/(N)
x         = np.linspace(x_inicial,x_final,N+1)

#criando a malha temporal
t_final           = 270.0 # tempo físico
Passo_Tempo_Total = 500 # total de pontos da malha temporal
Delta_t           = t_final/Passo_Tempo_Total

#propriedades especificas do modelo
cc     = 1.0 # dado já fornecido
alfa   = (cc*Delta_t)/(Delta_x**2) #de acordo com a demonstração algebrica
Lambda = 1.0 - 2*alfa # de acordo com a demonstração algebrica

# Verificando a condição CFL
CFL = (2*(cc*Delta_t))/Delta_x**2
print(CFL)

#Construção dos vetores da malha com variação Delta T e DElta X
u          = np.zeros(N+1) # declaração de u como vetor
u_ant      = np.zeros(N+1) # declaração de u como vetor
u0         = np.zeros(N+1) # declaração de u0 como vetor

# Condições de contorno de Dirichlet
u_left   = 30.0 #valores já dado pelo enunciado
u_right  = 67.0 #valores já dado pelo enunciado

# Condição inicial (C.I.)
for i in range(0,N+1):
    u0[i] = 30.0

u0[N] = u_right

# Inicializando os vetores u e u_ant com a C.I.
for i in range(0,N+1):
    u[i]          = u0[i]
    u_ant[i]      = u0[i]

# Solução numérica pelo método de Euler explícito
for i in range (1,Passo_Tempo_Total+1):
    for i in range(0,N):
        u_ant[i] = u[i] # atualização da solução no tempo
    # aplicando as condições de contorno
    u[1]   = alfa*u_left + Lambda*u_ant[1] + alfa*u_ant[2]
    for i in range(2,N-1):
        u[i] = alfa*u_ant[i-1] + Lambda*u_ant[i] + alfa*u_ant[i+1]
    u[N-1] = alfa*u_ant[N-2] + Lambda*u_ant[N-1] + alfa*u_right

#Gráficos
plt.plot(x,u0,'-*',x,u,'-*')
plt.legend(['Condicao Inicial','Tempo Final'])
plt.title('Equação de difusão de calor (Explícito)')
plt.show()

# em cores

u_graph = np.array([u])
heat_map = sb.heatmap(u_graph, cmap="rainbow", xticklabels = False, yticklabels = False)
plt.title('        Equação de difusão de calor (Explícito)')
plt.xlabel("Gradação visual de temperatura (°C)")
plt.show()
