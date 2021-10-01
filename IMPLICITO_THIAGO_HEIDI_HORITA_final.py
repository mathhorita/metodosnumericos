print('TRABALHO FINAL, THIAGO HEIDI HORITA')
print('AGRADEÇO A DEUS, PELO MILAGRE')
print('AO PROFESSOR FABIO, QUE ADMIRO MUITO E A DAVID, UM IRMÃO')

                            #EULER IMPLICITO

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
Passo_Tempo_Total = 1 # total de pontos da malha temporal
Delta_t           = t_final/Passo_Tempo_Total

#propriedades especificas do modelo
cc     = 1.0 # dado já fornecido
alfa   = (cc*Delta_t)/(Delta_x**2) #de acordo com a demonstração algebrica
Lambda = 1.0 + 2*alfa # de acordo com a demonstração algebrica

# Verificando a condição CFL
CFL = (2*(cc*Delta_t))/Delta_x**2
print(CFL)

#Construção dos vetores da malha com variação Delta T e DElta X
u          = np.zeros(N+1) # declaração de u como vetor
u_ant      = np.zeros(N+1) # declaração de u como vetor
u_iter     = np.zeros(N+1) # declaração de u_post como vetor
u0         = np.zeros(N+1) # declaração de u0 como vetor

#construção dos vetore do sistema matricial Ax=b da tridiagonal
a          = np.zeros(N+1) # diagonal principal da matriz do sistema
d          = np.zeros(N+1) # diagonal secundaria superior da matriz do sistema
c          = np.zeros(N+1) # diagonal secundaria inferior da matriz do sistema
b          = np.zeros(N+1) # lado direito do sistema

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

# Construção da Matriz "A" do sistema tridiagonal
for i in range(0,N+1):
    a[i] = Lambda     # Diagonal principal, valores já demonstrado algebricamente
    d[i] = -alfa      # Diagonal secundário superior, valores já demonstrado algebricamente
    c[i] = -alfa      # diagonal secundária inferior, valores já demonstrado algebricamente

# Construção do LADO DIREITO "B" do sistema tridiagonal
b[1]   = u_ant[1] + alfa*u_left
for i in range(2,N-1):
    b[i] = u_ant[i]
b[N-1] = u_ant[N-1] + alfa*u_right

# EVOLUÇÃO NO TEMPO
for t in range(1, Passo_Tempo_Total + 1):

    # Solução do sistema de equações lineares
    for i in range(0, N - 1):
        u_iter[i] = 1.0

    while (max(abs(u - u_iter)) > 0.001):

        # Atualização do Gauss Jacobi
        for i in range(0, N + 1):
            u_iter[i] = u[i]

        u[1] = (1 / a[1]) * (b[1] - d[1] * u_iter[2])
        for i in range(2, N - 1):
            u[i] = (1 / a[i]) * (b[i] - c[i - 1] * u_iter[i - 1] - d[i] * u_iter[i + 1])
        u[N - 1] = (1 / a[N - 1]) * (b[N - 1] - c[N - 2] * u[N - 2])

    # ATUALIZAÇÃO da Evolução no Tempo
    for i in range(0, N + 1):
        u_ant[i] = u[i]

    # Construção do NOVO LADO DIREITO do sistema
    b[1] = u_ant[1] + alfa * u_left
    for i in range(2, N - 1):
        b[i] = u_ant[i]
    b[N - 1] = u_ant[N - 1] + alfa * u_right

#gráficos
plt.plot(x, u0, '-*', x, u, '-*')
plt.legend(['Condicao Inicial', 'Tempo Final'])
plt.title('Equação de difusão de calor (Implícito)')
plt.show()

# em cores
u_graph = np.array([u])
heat_map = sb.heatmap(u_graph, cmap="rainbow", xticklabels=False, yticklabels=False)
plt.title('        Equação de difusão de calor (Implícito)')
plt.xlabel("Gradação visual de temperatura (°C)")
plt.show()

