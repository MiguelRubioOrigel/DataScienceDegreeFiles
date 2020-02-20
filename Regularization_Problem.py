import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from sklearn.model_selection import KFold

np.random.seed(100)

theta = [i+1 for i in range(10)]

print("Valores originales Theta")
print(theta)

m = 5
b = 2
mu = 0
sigma = 3
lamda = 100

x = np.arange(0,11,1)
X = [x, x**2, x**3, x**4, x**5, x**6]

xForPred = np.arange(0,10.1,.1)
XForPred = [xForPred, xForPred**2, xForPred**3, xForPred**4, xForPred**5, xForPred**6]

#-------------------------------Simulacion de Datos aleatorios

epsilon = np.random.normal(mu, sigma, (len(x)))
y = m * x + b + epsilon

plt.scatter(x,y)
plt.title("Datos Aleatorios Simulados - error Gauss(mu="+str(mu)+",sigma="+str(sigma)+")")
plt.show()

#-------------------------------Regresion Lineal

def l_function(Theta):
    return sum((y - np.dot(Theta[0:-1],X[0:1]) - Theta[-1])**2)

res = minimize(l_function, theta[0:2], method='BFGS')

ypred = res.x[0] * xForPred + res.x[1] 
plt.scatter(x,y)
plt.plot(xForPred,ypred)
plt.title("Regresion Lineal")
plt.show()

print("Valores de Theta en Ho - polinomio grado 1")
print(res.x)

#-------------------------------Hipotesis A
#-------------------------------Regresion Cuadratica

def ha_function(Theta):
    return sum((y - np.dot(Theta[0:-1],X[0:2]) - Theta[-1])**2)

resHa = minimize(ha_function, theta[0:3], method='BFGS')

ypredHa = np.dot(resHa.x[0:-1],XForPred[0:2]) + resHa.x[-1]
plt.scatter(x,y)
plt.plot(xForPred, ypredHa)
plt.title("Regresion Cuadratica")
plt.show()

print("Valores de Theta en Ha - polinomio grado 2")
print(resHa.x)

#-------------------------------Hipotesis B
#-------------------------------Regresion Cubica

def hb_function(Theta):
    #return sum((y - thetas[0]*x - thetas[1]*x**2 - thetas[2]*x**3 - thetas[3])**2)
    return sum((y - np.dot(Theta[0:-1],X[0:3]) - Theta[-1])**2)

resHb = minimize(hb_function, theta[0:4], method='BFGS')

ypredHb = np.dot(resHb.x[0:-1],XForPred[0:3]) + resHb.x[-1]
plt.scatter(x,y)
plt.plot(xForPred, ypredHb)
plt.title("Regresion Cubica")
plt.show()

print("Valores de Theta en Hb - polinomio grado 3")
print(resHb.x)

#-------------------------------Hipotesis C
#-------------------------------Regresion Grado Superior

def hc_function(Theta):
    return sum((y - np.dot(Theta[0:-1],X[0:6]) - Theta[-1])**2)

resHc = minimize(hc_function, theta[0:7], method='BFGS')

ypredHc = np.dot(resHc.x[0:-1],XForPred[0:6]) + resHc.x[-1]
plt.scatter(x,y)
plt.plot(xForPred, ypredHc)
plt.title("Regresion Grado 6")
plt.show()

print("Valores de Theta en Hc - polinomio grado 6")
print(resHc.x)

#-------------------------------Regularizacion Ridge

def ridge_function(Theta):
    return sum((y - np.dot(Theta[0:-1],X[0:6]) - Theta[-1])**2) + lamda*sum(Theta**2)

resRidge = minimize(ridge_function, theta[0:7], method='BFGS')

ypredRidge = np.dot(resRidge.x[0:-1],XForPred) + resRidge.x[-1]
plt.scatter(x,y)
plt.plot(xForPred, ypredHc,'b-')
plt.plot(xForPred,ypredRidge,'r-')
plt.title("Regularizacion Grado 6 - Ridge lambda="+str(lamda))
plt.show()

print("Valores de Theta Regularizados por Ridge - lambda="+str(lamda)+"\n")
print(resRidge.x)

#-------------------------------Regularizacion Lasso

def lasso_function(Theta):
    return sum((y - np.dot(Theta[0:-1],X) - Theta[-1])**2) + lamda*sum(abs(Theta))

resLasso = minimize(lasso_function, theta[0:7], method='BFGS')

ypredLasso = np.dot(resLasso.x[0:-1],XForPred) + resLasso.x[-1]
plt.scatter(x,y)
plt.plot(xForPred, ypredHc,'b-')
plt.plot(xForPred,ypredLasso,'r-')
plt.title("Regularizacion Grado 6 - Lasso lambda="+str(lamda))
plt.show()

print("Valores de Theta Regularizados por Lasso - lambda="+str(lamda)+"\n")
print(resLasso.x)

#-------------------------------Lambda optimo Ridge

lambda_candidates = [10.0**i for i in np.arange(-10,11,1)]
lambda_candidates = np.array(lambda_candidates)
resRidge = [0 for i in range(9)]

def ridge_function_CV(Theta):
  return sum((y - np.dot(Theta[0:-1],X_train) - Theta[-1])**2) + lamda*sum(Theta**2)

for i in range(len(lambda_candidates)):
  lamda = lambda_candidates[i]
  for k in range(2,3,1):
    print("k="+str(k))
    kf = KFold(n_splits=k)
    for train_index, test_index in kf.split(x):
      #print (train_index, test_index)
      x_train, x_test = x[train_index], x[test_index]
      y_train, y_test = y[train_index], y[test_index]

    print(x_train)
    X_train = [x_train, x_train**2, x_train**3, x_train**4, x_train**5, x_train**6]
    X_test  = [x_test,  x_test**2,  x_test**3,  x_test**4,  x_test**5,  x_test**6]

    #resRidge[k-2] = minimize(ridge_function_CV, theta[0:7], method='BFGS')
    if   k==2: resRidge[k-2] = minimize(ridge_function_CV, theta[0:6], method='BFGS')
    #elif k==3: resRidge[k-2] = minimize(ridge_function_CV, theta[0:6], method='BFGS')
    else:      resRidge[k-2] = minimize(ridge_function_CV, theta[0:7], method='BFGS')
    #print(resRidge[k-2])
    #ypredRidge_CV = np.dot(resRidge[k-2].x[0:-1],X_test) + resRidge[k-2].x[-1]
    #SSR[k-2] = sum((ypredRidge_CV-y_test)**2)    
  #SSR_lamda[i] = (1/9)*sum(SSR)

#lamda_optimo = lambda_candidates[np.argmin(SSR_lamda)]
