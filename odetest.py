from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def pend(y, t, beta, alpha):

	dydt = np.array(beta - alpha * y)

	return dydt
	
beta = 2
alpha = 5

y0 = np.array([0])
t = np.linspace(0, 0.5, 101)			#(start, finish, number of points in between)
sol = odeint(pend, y0, t, args=(beta,alpha))


#[:,0] refers to [row,column]
#so plot the column values against time. 0 refers to first column, 1 to second. etc
plt.plot(t, sol[:, 0], 'b', label='0')	
plt.legend(loc='best')					
plt.xlabel('t')
plt.grid()
plt.show()
