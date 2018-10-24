
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as ticker 
import pickle 


def pcc(X, Y):
   ''' Compute Pearson Correlation Coefficient. '''
   # Normalise X and Y
   X -= X.mean(0)
   Y -= Y.mean(0)
   # Standardise X and Y
   X /= X.std(0)
   Y /= Y.std(0)
   # Compute mean product
   return np.mean(X*Y)

#parameters
Tau_1 = 500
Tau_2 = 2500
# TAU_2 =np.int64(np.linspace(100,8000,5))
Tau_3 = 500
Tau_4 = 2000.
# TAU_4 = np.int64(np.linspace(500,5000,15))
Tau_5 = 100
timestep = 10000
noise_std_1 = 0.001
noise_std_2 = 0.001
noise_std_3 = 0.001

Pi_1 = 1.
Pi_3 = 1.
Pi_2 = np.linspace(0.,2.,10)
z_threshold = 2000

e_1 = np.zeros((timestep,1),dtype = np.float32)
x = np.zeros((timestep,1),dtype = np.float32)
z = np.zeros((timestep,1),dtype = np.float32)
u = np.zeros((timestep,1),dtype = np.float32)
e_2 = np.zeros((timestep,1),dtype = np.float32)
v = np.zeros((timestep,1),dtype = np.float32)


Rho = np.zeros((10,10),dtype=  np.float32)

Z = np.random.normal(0.75,.25,50)

for m in range(10):
	for n in range(10):
		print('process:'+str(m))
		# z_threshold = Z_threshold[m]
		# Tau_2 = TAU_2[m]
		U = []
		V = []
		for l  in range(50):
			flag = 0
			e_1 = np.zeros((timestep,1),dtype = np.float32)
			x = np.zeros((timestep,1),dtype = np.float32)
			z = np.zeros((timestep,1),dtype = np.float32)
			u = np.zeros((timestep,1),dtype = np.float32)
			e_2 = np.zeros((timestep,1),dtype = np.float32)
			v = np.zeros((timestep,1),dtype = np.float32)
			z[:np.int64(timestep*.25),0] = Z[l]
			for i in range(1,timestep):	
				e_1[i,0] = np.abs(x[i,0] - z[i,0])
				u[i,0] = (1 - 1./Tau_1)*u[i-1] + Pi_1*e_1[i-1,0]*1./Tau_1 + np.random.normal(0,noise_std_1)
				if i > Tau_2:
					e_2[i,0] = u[i-Tau_2,0] 
					v[i,0] = (1 - 1./Tau_3)*v[i-1] + Pi_2[m]*e_2[i-1,0]*1./Tau_3 + Pi_3*z[i-1,0]*1./Tau_3+np.random.normal(0,noise_std_2)
				else:
					v[i,0] = (1 - 1./Tau_3)*v[i-1] + Pi_3*z[i-1,0]*1./Tau_3 + np.random.normal(0,noise_std_2)
				if np.sum(z[:i,:]) < z_threshold:
					if i < Tau_5:
						z[i,0] = z[i-1,0]*(1 - 1./Tau_4) + np.random.normal(0,noise_std_3)
					else:
						z[i,0] = z[i-1,0]*(1 - 1./Tau_4) + 1./Tau_4*x[i-Tau_5,0] + np.random.normal(0,noise_std_3)
				else:
					flag = 1
					z[i,0] = 0
			if flag > 0:
				U.append(np.sum(u[np.where(z>0)]))
				V.append(np.sum(v[np.where(z<0.01)]))
		U_ = np.asarray(U)
		V_  = np.asarray(V)
		Rho[m,n] = pcc(U_,V_)
# print(pcc(U_,V_))
# j  = np.max(np.where(z>0))
# t = (np.arange(timestep)-j)/1000.
# plt.plot(t,u,c = 'b',label = 'u')
# plt.plot(t,v,c = 'r',label = 'v')
# plt.plot(t,z,c = 'g',label = 'z')
# plt.plot(t,x, c= 'm',label = 'x')
# plt.legend()
# plt.show()
# ticks_x = ticker.FuncFormatter(lambda t, pos: '{0:g}'.format(t))
# plt.xlabel("seconds")

# plt.legend()
	# t = str(x_) + '.png'
	# plt.savefig(t,dpi=fig.dpi)`

# fs = plt.scatter(U_,V_,c =Z,cmap = 'rainbow')
# plt.xlabel('pre S1 '+r'$\gamma$'+' power')
# plt.ylabel('post ACC '+r'$\beta$'+' power')
# plt.show()

# cbar = plt.colorbar(fs)
# cbar.ax.set_ylabel('z(0)')
# plt.tight_layout()
with open('Rho7-3.pkl','wb') as f:
	pickle.dump([Rho],f,protocol =2)

