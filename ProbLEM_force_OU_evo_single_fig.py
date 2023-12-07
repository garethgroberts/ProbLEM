import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition
import os

# gareth roberts, 12 may 2023, caltech & imperial college london
f=open('rivevo.csv','w')
f=open('rivevo.csv','a')
f2=open('rivevo_exp_var.csv','w')
f2=open('rivevo_exp_var.csv','a')
f3=open('rivevo_early.csv','w')
f3=open('rivevo_early.csv','a')
f4=open('xprob.csv','w')
f4=open('xprob.csv','a')
f5=open('xforce.csv','w')
f5=open('xforce.csv','a')
f6=open('txforce.csv','w')
f6=open('txforce.csv','a')
f7=open('txexpd.csv','w')
f7=open('txexpd.csv','a')

print("=============================================================================================")
N = 1000 		#number of blocks (elements in spatial array) along river
dx = 1			#width of each block
dt = 1			#time step length
times = np.arange(0, N, dt)

# critical value(s) above which erosion occurs
cmean = 2					#mean
cstd = 0					#standard deviation 
cmean = np.full(N,cmean)	#define along entire profile
cstd = np.full(N,cstd)		#define along entire profile
cvar = cstd**2				#variance

x = np.arange(0,N,dx)		#set up spatial grid
z = np.arange(N,0,-1)	    #set up initial elevation without dz

#assume force is determined by Ornstein-Uhlenbeck process
theta=1e-3	#drift/reversion coefficient
sigma=2e-2	#volatility
mu=3		#long term mean
xsize = np.size(x) 
dw = np.random.normal(0, scale = np.sqrt(dx), size = xsize)
force = np.zeros(xsize)
initialforce= 1.9 #initial value for force (at river head)
force[0] = initialforce 		

for i in range(1, xsize):
	force[i] = force[i-1] + theta*(mu - force[i-1])*dx + sigma*dw[i] #Euler-M scheme: OU process

var = ((sigma**2) /(2*theta)) * (1-np.exp(-2*theta*x))	#Gardiner04 Eqn 4.4.29, with zero variance for x(0)
exp = mu + (force[0] - mu)*np.exp(-theta * x)
prob = 1 - norm.cdf((cmean-exp)/np.sqrt(var))	
probflip=np.flip(prob)
rundispa = np.cumsum(probflip) * dx					#cumulative expected displacement
runvar = (dx**2) * np.cumsum(probflip*(1-probflip))	#cumulative variance
 
xprob = np.array([x,prob])
xprob = xprob.T
np.savetxt(f4, xprob, delimiter=' ')

# fig=plt.figure(figsize=(6,8), facecolor='white')
# 
# ax = fig.add_subplot(211)
# plt.text(25, 1000, "Time step 0", size=14, horizontalalignment='left')
# plt.bar(x,z,color='gray',width=1)
# plt.xlabel('Distance, [L]')
# plt.ylabel('Elevation, [L]')
# plt.ylim(ymin=-50,ymax=1100)
# plt.plot(x, 400*prob, color='black', label='P(F>c)', linewidth=3)
# plt.text(1015, 400, "- p = 1", size=8, horizontalalignment='left')
# plt.text(1015, 0, "- p = 0", size=8, horizontalalignment='left')
# plt.xlim([-15, 1015])
# plt.errorbar(N - rundispa[0], 0, xerr=runvar[0], color='black')		#plot variance
# plt.scatter(N - rundispa[0], 0, color='red', edgecolors='black', zorder = 10)	#plot expected values]
# 
# axs_ins = inset_axes(ax, 
#                     width="25%", 	# width relative to width of parent_bbox
#                     height=1, 		# height in inches
#                     loc=1)
# plt.xlabel('$\Delta z$, [L]')
# plt.ylabel('PDF')
# plt.xlim(xmin=0,xmax=99)
# plt.ylim(ymin=0,ymax=0.1)
# dz = z[:-1] - z[1:] #actual dz
# binwidth=1
# plt.hist(dz, bins=np.arange(min(dz), max(dz) + binwidth, binwidth), color='gray', density=True)
# 
# ax = fig.add_subplot(212)
# plt.ylim([-0.2,5])
# plt.plot(x, cmean, color='orange', label='c, [L]', linewidth=2)
# plt.plot(x, exp, color='grey', linewidth=3, label='$<F_x>$')
# plt.plot(x, exp+np.sqrt(var), color='grey', linestyle='dashed')
# plt.plot(x, exp-np.sqrt(var), color='grey', linestyle='dashed')
# plt.plot(x, force, color='black', label='F, [M][L]/[T]$^2$', linewidth=1)
# plt.xlabel('Distance, [L]')
# plt.ylabel('Force, [L][M]/[T]$^2$')
# plt.legend()
# #plt.savefig('./riv_evo_OU/fig0.png', dpi=300)
# plt.show()
# plt.close()
# exit()

xforce = np.array([x,cmean, exp, exp+np.sqrt(var), exp-np.sqrt(var), force])
#xforce = np.array([x,cmean, force])
xforce = xforce.T
np.savetxt(f5, xforce, delimiter=' ')


#output starting solution for external plotting
time=np.full(N, times[0])
ssol = np.array([time, x, z])
ssol = ssol.T
np.savetxt('ssol.csv', ssol, delimiter=' ')

np.savetxt(f2, np.c_[time[0], N - rundispa[0], runvar[0]], delimiter=' ')

txexpd = np.array([times, rundispa, runvar])
txexpd = txexpd.T
np.savetxt(f7, txexpd, delimiter=' ')


runs = 1000			#number of time steps that the model will take

j=50
for i in range(runs):
	dw = np.random.normal(0, scale = np.sqrt(dx), size = xsize)
	force = np.zeros(xsize)
	force[0] = initialforce #initial value for force (at river head)		

	for l in range(1, xsize):
		force[l] = force[l-1] + theta*(mu - force[l-1])*dx + sigma*dw[l] #Euler-M scheme: OU process

	for k in range(N-1):
		if force[k] > cmean[k]:
			z[k] = z[k+1]
			
	if i >=590 and i <=610: 
		time=np.full(N, times[i])
		rivevo = np.array([time, x, z])
		rivevo = rivevo.T
		np.savetxt(f3, rivevo, delimiter=' ')

	if i == j :	
		dz = z[:-1] - z[1:] #actual dz
	
# 		fig=plt.figure(figsize=(6,8), facecolor='white')
# 
# 		ax = fig.add_subplot(211)
# 		plt.text(25, 1000, "Time step {}".format(i), size=14, horizontalalignment='left')
# 		plt.bar(x,z,color='gray',width=1)
# 		plt.xlabel('Distance, [L]')
# 		plt.ylabel('Elevation, [L]')
# 		plt.plot(x, 400*prob, color='black', label='P(F>c)', linewidth=3)
# 		plt.text(1015, 400, "- p = 1", size=8, horizontalalignment='left')
# 		plt.text(1015, 0, "- p = 0", size=8, horizontalalignment='left')
# 		plt.xlim([-15, 1015])
# 		plt.ylim(ymin=-50,ymax=1100)
# 		plt.errorbar(N - rundispa[i], 0, xerr=runvar[i], color='black')		#plot variance
# 		plt.scatter(N - rundispa[i], 0, color='red', edgecolors='black', zorder = 10)	#plot expected values]
# 
# 		
# 		axs_ins = inset_axes(ax, 
#                     	width="25%", 	# width relative to width of parent_bbox
#                     	height=1, 		# height in inches
#                     	loc=1)
# 		plt.xlabel('$\Delta z$, [L]')
# 		plt.ylabel('PDF')
# 		plt.xlim(xmin=0,xmax=99)
# 		plt.ylim(ymin=0,ymax=0.1)
# 		dz = z[:-1] - z[1:] #actual dz
# 		binwidth=1
# 		plt.hist(dz, bins=np.arange(min(dz), max(dz) + binwidth, binwidth), color='gray', density=True)
# 
# 		ax = fig.add_subplot(212)
# 		plt.ylim([-0.2,5])
# 		plt.plot(x, cmean, color='orange', label='c, [L]', linewidth=2)
# 		plt.plot(x, exp, color='grey', linewidth=3, label='$<F_x>$')
# 		plt.plot(x, exp+np.sqrt(var), color='grey', linestyle='dashed')
# 		plt.plot(x, exp-np.sqrt(var), color='grey', linestyle='dashed')
# 		plt.plot(x, force, color='black', label='F, [M][L]/[T]$^2$', linewidth=1)
# 		plt.xlabel('Distance, [L]')
# 		plt.ylabel('Force, [L][M]/[T]$^2$')
# #		plt.legend()
# 
# 
# #		plt.show()
# 		plt.savefig('./riv_evo_OU/fig%s.png' % i, dpi=300)
# 		plt.close()
# 		print('fig%s.png' % i)
		
		
		time=np.full(N, times[i])
		
		rivevo = np.array([time, x, z])
		rivevo = rivevo.T
		np.savetxt(f, rivevo, delimiter=' ')

		np.savetxt(f2, np.c_[time[i], N - rundispa[i], runvar[i]], delimiter=' ')
		
		txforce = np.array([time,x,force])
		txforce = txforce.T
		np.savetxt(f6, txforce, delimiter=' ')

		j+=50
		
exit()

			
			
	