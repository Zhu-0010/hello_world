import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyrenn_bak as prn
import json
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import MultipleLocator

df=pd.read_csv('acc_hidden_neuron.csv')

# 1 hidden layer accuracy test
P1_=df.iloc[0:10, 0:3].values
P=P1_.astype(np.float64)
print('P: ', P)

X=P[:,0]
print('X:', X)
Y_CNN=P[:,1]
Y_RNN=P[:,2]
print('Y_CNN:', Y_CNN)
print('Y_RNN:', Y_RNN)

def plt_config():
	#plt.xticks(x0, labels=labels, rotation=90, fontsize=5, weight='bold')
	#plt.ylabel('Rs', fontdict=font2 )
	#plt.ylim(-60, 80)
	#plt.yticks(np.arange(-60,81,20), fontsize=5, weight='bold')
	plt.axhline(y=0, ls="--",c="black", linewidth=0.5)
	ax=plt.gca()
	ax.spines['left'].set_linewidth(1)
	ax.spines['right'].set_linewidth(0)
	ax.spines['bottom'].set_linewidth(1)
	ax.spines['top'].set_linewidth(0)
	#plt.legend()
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.5, bottom=0.1, top=0.9)

font1={'family': 'Times New Roman',
'weight': 'bold',
'size': 5}
font2={'family': 'Times New Roman',
'weight': 'bold',
'size': 6}
font3={'family': 'Times New Roman',
'weight': 'bold',
'size': 7}

fig=plt.figure(figsize=(3.5,2), dpi=300)
ax0=fig.add_subplot(1,2,1)
#ax1=fig.add_subplot(1,2,2, projection='3d')

ax0.set_title('C          Test Accuracy \n           of NN and RNN', font3, x=0.35, )
ax0.plot(X, Y_CNN, color='r', label='NN')
ax0.plot(X, Y_RNN, color='b', label='RNN')
ax0.legend(loc='best', fontsize=5)
xticks=list(range(0,11,2))
plt.xticks(xticks, rotation=0, fontsize=5, weight='bold')
plt.ylim(0.45, 0.75)
plt.yticks(np.arange(0.45, 0.76, 0.05), fontsize=5, weight='bold')
ax0.set_xlabel('Hidden Neuron Number', font1)
ax0.set_ylabel('Accuracy', font2)
plt.tight_layout()
plt.subplots_adjust(hspace=0.5, bottom=0.1, top=0.8)

# 2 hidden layers accuracy test
ax1=fig.add_subplot(1,2,2)
P2_=df.iloc[12:18, 0:6].values
P2=P2_.astype(np.float64)
print('P2: ', P2)

X1_=P2[1:,0]
print('X1_:', X1_)
X2_=P2[0,1:]
Y=P2[1:,1:]
print('X2_', X2_)
print('Y:', Y)
X1, X2= np.meshgrid(X1_, X2_)
print('X1:', X1)
print('X2:', X2)
#ax1.scatter(X1,X2,Y)
#ax1.plot_trisurf(X1,X2,Y)
#ax1.plot_wireframe(X1, X2, Y)
#ax1.plot_surface(X1,X2,Y)
acc_2h=ax1.contourf(X1, X2, Y, 10)
cbar = plt.colorbar(acc_2h) 
cbar.ax.tick_params(labelsize=4)
plt.title('D        Test Accuracy of 2 \n         hidden layers of NN', font3, x=0.45)
plt.xticks(fontsize=5, weight='bold')
plt.yticks(fontsize=5, weight='bold')
ax1.set_xlabel('Neuron numbers in hidden_layer 2', font1)
ax1.set_ylabel('Neuron numbers in hidden_layer 1', font1)
#cbar.set_label('colorbar', fontdict=font1)
ax1.xaxis.set_major_locator(MultipleLocator(1))
ax1.yaxis.set_major_locator(MultipleLocator(1))
plt.tight_layout()
plt.subplots_adjust(hspace=0.1, bottom=0.2, top=0.8)

plt.show()

