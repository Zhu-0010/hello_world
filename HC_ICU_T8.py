import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyrenn_bak as prn
import json
from matplotlib.pyplot import MultipleLocator

#df=pd.ExcelFile('cell_type_percent.csv').parse('Sheet0')
df=pd.read_csv('HC_ICU_T8.csv')
labels=list(df.columns.values)
print(labels)
P_=df.iloc[0:26, 2:41].values
P=P_.astype(np.float64)
print('P: ', P)

#P_ext=prn.Matrix_extend(P)
P_ext=P
np.random.shuffle(P_ext)
n=int(np.shape(P_ext)[1])
#print('P_ext: ', P_ext)
print('P_ext.shape[1]:', n)

X=P_ext[0:20,1:]
X=X.T
#print('X.shape:', X.shape)
#print(X)
Y=P_ext[0:20,0]
Y=Y.T
#print('Y.shape:', Y.shape)
#print(Y)

Xt=P_ext[20:,1:]
Xt=Xt.T
#print('Xt.shape:', Xt.shape)
#print(Xt)
Yt=P_ext[20:,0]
Yt=Yt.T
#print('Yt.shape:', Yt.shape)

net=prn.CreateNN([n-1,4,3,1], dIn=[0], dIntern=[], dOut=[])

net=prn.train_LM(X, Y, net, verbose=True, k_max=10000, E_stop=1e-4)

y=prn.NNOut(X, net)
yt=prn.NNOut(Xt, net)

Yt_delta=Yt-yt
e=np.reshape(Yt_delta, (1, np.size(Yt_delta)), order='F')[0]
MSEt=np.dot(e, e.transpose())
print('MSEt=', MSEt)

P_0=P_ext[:,0]
accuracy_cnt=0
for i in range(len(Yt)):
	if Yt[i] == min(1,2,3):
		if yt[i]-Yt[i] < 0.5:
			accuracy_cnt += 1
	elif Yt[i] == max(1,2,3):
		if yt[i]-Yt[i] > -0.5:
			accuracy_cnt += 1
	else:
		if abs(yt[i]-Yt[i]) <= 0.5:
			accuracy_cnt += 1
Accuracy=float(accuracy_cnt)/len(Yt)
print('Accuracy: ', Accuracy)

IW,LW,b=prn.w2Wb(net)
#print('IW:', IW)
#print('LW:', LW)
#print('b:', b)
#print(type(IW))
IW_110=IW[(1,1,0)]
#print('IW[(1,1,0)]:', IW[(1,1,0)])
LW_210=LW[(2,1,0)]
#print('LW_210:', LW_210)
LW_320=LW[(3,2,0)]
#iw_110=json.dumps(IW_110.tolist())
#lw_210=json.dumps(LW_210.tolist())

X_weight=[]
for i in range(0, n-1):
	#for j in range(0, hl_num):
	xi_1=np.dot(IW_110[:, i].T, LW_210.T)
	xi_2=np.dot(xi_1.T, LW_320.T)
	#xi=xi.tolist()
	xi=float(xi_2)
	X_weight.append(xi)
print('X_weight:', X_weight)
#print(type(X_weight))

data=[Accuracy, X_weight]
data_np=np.array(data)

filename='HC_ICU_T8.json'
try:
	with open(filename) as f_obj:
		data_old=json.load(f_obj)
		#print('data_old:', data_old)
except FileNotFoundError:
	with open(filename, 'w') as f_obj:
		json.dump(data, f_obj)
else:
	data_old_np=np.array(data_old)
	data_new_np=np.row_stack((data_old_np, data_np))
	data_new=data_new_np.tolist()
	#print('data_new:', data_new)
	with open(filename, 'w') as f_obj:
		json.dump(data_new, f_obj)

Accuracy_Mean=np.mean(data_new_np[:,0], axis=0)
print('Accuracy_Mean: ', Accuracy_Mean)

m=data_new_np.shape
print(m)
weight=[]
for j in range(0, m[0]):
	for k in data_new_np[j,1:]:
		#print(k.shape)
		#x=range(1,29)
		weight.append(k)
		plt.title('FR-FCM-Z3WR_panel_all')
		plt.plot(labels[3:-3], k, 'o:')
		#plt.ylim(-0.3, 0.3)
		y_major_locator=MultipleLocator(0.25)
		ax=plt.gca()
		ax.yaxis.set_major_locator(y_major_locator)
		plt.grid(True)
		#plt.grid(linestyle='--')
		plt.xticks(rotation=60, fontsize=10)
		plt.ylabel("Weight", fontsize=12)
		#plt.figure().autofmt_xdate()
weight_np=np.array(weight)
weight_mean=np.mean(weight_np, axis=0)
plt.plot(labels[3:-3], weight_mean, '-rs', markersize=10)

fig=plt.figure(figsize=(15,10))
ax0=fig.add_subplot(2,1,1)
ax1=fig.add_subplot(2,1,2)

ax0.set_title('Train Data')
ax0.plot(y, color='b', label='NN Output')
ax0.plot(Y, color='r', linestyle=':')
ax0.legend(loc='upper left')
ax0.set_ylabel('Severity')

ax1.set_title('Test Data')
ax1.plot(yt, color='b', label='NN Output')
ax1.plot(Yt, color='r', linestyle=':')
ax1.legend(loc='upper left')
ax1.set_ylabel('Severity')

#fig.autofmt_xdate()
#plt.show()

