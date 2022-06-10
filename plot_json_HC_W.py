import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.pyplot import MultipleLocator

df=pd.read_csv('./PhenoGraph_Acsinh_csvs/PhenoGraph1_Acsinh_Expr.csv')
labels=list(df.columns.values)
print(labels[2:])

read_path='PhenoGraph_Acsinh_csvs_HC_W'
csv_read=glob.glob(os.path.join(read_path,'*.json'))
#cluster_ID=6
#y_major_locator=MultipleLocator(1.0)
ax=plt.gca()

#fig=plt.figure(figsize=())
#ax0=fig.add_subplot(3,1,1)
#ax1=fig.add_subplot(3,1,2)

for k,path in enumerate(csv_read):
	#df_k=pd.read_csv(path)
	with open(path) as f_obj:
		data=json.load(f_obj)
		data_np=np.array(data)
		m=data_np.shape
		print(m)
		#MSE=data_np[:,0]
		weight_k=[]
		filename = os.path.splitext(path)[0]
		fig_k=plt.figure(figsize=(15,10))
		ax_0=fig_k.add_subplot(1,1,1)
		#ax_1=fig_k.add_subplot(2,1,2)
		for i in range(0, m[0]):
			for j in data_np[i,1:]:
				weight_k.append(j)
				#print(k.shape)
				#x=range(1,29)
				#fig=plt.figure()
				#ax_i=fig.add_subplot(6, 6, k%6+1)
				ax_0.set_title(filename, fontsize=20)
				ax_0.plot(labels[2:], j, 'o:')
				y_major_locator=MultipleLocator(1.0)
				ax.yaxis.set_major_locator(y_major_locator)
				ax_0.grid(True)
				ax_0.set_xticklabels(labels[2:], rotation=45, fontsize=14)
				#ax_0.set_ylabel("Weight", fontsize=18)
				#ax_1.plot(MSE)
				#ax_1.set_title('MSE', fontsize=20)
				#ax_1.set_xlabel('Repeat_Times', fontsize=18)
				#ax_1.set_ylabel('MSE', fontsize=18)
			#fig=plt.figure()
			#ax_k=fig.add_subplot(6, 6, k%6+1)
		weight_np=np.array(weight_k)
		weight_mean=np.mean(weight_np, axis=0)
		#print(weight_np)
		#print(weight_mean)
		#filename = os.path.splitext(path)[0]
		#newfilename = '%s.json' %filename
		ax_0.plot(labels[2:], weight_mean, '-rs', markersize=10)
		
		Accuracy_Mean=np.mean(data_np[:,0], axis=0)
		#print('Accuracy_Mean: ', Accuracy_Mean)
		ax_0.set_ylabel("Weight/Accuracy=%s"%Accuracy_Mean, fontsize=18)
		
		ax_0.figure.savefig('%s.png' %filename)
		#ax_1.figure.savefig('%s_MSE.jpg' %filename)
#plt.show()
