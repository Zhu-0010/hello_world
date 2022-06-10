import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyrenn_bak as prn
import json
import glob 
import os
import csv
import re
"""
#test
df_i=pd.read_csv('./PhenoGraph_Acsinh_csvs_HC_ICU/PhenoGraph1_Acsinh_Expr.csv')
P_i=df_i.iloc[[0,1,2,3,4,5,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48], 1:15].values
labels=list(df_i.columns.values)
Pi=P_i.astype(np.float64)
n=int(np.shape(Pi)[1])
print('Pi: ', Pi)
print('Pi.shape:', Pi.shape)
print(labels)
"""
read_path='PhenoGraph_Acsinh_csvs_HC_ICU'
csv_read=glob.glob(os.path.join(read_path,'*.csv'))
#count=0
active=0
while active<20:
	for i,path in enumerate(csv_read):
		df_i=pd.read_csv(path)
		#df_i=pd.read_csv('PhenoGraph_Acsinh_csvs/PhenoGraph1_Acsinh_Expr.csv')
		P_i=df_i.iloc[[0,1,2,3,4,5,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48], 1:15].values
		labels=list(df_i.columns.values)
		Pi=P_i.astype(np.float64)
		n=int(np.shape(Pi)[1])
		#print('Pi: ', Pi)
		#print('Pi.shape:', Pi.shape)
		#P_ext=prn.Matrix_extend(P)
		np.random.shuffle(Pi)
		#print('P_ext: ', P_ext)
		#print('P_ext.shape:', P_ext.shape)

		Xi=Pi[0:18,1:]
		Xi=Xi.T
		#print('Xi.shape:', Xi.shape)
		#print(Xi)
		Yi=Pi[0:18,0]
		Yi=Yi.T
		#print('Yi.shape:', Yi.shape)
		#print(Yi)

		Xti=Pi[18:,1:]
		Xti=Xti.T
		#print('Xti.shape:', Xti.shape)
		#print(Xti)
		Yti=Pi[18:,0]
		Yti=Yti.T
		#print('Yti.shape:', Yti.shape)

		net=prn.CreateNN([n-1,2,2,1], dIn=[0], dIntern=[], dOut=[])

		net=prn.train_LM(Xi, Yi, net, verbose=True)

		yi=prn.NNOut(Xi, net)
		yti=prn.NNOut(Xti, net)

		Yti_delta=Yti-yti
		ei=np.reshape(Yti_delta, (1, np.size(Yti_delta)), order='F')[0]
		MSEti=np.dot(ei, ei.transpose())
		#print('MSEti=', MSEti)
		
		accuracy_cnt=0
		for k in range(len(Yti)):
			if Yti[k] == min(1, 2, 3):
				if yti[k]-Yti[k] < 0.5:
					accuracy_cnt += 1
			elif Yti[k] == max(1, 2, 3):
				if yti[k]-Yti[k] > -0.5:
					accuracy_cnt += 1
			else:
				if abs(yti[k]-Yti[k]) <= 0.5:
					accuracy_cnt += 1
		Accuracy=float(accuracy_cnt)/len(Yti)
		print('Accuracy: ', Accuracy)


		IWi,LWi,bi=prn.w2Wb(net)
		#print('IW:', IW)
		#print('LW:', LW)
		#print('b:', b)
		#print(type(IW))
		IWi_110=IWi[(1,1,0)]
		#print('IW[(1,1,0)]:', IW[(1,1,0)])
		LWi_210=LWi[(2,1,0)]
		#print('LW_210:', LW_210)
		LWi_320=LWi[(3,2,0)]
		#iw_110=json.dumps(IW_110.tolist())
		#lw_210=json.dumps(LW_210.tolist())

		X_weight=[]
		for j in range(0, n-1):
			xj_=np.dot(IWi_110[:, j].T, LWi_210.T)
			xj=np.dot(xj_.T, LWi_320.T)
			#xi=xi.tolist()
			xj=float(xj)
			X_weight.append(xj)
		print('X_weight:', X_weight)
		#print(type(X_weight))

		datai=[Accuracy, X_weight]
		data_npi=np.array(datai)

		filename = os.path.splitext(path)[0]
		newfilename = '%s.json' %filename
		#fni=re.findall(r'\\(.+?)\.', path)
		#filenamei=str(fni)+'.json'
		try:
			with open(newfilename) as f_obji:
				data_oldi=json.load(f_obji)
				#print('data_old:', data_old)
		except FileNotFoundError:
			with open(newfilename, 'w') as f_obji:
				json.dump(datai, f_obji)
		else:
			data_old_npi=np.array(data_oldi)
			data_new_npi=np.row_stack((data_old_npi, data_npi))
			data_newi=data_new_npi.tolist()
			#print('data_new:', data_new)
			with open(newfilename, 'w') as f_obji:
				json.dump(data_newi, f_obji)
	active += 1

Accuracy_Mean=np.mean(data_new_np[:,0], axis=0)
print('Accuracy_Mean: ', Accuracy_Mean)


