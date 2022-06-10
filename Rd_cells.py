import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyrenn_bak as prn
import json
from matplotlib.pyplot import MultipleLocator
from scipy.optimize import leastsq
from operator import itemgetter, attrgetter
from numpy import *

F_T1='HC_W_T1.json'
with open(F_T1) as T1:
	data_T1=json.load(T1)
	#print('data:', data)
data_np_T1=np.array(data_T1)
Accuracy_Mean_T1=np.mean(data_np_T1[:,0], axis=0)
#print('Accuracy_Mean_T1: ', Accuracy_Mean_T1)
m_T1=data_np_T1.shape
#print(m)
weight_T1=[]
for j in range(0, m_T1[0]):
	for k in data_np_T1[j,1:]:
		weight_T1.append(k)
weight_np_T1=np.array(weight_T1)
weight_mean_T1=np.mean(weight_np_T1, axis=0)
#print("weight_mean_T1: ", weight_mean_T1)
#print("weight_mean_T1.shape: ", weight_mean_T1.shape)


F_T4='HC_W_T4.json'
with open(F_T4) as T4:
	data_T4=json.load(T4)
	#print('data:', data)
data_np_T4=np.array(data_T4)
Accuracy_Mean_T4=np.mean(data_np_T4[:,0], axis=0)
#print('Accuracy_Mean_T4: ', Accuracy_Mean_T4)
m_T4=data_np_T4.shape
#print(m)
weight_T4=[]
for j in range(0, m_T4[0]):
	for k in data_np_T4[j,1:]:
		weight_T4.append(k)
weight_np_T4=np.array(weight_T4)
weight_mean_T4=np.mean(weight_np_T4, axis=0)
#print("weight_mean_T4: ", weight_mean_T4)
#print("weight_mean_T4.shape: ", weight_mean_T4.shape)


F_T8='HC_W_T8.json'
with open(F_T8) as T8:
	data_T8=json.load(T8)
	#print('data:', data)
data_np_T8=np.array(data_T8)
Accuracy_Mean_T8=np.mean(data_np_T8[:,0], axis=0)
#print('Accuracy_Mean_T8: ', Accuracy_Mean_T8)
m_T8=data_np_T8.shape
#print(m)
weight_T8=[]
for j in range(0, m_T8[0]):
	for k in data_np_T8[j,1:]:
		weight_T8.append(k)
weight_np_T8=np.array(weight_T8)
weight_mean_T8=np.mean(weight_np_T8, axis=0)
#print("weight_mean_T8: ", weight_mean_T8)
#print("weight_mean_T8.shape: ", weight_mean_T8.shape)

weight_mean_T0=np.zeros((38,), dtype=float)
#print(weight_mean_T0)

#data_np=np.row_stack((weight_mean_T0, weight_mean_T1, weight_mean_T4, weight_mean_T8))
data_W_np=np.row_stack((weight_mean_T1, weight_mean_T4, weight_mean_T8))
#print('data_W_np: ', data_W_np)

F_T1='HC_ICU_T1.json'
with open(F_T1) as T1:
	data_T1=json.load(T1)
	#print('data:', data)
data_np_T1=np.array(data_T1)
Accuracy_Mean_T1=np.mean(data_np_T1[:,0], axis=0)
#print('Accuracy_Mean_T1: ', Accuracy_Mean_T1)
m_T1=data_np_T1.shape
#print(m)
weight_T1=[]
for j in range(0, m_T1[0]):
	for k in data_np_T1[j,1:]:
		weight_T1.append(k)
weight_np_T1=np.array(weight_T1)
weight_mean_T1=np.mean(weight_np_T1, axis=0)
#print("weight_mean_T1: ", weight_mean_T1)
#print("weight_mean_T1.shape: ", weight_mean_T1.shape)


F_T4='HC_ICU_T4.json'
with open(F_T4) as T4:
	data_T4=json.load(T4)
	#print('data:', data)
data_np_T4=np.array(data_T4)
Accuracy_Mean_T4=np.mean(data_np_T4[:,0], axis=0)
#print('Accuracy_Mean_T4: ', Accuracy_Mean_T4)
m_T4=data_np_T4.shape
#print(m)
weight_T4=[]
for j in range(0, m_T4[0]):
	for k in data_np_T4[j,1:]:
		weight_T4.append(k)
weight_np_T4=np.array(weight_T4)
weight_mean_T4=np.mean(weight_np_T4, axis=0)
#print("weight_mean_T4: ", weight_mean_T4)
#print("weight_mean_T4.shape: ", weight_mean_T4.shape)


F_T8='HC_ICU_T8.json'
with open(F_T8) as T8:
	data_T8=json.load(T8)
	#print('data:', data)
data_np_T8=np.array(data_T8)
Accuracy_Mean_T8=np.mean(data_np_T8[:,0], axis=0)
#print('Accuracy_Mean_T8: ', Accuracy_Mean_T8)
m_T8=data_np_T8.shape
#print(m)
weight_T8=[]
for j in range(0, m_T8[0]):
	for k in data_np_T8[j,1:]:
		weight_T8.append(k)
weight_np_T8=np.array(weight_T8)
weight_mean_T8=np.mean(weight_np_T8, axis=0)
#print("weight_mean_T8: ", weight_mean_T8)
#print("weight_mean_T8.shape: ", weight_mean_T8.shape)

weight_mean_T0=np.zeros((38,), dtype=float)
#print(weight_mean_T0)

#data_np=np.row_stack((weight_mean_T0, weight_mean_T1, weight_mean_T4, weight_mean_T8))
data_ICU_np=np.row_stack((weight_mean_T1, weight_mean_T4, weight_mean_T8))
#print('data_ICU_np: ', data_ICU_np)

T=[1,4,8]
cell_type=['neu-1','neu-2(CCR6+)','iEos','neu-3','neu-4','neu-5','plasma B','neu-6','neu-7',
'NA-1','NA-2','naive B','CD161+ effector \nmemory CD8 T','neu-8','neu-9','naive CD8 T',
'neu-10','neu-11','neu-12','neu-13','CD57hi memory CD4 T','Basophils','Classical mono','neu-14',
'neu-15','neu-16','pDC','rEos','neu-17','gamma delta T','naive CD4 T','CXCR3+ CCR6- \nmemory CD4 T',
'CD16lo NK','effector memory CD8 T','NA-3','Non-classical mono','CD57hi CD8 TEMRA','NA-4']
#n=cell_type[0]
#print('n:', n)
B_W=data_W_np[:, (6,11)]
CD8_W=data_W_np[:, (12, 15, 33, 36)]
CD4_W=data_W_np[:, (20,29,30,31)]
Innate_W=data_W_np[:,(2,21,22,26,27,32,35)]
Neu_W=data_W_np[:, (0,1,3,4,5,7,8,13,14,16,17,18,19,23,24,25,28)]
#Sig_pos_W=data_W_np[:,(22,32,1,7,15,24,25,16)]
#Sig_neg_W=data_W_np[:,(14,5,12,21,17,11,4,28,29)]
#Early_pos_W=Sig_pos_W[0,:]
#Early_neg_W=Sig_neg_W[0,:]
#print('Sig_pos_early_W:', Early_pos_W)
#print('Sig_pos_early_W:', Early_neg_W)

B_ICU=data_ICU_np[:, (6,11)]
CD8_ICU=data_ICU_np[:, (12, 15, 33, 36)]
CD4_ICU=data_ICU_np[:, (20,29,30,31)]
Innate_ICU=data_ICU_np[:,(2,21,22,26,27,32,35)]
Neu_ICU=data_ICU_np[:, (0,1,3,4,5,7,8,13,14,16,17,18,19,23,24,25,28)]
#Sig_pos_ICU=data_ICU_np[:,(22,32,1,7,15,24,25,16)]
#Sig_neg_ICU=data_ICU_np[:,(14,5,12,21,17,11,4,28,29)]
#Early_pos_ICU=Sig_pos_ICU[0,:]
#Early_neg_ICU=Sig_neg_ICU[0,:]

B_labels=[cell_type[i] for i in (6,11)]
CD8_labels=[cell_type[i] for i in (12, 15, 33, 36)]
CD4_labels=[cell_type[i] for i in (20,29,30,31)]
Innate_labels=[cell_type[i] for i in (2,21,22,26,27,32,35)]
Neu_labels=[cell_type[i] for i in (0,1,3,4,5,7,8,13,14,16,17,18,19,23,24,25,28)]
#Sig_pos_labels=[cell_type[i] for i in (22,32,1,7,15,24,25,16)]
#Sig_neg_labels=[cell_type[i] for i in (14,5,12,21,17,11,4,28,29)]
#print(B_labels)
"""
#1. B cells
fig, axes=plt.subplots(1,2,sharex=True,sharey=True)
scatters=[]
labels=['moderate','severe']
for i in range(0, B_W.shape[1]):
	B_W_i=B_W[:,i]
	B_ICU_i=B_ICU[:,i]
	plt.subplot(1,2,i+1)
	scatter_m=plt.scatter(T, B_W_i, marker='o')
	if i == 0:
		plt.legend(scatter_m, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	scatter_s=plt.scatter(T, B_ICU_i, marker='^')
	if i == 0:
		plt.legend(scatter_s, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	plt.title(B_labels[i])
	#plt.yticks(np.linspace(-0.4, 0.4, 5))
	#plt.ylabel('Relevance score', fontsize=10, rotation=90)
	plt.xticks(np.linspace(0,8,3))
	#plt.xlabel('T/Days', fontsize=10)
	plt.ylim(-1.6, 1.6)
	plt.axhline(y=0,ls=":",c="black")
	#plt.axhline(y=-0.1,ls=":",c="black")
	#plt.axhline(y=0.1,ls=":",c="black")
	#plt.tight_layout()
	plt.subplots_adjust(wspace=0.5, hspace=0.5)
	#if i==0:
	#	scatters=[scatter_w, scatter_s]
	#	plt.legend(scatters, labels=['moderate', 'severe'], bbox_to_anchor=(8.1,-1.5), loc='upper left', ncol=1)
fig.text(0.5, 0.05, 'T/Days', fontsize=15, ha='center')
fig.text(0.08, 0.5, 'Rd', fontsize=15, va='center', rotation=90)

#2. CD8+ T cells
fig, axes=plt.subplots(1,4,sharex=True,sharey=True)
scatters=[]
labels=['moderate','severe']
for i in range(0, CD8_W.shape[1]):
	CD8_W_i=CD8_W[:,i]
	CD8_ICU_i=CD8_ICU[:,i]
	plt.subplot(1,4,i+1)
	scatter_m=plt.scatter(T, CD8_W_i, marker='o')
	if i == 0:
		plt.legend(scatter_m, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	scatter_s=plt.scatter(T, CD8_ICU_i, marker='^')
	if i == 0:
		plt.legend(scatter_s, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	plt.title(CD8_labels[i])
	plt.xticks(np.linspace(0,8,3))
	plt.ylim(-1.6, 1.6)
	plt.axhline(y=0,ls=":",c="black")
	plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.text(0.5, 0.05, 'T/Days', fontsize=15, ha='center')
fig.text(0.08, 0.5, 'Rd', fontsize=15, va='center', rotation=90)

#3. CD4+ T cells
fig, axes=plt.subplots(1,2,sharex=True,sharey=True)
scatters=[]
labels=['moderate','severe']
for i in range(0, CD4_W.shape[1]):
	CD4_W_i=CD4_W[:,i]
	CD4_ICU_i=CD4_ICU[:,i]
	plt.subplot(1,4,i+1)
	scatter_m=plt.scatter(T, CD4_W_i, marker='o')
	if i == 0:
		plt.legend(scatter_m, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	scatter_s=plt.scatter(T, CD4_ICU_i, marker='^')
	if i == 0:
		plt.legend(scatter_s, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	plt.title(CD4_labels[i])
	plt.xticks(np.linspace(0,8,3))
	plt.ylim(-1.6, 1.6)
	plt.axhline(y=0,ls=":",c="black")
	plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.text(0.5, 0.05, 'T/Days', fontsize=15, ha='center')
fig.text(0.08, 0.5, 'Rd', fontsize=15, va='center', rotation=90)

#4. Innate Immune cells
fig, axes=plt.subplots(1,2,sharex=True,sharey=True)
scatters=[]
labels=['moderate','severe']
for i in range(0, Innate_W.shape[1]):
	Innate_W_i=Innate_W[:,i]
	Innate_ICU_i=Innate_ICU[:,i]
	plt.subplot(2,4,i+1)
	scatter_m=plt.scatter(T, Innate_W_i, marker='o')
	if i == 0:
		plt.legend(scatter_m, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	scatter_s=plt.scatter(T, Innate_ICU_i, marker='^')
	if i == 0:
		plt.legend(scatter_s, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	plt.title(Innate_labels[i])
	plt.xticks(np.linspace(0,8,3))
	plt.ylim(-1.6, 1.6)
	plt.axhline(y=0,ls=":",c="black")
	plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.text(0.5, 0.05, 'T/Days', fontsize=15, ha='center')
fig.text(0.08, 0.5, 'Rd', fontsize=15, va='center', rotation=90)

#5. Neutrophils
fig, axes=plt.subplots(1,2,sharex=True,sharey=True)
scatters=[]
labels=['moderate','severe']
for i in range(0, Neu_W.shape[1]):
	Neu_W_i=Neu_W[:,i]
	Neu_ICU_i=Neu_ICU[:,i]
	plt.subplot(3,6,i+1)
	scatter_m=plt.scatter(T, Neu_W_i, marker='o')
	if i == 0:
		plt.legend(scatter_m, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	scatter_s=plt.scatter(T, Neu_ICU_i, marker='^')
	if i == 0:
		plt.legend(scatter_s, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	plt.title(Neu_labels[i])
	plt.xticks(np.linspace(0,8,3))
	plt.ylim(-1.6, 1.6)
	plt.axhline(y=0,ls=":",c="black")
	plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.text(0.5, 0.05, 'T/Days', fontsize=15, ha='center')
fig.text(0.08, 0.5, 'Rd', fontsize=15, va='center', rotation=90)

#6. Sig_pos
fig, axes=plt.subplots(1,2,sharex=True,sharey=True)
scatters=[]
labels=['moderate','severe']
for i in range(0, Sig_pos_W.shape[1]):
	Sig_pos_W_i=Sig_pos_W[:,i]
	Sig_pos_ICU_i=Sig_pos_ICU[:,i]
	plt.subplot(2,4,i+1)
	scatter_m=plt.scatter(T, Sig_pos_W_i, marker='o')
	if i == 0:
		plt.legend(scatter_m, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	scatter_s=plt.scatter(T, Sig_pos_ICU_i, marker='^')
	if i == 0:
		plt.legend(scatter_s, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	plt.title(Sig_pos_labels[i])
	plt.xticks(np.linspace(0,8,3))
	plt.ylim(-1.6, 1.6)
	plt.axhline(y=0,ls=":",c="black")
	plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.text(0.5, 0.05, 'T/Days', fontsize=15, ha='center')
fig.text(0.08, 0.5, 'Rd', fontsize=15, va='center', rotation=90)

#6. Sig_neg
fig, axes=plt.subplots(1,2,sharex=True,sharey=True)
scatters=[]
labels=['moderate','severe']
for i in range(0, Sig_neg_W.shape[1]):
	Sig_neg_W_i=Sig_neg_W[:,i]
	Sig_neg_ICU_i=Sig_neg_ICU[:,i]
	plt.subplot(2,5,i+1)
	scatter_m=plt.scatter(T, Sig_neg_W_i, marker='o')
	if i == 0:
		plt.legend(scatter_m, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	scatter_s=plt.scatter(T, Sig_neg_ICU_i, marker='^')
	if i == 0:
		plt.legend(scatter_s, labels=labels, bbox_to_anchor=(11.8,-1.5), loc='upper left', ncol=1)
	plt.title(Sig_neg_labels[i])
	plt.xticks(np.linspace(0,8,3))
	plt.ylim(-1.6, 1.6)
	plt.axhline(y=0,ls=":",c="black")
	plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.text(0.5, 0.05, 'T/Days', fontsize=15, ha='center')
fig.text(0.08, 0.5, 'Rd', fontsize=15, va='center', rotation=90)

#7. Early stage of moderate patients
fig=plt.figure(figsize=(15,10))
ax0=fig.add_subplot(1,2,1)
ax1=fig.add_subplot(1,2,2)
total_width, n = 0.8, 2
width = total_width / n
x0=list(range(len(Early_pos_W)))
ax0.bar(x0, Early_pos_W, width=width, label='Moderate patients', fc='b')
for i in range(len(x0)):
	x0[i]=x0[i]+width
ax0.bar(x0, Early_pos_ICU, width=width, label='Severe patients', fc='r', tick_label=Sig_pos_labels)

#print('Sig_pos_labels:', Sig_pos_labels)
#print('len(Sig_pos_early_W:', len(Early_pos_W))
#x2=range(len(Early_pos_W), len(Early_neg_W)+len(Early_pos_W))
#plt.bar(x2, Early_neg_W, tick_label=Sig_neg_labels)
x1=list(range(len(Early_neg_W)))
ax1.bar(x1, Early_neg_W, width=width, label='Moderate patients', fc='b')
for i in range(len(x1)):
	x1[i]=x1[i]+width
ax1.bar(x1, Early_neg_ICU, width=width, label='Severe patients', fc='r', tick_label=Sig_neg_labels)
plt.legend()
"""
font1={'family': 'Times New Roman',
'weight': 'bold',
'size': 5.5}

font2={'family': 'Times New Roman',
'weight': 'bold',
'size': 7}

#y_major_locator=MultipleLocator(0.2)
#y_locator=MultipleLocator(0.5)

#8. Early stage of moderate patients
# sorted
key_value_T1={}
for i in range(0,data_W_np.shape[1]):
	key_value_T1[cell_type[i]]=[data_W_np[0,i], data_ICU_np[0,i]]
#print('key_value_T1:', key_value_T1)
dic_T1=sorted(key_value_T1.items(), key=itemgetter(1))
#print(dic_T1)
dic_np_T1=np.array(dic_T1)
#print(dic_np_T1.shape)
cell_type_T1=dic_np_T1[:,0]
#print('cell_type_T1:', cell_type_T1 )
#print(cell_type_T1.shape)
R_T1=[]
for j in range(0, dic_np_T1.shape[0]):
	for k in dic_np_T1[j,1:]:
		R_T1.append(k)
R_np_T1=np.array(R_T1)
#print('R_np_T1:', R_np_T1)
#print(R_np_T1.shape)

total_width, n = 0.8, 2
width = total_width / n
Early_pos_W=R_np_T1[-8:,0][::-1]
Early_pos_ICU=R_np_T1[-8:,1][::-1]
Sig_pos_labels=cell_type_T1[-8:][::-1]
Early_neg_W=R_np_T1[0:8,0]
Early_neg_ICU=R_np_T1[0:8,1]
Sig_neg_labels=cell_type_T1[0:8]

plt.figure(figsize=(3.5,2.5), dpi=300)
plt.suptitle('A                                        Early Stage', x=0.35, fontsize=7, fontweight='bold')
ax0=plt.subplot(1,2,1)
x0=list(range(len(Early_pos_W)))
ax0.bar(x0, Early_pos_W, width=width, label='Moderate', tick_label = Sig_pos_labels, fc='b')
for i in range(len(x0)):
	x0[i]=x0[i]+width/2
xticks=x0
ax0.set_xticks(xticks)
ax0.set_xticklabels(rotation=90, labels=Sig_pos_labels, fontdict=font1)
ax0.set_ylabel('RS', fontdict=font2 )
ax0.set_ylim(-0.2, 1.6)
plt.yticks(np.arange(-0.2, 1.8, 0.2), fontsize=5, weight='bold')
plt.axhline(y=0,ls="-",c="black", linewidth=1)
ax=plt.gca()
ax.spines['left'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
#ax.yaxis.set_major_locator(y_locator)
for i in range(len(x0)):
	x0[i]=x0[i]+width/2
ax0.bar(x0, Early_pos_ICU, width=width, label='Severe', fc='r')
ax0.legend(fontsize=5)
plt.tight_layout()

ax1=plt.subplot(1,2,2)
x1=list(range(len(Early_neg_W)))
ax1.bar(x1, Early_neg_W, width=width, label='Moderate', tick_label=Sig_neg_labels, fc='b')
for i in range(len(x1)):
	x1[i]=x1[i]+width/2
xticks=x1
ax1.set_xticks(xticks)
ax1.set_xticklabels(rotation=90, labels=Sig_neg_labels, fontdict=font1)
ax1.set_ylim(-1.6, 0.2)
#ax1.set_ylabel('Rd', fontdict=font2)
plt.axhline(y=0,ls="-",c="black", linewidth=1)
plt.yticks(np.arange(-1.6, 0.4, 0.2), fontsize=5, weight='bold')
ax=plt.gca()
ax.spines['left'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
#ax.yaxis.set_major_locator(y_locator)
for i in range(len(x1)):
	x1[i]=x1[i]+width/2
ax1.bar(x1, Early_neg_ICU, width=width, label='Severe', fc='r')
plt.tight_layout()
plt.subplots_adjust(top=0.9)

#9. Middle stage of moderate patients
key_value_T4={}
for i in range(0, data_W_np.shape[1]):
	key_value_T4[cell_type[i]]=[data_W_np[1,i], data_ICU_np[1,i]]
dic_T4=sorted(key_value_T4.items(), key=itemgetter(1))
dic_np_T4=np.array(dic_T4)
cell_type_T4=dic_np_T4[:,0]
R_T4=[]
for j in range(0, dic_np_T4.shape[0]):
	for k in dic_np_T4[j,1:]:
		R_T4.append(k)
R_np_T4=np.array(R_T4)

total_width, n = 0.8, 2
width = total_width / n
Middle_pos_W=R_np_T4[-8:,0][::-1]
Middle_pos_ICU=R_np_T4[-8:,1][::-1]
Middle_pos_labels=cell_type_T4[-8:][::-1]
Middle_neg_W=R_np_T4[0:8,0]
Middle_neg_ICU=R_np_T4[0:8,1]
Middle_neg_labels=cell_type_T4[0:8]

fig=plt.figure(figsize=(3.5,2.5),dpi=300)
plt.suptitle('B                                        Middle Stage', x=0.35, fontsize=7,fontweight='bold')
ax0=fig.add_subplot(1,2,1)
x0=list(range(len(Middle_pos_W)))
plt.bar(x0, Middle_pos_W, width=width, label='Moderate', tick_label = Middle_pos_labels, fc='b')
for i in range(len(x0)):
	x0[i]=x0[i]+width/2
xticks=x0
ax0.set_xticks(xticks)
ax0.set_xticklabels(rotation=90, labels=Middle_pos_labels, fontdict=font1)
ax0.set_ylim(-0.2,0.6)
plt.ylabel('RS', fontdict=font2 )
plt.yticks(np.arange(-0.2,0.8,0.2),fontsize=5, weight='bold')
plt.axhline(y=0,ls="-",c="black", linewidth=1)
ax=plt.gca()
ax.spines['left'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
#ax.yaxis.set_major_locator(y_major_locator)

for i in range(len(x0)):
	x0[i]=x0[i]+width/2
plt.bar(x0, Middle_pos_ICU, width=width, label='Severe patients', fc='r')
#plt.legend(prop=font1)
plt.tight_layout()

ax1=fig.add_subplot(1,2,2)
x1=list(range(len(Middle_neg_W)))
ax1.bar(x1, Middle_neg_W, width=width, label='Moderate', tick_label=Middle_neg_labels, fc='b')
for i in range(len(x1)):
	x1[i]=x1[i]+width/2
xticks=x1
ax1.set_xticks(xticks)
ax1.set_xticklabels(rotation=90, labels=Middle_neg_labels, fontdict=font1)
ax1.set_ylim(-0.6,0.2)
#ax1.set_ylabel('Rd', fontdict=font2)
plt.axhline(y=0,ls="-",c="black", linewidth=1)
plt.yticks(np.arange(-0.6,0.4,0.2),fontsize=5, weight='bold')
ax=plt.gca()
ax.spines['left'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
#ax.yaxis.set_major_locator(y_major_locator)

for i in range(len(x1)):
	x1[i]=x1[i]+width/2
ax1.bar(x1, Middle_neg_ICU, width=width, label='Severe patients', fc='r')
plt.tight_layout()
plt.subplots_adjust(top=0.9)

#10. Late stage of moderate patients
key_value_T8={}
for i in range(0, data_W_np.shape[1]):
	key_value_T8[cell_type[i]]=[data_W_np[2,i], data_ICU_np[2,i]]
dic_T8=sorted(key_value_T8.items(), key=itemgetter(1))
dic_np_T8=np.array(dic_T8)
cell_type_T8=dic_np_T8[:,0]
R_T8=[]
for j in range(0, dic_np_T8.shape[0]):
	for k in dic_np_T8[j,1:]:
		R_T8.append(k)
R_np_T8=np.array(R_T8)

total_width, n = 0.8, 2
width = total_width / n
Late_pos_W=R_np_T8[-8:,0][::-1]
Late_pos_ICU=R_np_T8[-8:,1][::-1]
Late_pos_labels=cell_type_T8[-8:][::-1]
Late_neg_W=R_np_T8[0:8,0]
Late_neg_ICU=R_np_T8[0:8,1]
Late_neg_labels=cell_type_T8[0:8]

plt.figure(figsize=(3.5,2.5), dpi=300)
plt.suptitle('C                                        Late Stage', x=0.35, fontsize=7, fontweight='bold')
ax0=plt.subplot(1,2,1)
x0=list(range(len(Late_pos_W)))
plt.bar(x0, Late_pos_W, width=width, label='Moderate', tick_label = Late_pos_labels, fc='b')
for i in range(len(x0)):
	x0[i]=x0[i]+width/2
xticks=x0
ax0.set_xticks(xticks)
ax0.set_xticklabels(rotation=90, labels=Late_pos_labels, fontdict=font1)
ax0.set_ylim(-0.4,0.6)
plt.ylabel('RS', fontdict=font2 )
plt.yticks(np.arange(-0.4,0.7,0.2),fontsize=5, weight='bold')
plt.axhline(y=0,ls="-",c="black", linewidth=1)
ax=plt.gca()
ax.spines['left'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
#ax.yaxis.set_major_locator(y_major_locator)

for i in range(len(x0)):
	x0[i]=x0[i]+width/2
ax0.bar(x0, Late_pos_ICU, width=width, label='Severe', fc='r')
#plt.legend(prop=font1)
plt.tight_layout()

ax1=plt.subplot(1,2,2)
x1=list(range(len(Late_neg_W)))
ax1.bar(x1, Late_neg_W, width=width, label='Moderate', tick_label=Late_neg_labels, fc='b')
for i in range(len(x1)):
	x1[i]=x1[i]+width/2
xticks=x1
ax1.set_xticks(xticks)
ax1.set_xticklabels(rotation=90, labels=Late_neg_labels, fontsize=5, weight='bold')
ax1.set_ylim(-0.6,0.4)
#ax1.set_ylabel('Rd', fontdict=font2)
plt.axhline(y=0,ls="-",c="black", linewidth=1)
plt.yticks(np.arange(-0.6,0.6,0.2),fontsize=5, weight='bold')
ax=plt.gca()
ax.spines['left'].set_linewidth(1)
ax.spines['right'].set_linewidth(1)
ax.spines['bottom'].set_linewidth(1)
ax.spines['top'].set_linewidth(1)
#ax.yaxis.set_major_locator(y_major_locator)

for i in range(len(x1)):
	x1[i]=x1[i]+width/2
ax1.bar(x1, Late_neg_ICU, width=width, label='Severe', fc='r')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
#plt.subplots_adjust(wspace=0.15, left=0.15, bottom=0.4, top=0.9)

plt.show()
