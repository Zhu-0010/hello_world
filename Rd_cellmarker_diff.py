import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyrenn_bak as prn
import json
from matplotlib.pyplot import MultipleLocator
from scipy.optimize import leastsq
from operator import itemgetter, attrgetter
from numpy import *
import glob
import os
import re
from scipy import stats
#import seaborn as sns

df=pd.read_csv('./PhenoGraph_Acsinh_csvs/PhenoGraph1_Acsinh_Expr.csv')
labels_=list(df.columns.values)
labels=labels_[2:]
#print(labels)

cell_type=['Treg','Treg','DNT','CD8+T','classical mono','classical mono','Tconv','NA','B','classical mono',
'Treg','Treg','Tconv','B','B','Tconv','B','Treg','Treg','Treg',
'classical mono','classical mono','CD8+T','CD8+T','Treg','Treg','classical mono','Treg','B','CD8+T',
'NA','classical mono','classical mono','CD8+T','B','classical mono','NA','NA','B','Treg',
'classical mono','B','classical mono','Treg','Tconv','Treg','B','Tconv','CD8+T','classical mono',
'NA','B']

#1 Moderate Patients
read_path='PhenoGraph_Acsinh_csvs_HC_W'
csv_read=glob.glob(os.path.join(read_path,'*.json'))

dic_weight_W={}
for k,path in enumerate(csv_read):
	with open(path) as f_obj:
		data=json.load(f_obj)
		data_np=np.array(data)
		m=data_np.shape
		#print(m)
		weight_k=[]
		filename = os.path.splitext(path)[0]
		regex = re.compile(r'\d+')
		num = int(max(regex.findall(filename)))
		#print(num)
		for i in range(0, m[0]):
			for j in data_np[i,1:]:
				weight_k.append(j)
		weight_np=np.array(weight_k)
		weight_mean=np.mean(weight_np, axis=0)
		#print(weight_mean)
		dic_weight_W[num]=weight_mean.tolist()
#print(dic_weight_mean)
#weight_mean_sorted=sorted(dic_weight_mean.items())
#print('weight_mean_sorted:', weight_mean_sorted)
Rd_W=[]
for key in sorted(dic_weight_W):
	Rd_W.append(dic_weight_W[key])
#print(Rd_W)
Rd_W_np=np.array(Rd_W)
#print(Rd_W_np.shape)

#2 Severe Patients
read_path='PhenoGraph_Acsinh_csvs_HC_ICU'
csv_read=glob.glob(os.path.join(read_path,'*.json'))

dic_weight_ICU={}
for k,path in enumerate(csv_read):
	with open(path) as f_obj:
		data=json.load(f_obj)
		data_np=np.array(data)
		m=data_np.shape
		#print(m)
		weight_k=[]
		filename = os.path.splitext(path)[0]
		regex = re.compile(r'\d+')
		num = int(max(regex.findall(filename)))
		#print(num)
		for i in range(0, m[0]):
			for j in data_np[i,1:]:
				weight_k.append(j)
		weight_np=np.array(weight_k)
		weight_mean=np.mean(weight_np, axis=0)
		#print(weight_mean)
		dic_weight_ICU[num]=weight_mean.tolist()
#print(dic_weight_mean)
#weight_mean_sorted=sorted(dic_weight_mean.items())
#print('weight_mean_sorted:', weight_mean_sorted)
Rd_ICU=[]
for key in sorted(dic_weight_ICU):
	Rd_ICU.append(dic_weight_ICU[key])
#print(Rd_ICU)
Rd_ICU_np=np.array(Rd_ICU)
#print(Rd_ICU_np.shape)

B_W=Rd_W_np[(8,13,14,16,28,34,38,41,46,51), :]
CD8_W=Rd_W_np[(3,22,23,29,33,48), :]
Treg_W=Rd_W_np[(0,1,10,11,17,18,19,24,25,27,39,43,45), :]
Tconv_W=Rd_W_np[(6,12,15,44,47), :]
DNT_W=Rd_W_np[2, :]
classical_mono_W=Rd_W_np[(4,5,9,20,21,26,31,32,35,40,42,49), :]

B_ICU=Rd_ICU_np[(8,13,14,16,28,34,38,41,46,51), :]
CD8_ICU=Rd_ICU_np[(3,22,23,29,33,48), :]
Treg_ICU=Rd_ICU_np[(0,1,10,11,17,18,19,24,25,27,39,43,45), :]
Tconv_ICU=Rd_ICU_np[(6,12,15,44,47), :]
DNT_ICU=Rd_ICU_np[2, :]
classical_mono_ICU=Rd_ICU_np[(4,5,9,20,21,26,31,32,35,40,42,49), :]
#print(B_W)

font1={'family': 'Times New Roman',
'weight': 'bold',
'size': 6}
font2={'family': 'Times New Roman',
'weight': 'bold',
'size': 8}
font3={'family': 'Times New Roman',
'weight': 'bold',
'size': 7}

width=0.4
legend_labels=['Moderate','Severe']
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
	plt.subplots_adjust(hspace=0.7, wspace=0.3, bottom=0.12, top=0.9)

def convert_pvalue(pvalue):
    if pvalue < 0.05:
        return "p<0.05"
    elif pvalue >= 0.05:
        return 'p=%.2f'%(pvalue)

def jitter_dots(dots):
    offsets = dots.get_offsets()
    jittered_offsets = offsets
    # only jitter in the x-direction
    jittered_offsets[:, 0] += np.random.uniform(-0.05, 0.05, offsets.shape[0])
    dots.set_offsets(jittered_offsets)


'''一、 CD8+T cells'''
plt.figure(figsize=(3.5, 6), dpi=300)
plt.subplot(2,2,1)
#plt.suptitle('A                                        Early Stage', x=0.35, fontsize=7, fontweight='bold')
x0=list(range(Rd_W_np.shape[1]))
CD8_W_mean=np.mean(CD8_W,axis=0)
CD8_W_std=np.std(CD8_W,axis=0)
CD8_ICU_mean=np.mean(CD8_ICU,axis=0)
CD8_ICU_std=np.std(CD8_ICU, axis=0)
CD8_W_max=np.max(CD8_W, axis=0)
CD8_ICU_max=np.max(CD8_ICU, axis=0)

Pvalue=[]
for i in range(0, Rd_W_np.shape[1]):
	CD8_W_i=CD8_W[:,i]
	CD8_ICU_i=CD8_ICU[:,i]
	stat, pvalue=stats.ttest_ind(CD8_W_i, CD8_ICU_i, equal_var = False)
	Pvalue.append(pvalue)

index_diff_CD8=[3,4,6]
x_diff_CD8=list(range(len(index_diff_CD8)))
print('x_diff_CD8:', x_diff_CD8)
CD8_W_diff=CD8_W[:,(3,4,6)]
CD8_ICU_diff=CD8_ICU[:,(3,4,6)]
CD8_W_mean_diff=[CD8_W_mean[i] for i in index_diff_CD8]
CD8_W_std_diff=[CD8_W_std[i] for i in index_diff_CD8]
CD8_ICU_mean_diff=[CD8_ICU_mean[i] for i in index_diff_CD8]
CD8_ICU_std_diff=[CD8_ICU_std[i] for i in index_diff_CD8]
Pvalue_CD8_diff=[Pvalue[i] for i in index_diff_CD8]
CD8_labels_diff=[labels[i] for i in index_diff_CD8]
CD8_W_max_diff=[CD8_W_max[i] for i in index_diff_CD8]
CD8_ICU_max_diff=[CD8_ICU_max[i] for i in index_diff_CD8]

# moderate
for i,y in enumerate(CD8_W_diff):
	scatter_W=plt.scatter(x_diff_CD8, y, color='', marker='o', edgecolors='b', linewidths=0.8, s=12)
	jitter_dots(scatter_W)
plt.errorbar(x=x_diff_CD8, y=CD8_W_mean_diff, yerr=CD8_W_std_diff, fmt='k_', elinewidth=1, capsize=2, capthick=1 )

# ICU
width=0.3
for i in range(len(index_diff_CD8)):
	x_diff_CD8[i]=x_diff_CD8[i]+width
	CD8_max_i=np.max([CD8_W_max_diff[i], CD8_ICU_max_diff[i]])
	plt.plot([x_diff_CD8[i]-width, x_diff_CD8[i]], [CD8_max_i+8., CD8_max_i+8.], lw=0.5, c="k")
	plt.text((x_diff_CD8[i]-width+x_diff_CD8[i])*.5, CD8_max_i+8., convert_pvalue(Pvalue_CD8_diff[i]), ha='center', va='bottom', color="k", fontsize=4.5)
for j,y in enumerate(CD8_ICU_diff):
	scatter_ICU=plt.scatter(x_diff_CD8, y, color='', marker='s', edgecolors='r', linewidths=0.8, s=12)
	jitter_dots(scatter_ICU)
plt.errorbar(x=x_diff_CD8, y=CD8_ICU_mean_diff, yerr=CD8_ICU_std_diff, fmt='k_',  elinewidth=1, capsize=2, capthick=1 )
plt.title('CD8+ T', fontdict=font3)

xticks=[x_diff_CD8[i]-width/2 for i in range(len(index_diff_CD8))]
plt.xticks(xticks, labels=CD8_labels_diff, rotation=30, fontsize=5.5, weight='bold')
plt.xlim(-0.5, 2.5)
plt.ylim(-25, 90)
plt.yticks(np.arange(-25,76,25), fontsize=5, weight='bold')
plt.legend(handles=[scatter_W, scatter_ICU], labels=legend_labels, loc='best', fontsize=4)
plt_config()
plt.ylabel('RS', fontdict=font2 )


'''二、Treg'''
plt.subplot(2,2,2)
Treg_W_mean=np.mean(Treg_W,axis=0)
Treg_W_std=np.std(Treg_W,axis=0)
Treg_ICU_mean=np.mean(Treg_ICU,axis=0)
Treg_ICU_std=np.std(Treg_ICU,axis=0)
Treg_W_max=np.max(Treg_W,axis=0)
Treg_ICU_max=np.max(Treg_ICU,axis=0)


Pvalue=[]
for i in range(0, Rd_W_np.shape[1]):
	Treg_W_i=Treg_W[:,i]
	Treg_ICU_i=Treg_ICU[:,i]
	stat, pvalue=stats.ttest_ind(Treg_W_i, Treg_ICU_i, equal_var = False)
	Pvalue.append(pvalue)

index_diff_Treg=[0,3,4]
x_diff_Treg=list(range(len(index_diff_Treg)))
#print('x_diff_CD8:', x_diff_CD8)
Treg_W_diff=Treg_W[:,(0,3,4,)]
Treg_ICU_diff=Treg_ICU[:,(0,3,4)]
Treg_W_mean_diff=[Treg_W_mean[i] for i in index_diff_Treg]
Treg_W_std_diff=[Treg_W_std[i] for i in index_diff_Treg]
Treg_ICU_mean_diff=[Treg_ICU_mean[i] for i in index_diff_Treg]
Treg_ICU_std_diff=[Treg_ICU_std[i] for i in index_diff_Treg]
Pvalue_Treg_diff=[Pvalue[i] for i in index_diff_Treg]
Treg_labels_diff=[labels[i] for i in index_diff_Treg]
Treg_W_max_diff=[Treg_W_max[i] for i in index_diff_Treg]
Treg_ICU_max_diff=[Treg_ICU_max[i] for i in index_diff_Treg]

# moderate
for i,y in enumerate(Treg_W_diff):
	scatter_W=plt.scatter(x_diff_Treg, y, color='', marker='o', edgecolors='b', linewidths=0.8, s=12)
	jitter_dots(scatter_W)
plt.errorbar(x=x_diff_Treg, y=Treg_W_mean_diff, yerr=Treg_W_std_diff, fmt='k_', elinewidth=1, capsize=2, capthick=1 )

# ICU
width=0.3
for i in range(len(index_diff_Treg)):
	x_diff_Treg[i]=x_diff_Treg[i]+width
	Treg_max_i=np.max([Treg_W_max_diff[i], Treg_ICU_max_diff[i]])
	plt.plot([x_diff_Treg[i]-width, x_diff_Treg[i]], [Treg_max_i+8., Treg_max_i+8.], lw=0.5, c="k")
	plt.text((x_diff_Treg[i]-width+x_diff_Treg[i])*.5, Treg_max_i+8., convert_pvalue(Pvalue_Treg_diff[i]), ha='center', va='bottom', color="k", fontsize=4.5)
for j,y in enumerate(Treg_ICU_diff):
	scatter_ICU=plt.scatter(x_diff_Treg, y, color='', marker='s', edgecolors='r', linewidths=0.8, s=12)
	jitter_dots(scatter_ICU)
plt.errorbar(x=x_diff_Treg, y=Treg_ICU_mean_diff, yerr=Treg_ICU_std_diff, fmt='k_',  elinewidth=1, capsize=2, capthick=1 )
plt.title('Treg', fontdict=font3)

xticks=[x_diff_Treg[i]-width/2 for i in range(len(index_diff_Treg))]
plt.xticks(xticks, labels=Treg_labels_diff, rotation=30, fontsize=5.5, weight='bold')
plt.xlim(-0.5, 2.5)
plt.ylim(-50, 50)
plt.yticks(np.arange(-50,51,25), fontsize=5, weight='bold')
plt_config()
#plt.tight_layout()

'''三、B cells'''
plt.subplot(2,2,3)
B_W_mean=np.mean(B_W,axis=0)
B_W_std=np.std(B_W,axis=0)
B_ICU_mean=np.mean(B_ICU,axis=0)
B_ICU_std=np.std(B_ICU,axis=0)
B_W_max=np.max(B_W,axis=0)
B_ICU_max=np.max(B_ICU,axis=0)


Pvalue=[]
for i in range(0, Rd_W_np.shape[1]):
	B_W_i=B_W[:,i]
	B_ICU_i=B_ICU[:,i]
	stat, pvalue=stats.ttest_ind(B_W_i, B_ICU_i, equal_var = False)
	Pvalue.append(pvalue)

index_diff_B=[3,10,11,12]
x_diff_B=list(range(len(index_diff_B)))
#x_diff_Tconv=[0.4, 0.6]
print('x_diff_B:', x_diff_B)
B_W_diff=B_W[:,(3,10,11,12)]
print('B_W_diff:', B_W_diff)
B_ICU_diff=B_ICU[:,(3,10,11,12)]
B_W_mean_diff=[B_W_mean[i] for i in index_diff_B]
B_W_std_diff=[B_W_std[i] for i in index_diff_B]
B_ICU_mean_diff=[B_ICU_mean[i] for i in index_diff_B]
B_ICU_std_diff=[B_ICU_std[i] for i in index_diff_B]
Pvalue_B_diff=[Pvalue[i] for i in index_diff_B]
B_labels_diff=[labels[i] for i in index_diff_B]
B_W_max_diff=[B_W_max[i] for i in index_diff_B]
B_ICU_max_diff=[B_ICU_max[i] for i in index_diff_B]

# moderate
for i,y in enumerate(B_W_diff):
	scatter_W=plt.scatter(x_diff_B, y, color='', marker='o', edgecolors='b', linewidths=0.8, s=12)
	jitter_dots(scatter_W)
plt.errorbar(x=x_diff_B, y=B_W_mean_diff, yerr=B_W_std_diff, fmt='k_', elinewidth=1, capsize=2, capthick=1 )

# ICU
width=0.35
for i in range(len(index_diff_B)):
	x_diff_B[i]=x_diff_B[i]+width
	B_max_i=np.max([B_W_max_diff[i], B_ICU_max_diff[i]])
	plt.plot([x_diff_B[i]-width, x_diff_B[i]], [B_max_i+8., B_max_i+8.], lw=0.5, c="k")
	plt.text((x_diff_B[i]-width+x_diff_B[i])*.5, B_max_i+8., convert_pvalue(Pvalue_B_diff[i]), ha='center', va='bottom', color="k", fontsize=4.5)
for j,y in enumerate(B_ICU_diff):
	scatter_ICU=plt.scatter(x_diff_B, y, color='', marker='s', edgecolors='r', linewidths=0.8, s=12)
	jitter_dots(scatter_ICU)
plt.errorbar(x=x_diff_B, y=B_ICU_mean_diff, yerr=B_ICU_std_diff, fmt='k_',  elinewidth=1, capsize=2, capthick=1 )
plt.title('B cells', fontdict=font3)

xticks=[x_diff_B[i]-width/2 for i in range(len(index_diff_B))]
plt.xticks(xticks, labels=B_labels_diff, rotation=30, fontsize=5.5, weight='bold')
plt.xlim(-0.5, 3.5)
plt.ylim(-50, 25)
plt.yticks(np.arange(-50,26,25), fontsize=5, weight='bold')
plt_config()
plt.ylabel('RS', fontdict=font2 )

'''四、Classical monocytes'''
plt.subplot(2,2,4)
classical_mono_W_mean=np.mean(classical_mono_W,axis=0)
classical_mono_W_std=np.std(classical_mono_W,axis=0)
classical_mono_ICU_mean=np.mean(classical_mono_ICU,axis=0)
classical_mono_ICU_std=np.std(classical_mono_ICU,axis=0)
classical_mono_W_max=np.max(classical_mono_W,axis=0)
classical_mono_ICU_max=np.max(classical_mono_ICU,axis=0)


Pvalue=[]
for i in range(0, Rd_W_np.shape[1]):
	classical_mono_W_i=classical_mono_W[:,i]
	classical_mono_ICU_i=classical_mono_ICU[:,i]
	stat, pvalue=stats.ttest_ind(classical_mono_W_i, classical_mono_ICU_i, equal_var = False)
	Pvalue.append(pvalue)

index_diff_classical_mono=[0,2,3,9,11,12]
x_diff_classical_mono=list(range(len(index_diff_classical_mono)))
print('x_diff_Tconv:', x_diff_classical_mono)
classical_mono_W_diff=classical_mono_W[:,(0,2,3,9,11,12)]
print('classical_mono_W_diff:', classical_mono_W_diff)
classical_mono_ICU_diff=classical_mono_ICU[:,(0,2,3,9,11,12)]
classical_mono_W_mean_diff=[classical_mono_W_mean[i] for i in index_diff_classical_mono]
classical_mono_W_std_diff=[classical_mono_W_std[i] for i in index_diff_classical_mono]
classical_mono_ICU_mean_diff=[classical_mono_ICU_mean[i] for i in index_diff_classical_mono]
classical_mono_ICU_std_diff=[classical_mono_ICU_std[i] for i in index_diff_classical_mono]
Pvalue_classical_mono_diff=[Pvalue[i] for i in index_diff_classical_mono]
classical_mono_labels_diff=[labels[i] for i in index_diff_classical_mono]
classical_mono_W_max_diff=[classical_mono_W_max[i] for i in index_diff_classical_mono]
classical_mono_ICU_max_diff=[classical_mono_ICU_max[i] for i in index_diff_classical_mono]

# moderate
for i,y in enumerate(classical_mono_W_diff):
	scatter_W=plt.scatter(x_diff_classical_mono, y, color='', marker='o', edgecolors='b', linewidths=0.8, s=12)
	jitter_dots(scatter_W)
plt.errorbar(x=x_diff_classical_mono, y=classical_mono_W_mean_diff, yerr=classical_mono_W_std_diff, fmt='k_', elinewidth=1, capsize=2, capthick=1 )

# ICU
width=0.45
for i in range(len(index_diff_classical_mono)):
	x_diff_classical_mono[i]=x_diff_classical_mono[i]+width
	classical_mono_max_i=np.max([classical_mono_W_max_diff[i], classical_mono_ICU_max_diff[i]])
	plt.plot([x_diff_classical_mono[i]-width, x_diff_classical_mono[i]], [classical_mono_max_i+8., classical_mono_max_i+8.], lw=0.5, c="k")
	plt.text((x_diff_classical_mono[i]-width+x_diff_classical_mono[i])*.5, classical_mono_max_i+8., convert_pvalue(Pvalue_classical_mono_diff[i]), ha='center', va='bottom', color="k", fontsize=4)
for j,y in enumerate(classical_mono_ICU_diff):
	scatter_ICU=plt.scatter(x_diff_classical_mono, y, color='', marker='s', edgecolors='r', linewidths=0.8, s=12)
	jitter_dots(scatter_ICU)
plt.errorbar(x=x_diff_classical_mono, y=classical_mono_ICU_mean_diff, yerr=classical_mono_ICU_std_diff, fmt='k_',  elinewidth=1, capsize=2, capthick=1 )
plt.title('Classical monocytes', fontdict=font3)

xticks=[x_diff_classical_mono[i]-width/2 for i in range(len(index_diff_classical_mono))]
plt.xticks(xticks, labels=classical_mono_labels_diff, rotation=30, fontsize=5.5, weight='bold')
plt.xlim(-0.5, 5.7)
plt.ylim(-60, 60)
plt.yticks(np.arange(-50,51,25), fontsize=5, weight='bold')
plt_config()
#plt.tight_layout()

plt.show()
