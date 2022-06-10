import glob
import os
import pandas as pd
import numpy as np
import csv

read_path="PhenoGraph_Acsinh_csvs"
csv_read=glob.glob(os.path.join(read_path,'*.csv'))

with open('2.0_cluster_Percentage with group_cp.csv') as f:
	reader=csv.reader(f)
	header_row=next(reader)
	Patient_Type=[row[2] for row in reader]
print(header_row)
print(Patient_Type)

"""
Patient_Type=pd.read_csv('2.0_cluster_Percentage with group_cp.csv').iloc[:,2]
Patient_Type_2=Patient_Type.T
print(Patient_Type_2.shape)
"""

for i,path in enumerate(csv_read):
	df_i=pd.read_csv(path, index_col=0)
	#df_i=pd.read_csv(path)
	df_T=df_i.T
	col_name=df_T.columns.tolist()
	print(col_name)
	col_name.insert(0,'Patient_Type')	
	df_T=df_T.reindex(columns=col_name)
	#print(df_T)
	df_T['Patient_Type']=Patient_Type
	df_T.fillna(0, inplace=True)
	print(path)
	print(df_T)
	df_T.to_csv('./Transform/'+path)
