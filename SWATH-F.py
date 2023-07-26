# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 18:47:26 2023

@author: wxb2
"""
import pandas as pd

#Function 1: Format conversion of CorrDec export files
#The msp format export from MSDial-CorrDec was converted to a csv file containing the m/z,rt,MS2 fragments, abundance for each peak.
def read_data(filename):
    f1=open(filename,"r")
    filename_output=filename.split(".txt")[0]+"output.csv"
    data=f1.readlines()
    f1.close()
    f2=open(filename_output,"a")
    print("PRECURSORMZ","RETENTIONTIME","MZ","height",sep=",",file=f2)
    
    flag=0
    mz=""
    height=""
    for i in range(len(data)):
        if(data[i][:14]=="RETENTIONTIME:"):
            rt=data[i][15:].split("\n")[0]
            continue
        if(data[i][:12]=="PRECURSORMZ:"):
            pc=data[i][13:].split("\n")[0]
            print(pc,rt,end=",",sep=",",file=f2)
            continue
        if(data[i][:10]=="Num Peaks:"):
            flag=1
            numpeak=int(data[i][11:].split("\n")[0])
            if(numpeak<1):
                flag=0
                print("",file=f2)
                continue
            continue
        if(data[i]=="\n" or i==len(data)-1):
            print(mz,end=",",file=f2)
            print(height,file=f2)
            flag=0
            mz=""
            height=""
        if(flag==1):
            mz+=data[i].split("\t")[0]+"/"
            height+=data[i].split("\t")[1]+"/"
            continue   
    f2.close()
    return filename_output
#Function 2: Determination of retention time difference of homologs
#Function 2.1: Determine whether the retention time (rt) of any three peaks (n) satisfies a stepwise relationship (as n increases, the increment decreases).
def class_rtjudge3(list_n,list_rt,rt_error=1):
    if(list_rt[1]<list_rt[0]-rt_error or list_rt[2]<list_rt[1]-rt_error):
        return 0
    #The value of the third point is predicted based on the linear relationship between the first two points.
    pre_rt3=((list_rt[1]-list_rt[0])/(list_n[1]-list_n[0]))*(list_n[2]-list_n[1])+list_rt[1]
    if(pre_rt3<list_rt[2]-rt_error):
        return 0
    return 1
#Function 2.2: According to function 2.1, return the index number of the peaks in the series that satisfy the retention time difference condition, note that min_n is set to 3 here.
def class_rtjudge(line_index,line_n,line_rt,rt_error=1):
    #The final output is the index number,rt and n that satisfy the condition.
    line_index_match=[]
    line_n_match=[]
    line_rt_match=[]
    for i in range(len(line_index)):
        if(line_n[-1]-line_n[i]<2):
            continue
        for j in range(i+1,len(line_index)):
            if(line_n[-1]-line_n[j]<1):
                break
            if(line_n[j]==line_n[i]):
                continue
            for k in range(j+1,len(line_index)):
                if(line_n[k]==line_n[j]):
                    continue
                if(class_rtjudge3([line_n[i],line_n[j],line_n[k]],[line_rt[i],line_rt[j],line_rt[k]],rt_error=rt_error)):
                    if(line_index[i] not in line_index_match):
                        line_index_match.append(line_index[i])
                        line_n_match.append(line_n[i])
                        line_rt_match.append(line_rt[i])
                    if(line_index[j] not in line_index_match):
                        line_index_match.append(line_index[j])
                        line_n_match.append(line_n[j])
                        line_rt_match.append(line_rt[j])
                    if(line_index[k] not in line_index_match):
                        line_index_match.append(line_index[k])
                        line_n_match.append(line_n[k])
                        line_rt_match.append(line_rt[k])
    return line_index_match,line_n_match

# Function 3: Screening PFAS homolog series
# Setting parameters:
# homolog series difference: homo;
# relative error: mass_error (ppm);
# minimum number of substances required for a homologue series (multiple peaks of the same mass number are not included here): min_n;
# maximum number of series: homo_maxnum;
# The n value for each peak in each series: line_n.
#def class_find(filename,homo=49.99681,mass_error=5,rt_error=1,min_n=3):
def class_find(filename, homo=65.99171, mass_error=5, rt_error=1, min_n=3):
    # read file
    data=pd.read_csv(filename)
    # mass defect，<0.15 or >0.85
    index=[]
    for i in range(len(data)):
        # mz_md=float(data["PRECURSORMZ"][i])/49.99681*50
        mz_md = float(data["PRECURSORMZ"][i]) / 65.99171 * 66
        md=mz_md-int(mz_md)
        if(md>0.85 or md<0.15):
            index.append(i)
    data_md=data.iloc[index].reset_index(drop=True)
    data_md=data_md.sort_values(by="PRECURSORMZ").reset_index(drop=True)
    homo_maxnum=int((data_md["PRECURSORMZ"][len(data_md)-1]-data_md["PRECURSORMZ"][0])/homo)
    # class number
    class_index=[]
    class_num=0
    # The index number of the peak already in the series.
    peaks_used=[]
    line_n=[]
    for i in range(len(data_md)):
        if(i in peaks_used):
            continue
        # lines1 records the index number contained in the current series;
        # homo_n records the index number corresponding to the peak in the series.
        lines1=[]
        homo_n=[]
        # Identify the series corresponding to the i data in data_md.
        for j in range(homo_maxnum):
            mz_n=data_md["PRECURSORMZ"][i]+j*homo
            if(mz_n>data_md["PRECURSORMZ"][len(data_md)-1]+1):
                break
            data_md1=data_md[(data_md["PRECURSORMZ"]>mz_n*(1-mass_error*0.000001))&(data_md["PRECURSORMZ"]<mz_n*(1+mass_error*0.000001))]
            if(len(data_md1)==0):
                continue
            index1=list(data_md1.index)
            index1=[m for m in index1 if(m>=i)]
            lines1+=index1
            homo_n+=len(index1)*[j]
        # Determines whether the number of peaks in the series is greater than min_n.
        if(len(set(homo_n))<min_n):
            continue
        # Screening retention time
        line_rt=list(data_md.iloc[lines1]["RETENTIONTIME"])
        lines1,homo_n=class_rtjudge(lines1,homo_n,line_rt,rt_error=rt_error)
        if(len(set(homo_n))<min_n):
            continue
        
        class_num+=1
        class_index+=len(lines1)*[class_num]
        peaks_used+=lines1
        line_n+=homo_n    
    data_output=data_md.iloc[peaks_used].reset_index(drop=True)     
    data_output["Class"]=class_index 
    data_output["n"]=line_n
    data_output=data_output.reindex(columns=["Class","n","PRECURSORMZ","RETENTIONTIME","MZ","height"],fill_value=1)
    data_output.to_csv(filename,index=False)
# Function 4: Fragments.
def fragment_mark(filename_data,filename_database,error_ms2=20):
    data=pd.read_csv(filename_data)
    database=pd.read_csv(filename_database,sep="\t")
    database=database.sort_values(by="m/z").reset_index(drop=True)
    frag=[]
    for i in range(len(data)):   
        frag1=""
        mz=data["MZ"][i].split("/")[:-1]
        if(len(mz)<1):
            frag.append("")
        mz=[eval(i) for i in mz]
        for j in mz:
            for m in range(len(database)):
                if(j<database["m/z"][m]+1):
                    if(abs(j-database["m/z"][m])/j<error_ms2*0.000001):
                        frag1+=database["ion"][m]+"/"
                        break
        frag.append(frag1)
    data["PFAS_frag"]=frag
    data.to_csv(filename_data,index=False)

# Data file name
filename_data="tdcorrdec.txt"
# Fragmentation database name
filename_database="fragment database.txt"
# MS1 mass error during series screening.
error_ms1=5
# MS2 fragment mass error in fragment matching.
error_ms2=20
# Retention time errors for homolog series.
error_rt=2
# Homolog series unit mass.
# homo=49.99681(CF2)/65.99171(CF2O)
homo=65.99171
# Minimum n value for homolog series.
min_n=3
data_filename=read_data(filename_data)
class_find(data_filename,homo=homo,mass_error=error_ms1,rt_error=error_rt,min_n=min_n)
fragment_mark(data_filename,filename_database,error_ms2)