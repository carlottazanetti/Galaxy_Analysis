import os
import copy
import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats.stats import pearsonr

###STEP 1
####extracting the parent sample 
names=['z','petroMag_u', 'petroMag_g','petroMag_r', 'petroMag_i', 'petroMag_z', 'h_alpha_flux','h_beta_flux', 'oiii_5007_flux', 'nii_6584_flux','lgm_tot_p50','sfr_tot_p50','absMagU', 'absMagG', 'absMagR', 'absMagI', 'absMagZ']
data_names=['Redshift ','Apparent Magnitude_u ','Apparent Magnitude_g ','Apparent Magnitude_r ','Apparent Magnitude_i ','Apparent Magnitude_z ','Flux HAlpha ','Flux HBeta ','[OIII] ','[NII] ','Mass ','Starformation Rate ','Absolute Magnitude_U ','Absolute Magnitude_G ','Absolute Magnitude_R ','Absolute Magnitude_I ','Absolute Magnitude_Z ']
parent_samp={}
def parent(dizionario,nomi):
    i=0
    while i< len(nomi):
        dizionario[nomi[i]]=values[names[i]]
        i+=1
    return dizionario

####extracting my subsample from the parent sample
My_ID={}
def sub(dizionario,subsample,nomi):
    i=0
    while i< len(nomi):
        subsample[nomi[i]]=dizionario[nomi[i]][my_nameid]
        i+=1
    return subsample

####writing the functions in order to clean the arrays
def nan_inf(array):
    mask=np.logical_and(~np.isnan(array),~np.isinf(array))
    array=array[mask]
    return array,mask

def sigma_clip(array):
    mask2=np.logical_and(array>np.mean(array)-4*np.std(array),array<np.mean(array)+4*np.std(array))
    array=array[mask2]
    return array,mask2

###STEP 2 
####writing a function for the mean value and the error of the redshift
def redshift(array):
    z=array
    z,mask=nan_inf(z)
    h=0
    while h<4:
        z,mask2=sigma_clip(z)
        h+=1
    n=len(z)
    i=0
    j=0
    somma=0
    somma2=0
    while i<n:
        somma=somma+z[i]
        i+=1
    m_z=somma/n
    while j<n:
        somma2=somma2+(z[j]-m_z)**(2)
        j+=1
    sd_z=(somma2/n)**(1/2)
    error=sd_z/((n)**(1/2))

####writing a function to collect the mean,median and standard deviation
def diz(dizionario,nomi):
    dizionario1=copy.deepcopy(dizionario)
    L1=[]
    L2=[]
    L3=[]
    for k in nomi:
        dizionario1[k],mask=nan_inf(dizionario1[k])
        h=0
        while h<4:
            dizionario1[k],mask2=sigma_clip(dizionario1[k])
            h+=1
        L1.append(np.mean(dizionario1[k]))
        L2.append(np.median(dizionario1[k]))
        L3.append(np.std(dizionario1[k]))
    dizionario['Mean ']=np.array(L1)
    dizionario['Median ']=np.array(L2)
    dizionario['Std ']=np.array(L3)
    return dizionario
    
####saving the output results for the mean, median standard deviation
def file(dizionario,filename,arr1,arr2,arr3):
    os.chdir(dir_out)
    if os.path.exists(filename+'.fits') is True:
        print('File '+filename+' already existing ')
    else:
        print('Saving '+filename+' to output')
        col1=pf.Column(name=arr1,format='D',array=dizionario[arr1])
        col2=pf.Column(name=arr2,format='D',array=dizionario[arr2])
        col3=pf.Column(name=arr3,format='D',array=dizionario[arr3])
        cols=pf.ColDefs([col1,col2,col3])
        tbhdu=pf.BinTableHDU.from_columns(cols)
        tbhdu.writeto(filename+'.fits')
    
####Gaussian
def Gauss(bins,mu,sigma):
    x=np.zeros(len(bins)-1)
    for j in range(len(x)):
        x[j]=(bins[j]+bins[j+1])/2
    return x,1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))

####figure setting
plt.rcParams["axes.facecolor"]="white"
plt.rcParams["axes.edgecolor"]="black"
plt.rcParams["axes.linewidth"]=0.8
plt.rcParams["axes.grid"]=True
plt.rcParams["axes.titlesize"]="large"
plt.rcParams["axes.labelsize"]="medium"
plt.rcParams["legend.fontsize"]="small"
plt.rcParams["grid.color"]="0.5"
plt.rcParams["grid.linestyle"]=":"
plt.rcParams["grid.linewidth"]=0.8
plt.rcParams["grid.alpha"]=0.75
plt.rcParams.update({'font.size': 10})

####plotting the histograms for my ID and the parent sample
def hist(dizionario, nomi,directory):
    dizionario1= copy.deepcopy(dizionario)
    k=0
    for i in nomi:
        dizionario1[i],mask=nan_inf(dizionario1[i])
        h=0
        while h<4:
            dizionario1[i],mask2=sigma_clip(dizionario1[i])
            h+=1
        plt.figure(1, figsize=(15,5))
        gs=gridspec.GridSpec(1, 3)
        #mean
        ax1=plt.subplot(gs[0])
        counts1, bins1, ign=ax1.hist(dizionario1[i],density = True, bins=100,histtype='bar',align='mid', edgecolor='k',label=i)
        xm1,mod1=Gauss(bins1,np.mean(dizionario1[i]),np.std(dizionario1[i]))
        ax1.plot(xm1,mod1,c='green',ls='--',lw=2,label='gaussian')
        ax1.axvline(np.mean(dizionario1[i]), color='r',label='mean')
        ax1.axvspan(np.mean(dizionario1[i])-np.std(dizionario1[i]),np.mean(dizionario1[i])+np.std(dizionario1[i]),color='0.3',alpha=0.25)
        ax1.set_xlabel(i)
        plt.legend(loc='upper left',fontsize='x-small')
        #median
        ax2=plt.subplot(gs[1])
        counts2,bins2,ign=ax2.hist(dizionario1[i],density=True,bins=100,histtype='bar',align='mid',edgecolor='k',label=i)
        xm2,mod2=Gauss(bins2,np.mean(dizionario1[i]),np.std(dizionario1[i]))
        ax2.plot(xm2,mod2,c='green',ls='--',lw=2,label='gaussian')
        ax2.axvline(np.median(dizionario1[i]),color='r',label='median')
        ax2.axvspan(np.mean(dizionario1[i])-np.std(dizionario1[i]),np.mean(dizionario1[i])+np.std(dizionario1[i]),color='0.3',alpha=0.25)
        ax2.set_xlabel(i)
        plt.legend(loc='upper left',fontsize='x-small')
        #residuals
        ax3=plt.subplot(gs[2])
        resid1=(counts1-mod1)
        ax3.scatter(xm1,resid1)
        ax3.axhline(resid1.mean(),c='0.3',ls='--')
        ax3.set_ylabel('Residuals')
        ax3.set_xlabel(i)
        os.chdir(directory)
        if os.path.exists(str(i)+'.png') is False :
            while k<1:
                print('Plotting the histograms ...')
                k+=1
            plt.savefig(str(i),dpi=150)
            plt.close()
        else:
            plt.close()

### STEP 3    
###plotting the scatters
data_names2=['Apparent Magnitude_z ','Flux HAlpha ','Mass ','Starformation Rate ','Absolute Magnitude_Z ']
def scatter(dizionario,keyx,keyy):
    k=0
    x=dizionario[keyx]
    y=dizionario[keyy]
    y,mask=nan_inf(y)
    x=x[mask]
    x,mask1=nan_inf(x)
    y=y[mask1]
    h=0
    while h<4:
        y,mask2=sigma_clip(y)
        x=x[mask2]
        x,mask3=sigma_clip(x)
        y=y[mask3]
        h+=1
    res=np.polyfit(x,y,1)
    x_bf=np.arange(x.min(),x.max(),0.01)
    y_bf=np.polyval(res,x_bf)
    plt.figure(1,figsize=(5,5))
    sc=plt.scatter(x,y,marker='o',color='#6666ff',edgecolors='k',zorder=2,label='Data')
    plt.plot(x_bf,y_bf,c='r',ls='--',lw=2)
    plt.ylabel(keyy)
    plt.xlabel(keyx)
    os.chdir(dir_scatter)
    if os.path.exists('Scatter '+str(keyy)+' , '+str(keyx)+'.png') is False:
        while k<1:
            print('Plotting the scatter for '+str(keyy)+'and '+str(keyx))
            k+=1
        plt.savefig('Scatter '+str(keyy)+' , '+str(keyx),dpi=150)
        plt.close()
    else:
        plt.close()

####analyzing correlations and slicing the arrays
subarray={}
def correlation(z,nomi,dizionario,subarray):
    for l in nomi:
        x=dizionario[z]
        j=dizionario[l]
        j,mask1=nan_inf(j)
        x = x[mask1]
        x,mask=nan_inf(x)
        j = j[mask]
        h=0
        while h<4:
            j,mask2=sigma_clip(j)
            x = x[mask2]
            x,mask3=sigma_clip(x)
            j = j[mask3]
            h+=1
        r=pearsonr(x,j)
        if r[0] >= 0.4 or r[0]<= -0.4:
            print('correlation between ',z,'and ',l,'found: Pearson coefficient= ',r[0])
            n=len(j)
            d={}
            if n%2==0:
                m=n/2
                m=int(m)
                j1=j[0:m]
                j2=j[m:n]
                if m%2==0:
                    p=m/2
                    p=int(p)
                    j3=j1[0:p]
                    j4=j1[p:m]
                    j5=j2[0:p]
                    j6=j2[p:m]
                    d[l]=j3,j4,j5,j6
                    for i in range(len(d[l])):
                        subarray['subarray '+str(l)+' '+str(i)]=d[l][i]
                else:
                    q=m/3
                    q=int(q)
                    j7=j1[0:q]
                    j8=j1[q:2*q]
                    j9=j1[2*q:m]
                    j10=j2[0:q]
                    j11=j2[q:2*q]
                    j12=j2[2*q:m]
                    d[l]=j7,j8,j9,j10,j11,j12
                    for i in range(len(d[l])):
                        subarray['subarray '+str(l)+' '+str(i)]=d[l][i]
            else:
                k=m/3
                k=int(k)
                j13=j1[0:k]
                j14=j1[k:2*k]
                j15=j1[2*k:n]
                d[l]=j13,j14,j15
                for i in range(len(d[l])):
                        subarray['subarray '+str(l)+' '+str(i)]=d[l][i]
        else:
            pass
    hist(subarray,subarray.keys(),dir_hist)    
    return subarray

####STEP 4
####defying the color for the colour-mass diagram
color_mass={}
BPT={}
SFR_mass={}
def col(dizionario,arrnew,arr1,arr2):
    dizionario[arrnew]=arr1-arr2
    return dizionario

####defying a function in order to create an array filled with random numbers for the theoretical relation
def random(arr):
    x=np.linspace(np.min(arr),np.max(arr),len(arr))
    return x

####defying the axes for the BPT diagram
def bptxy(diz,x,y,z,q):
    BPTx=[]
    BPTy=[]
    y,mask=nan_inf(y)
    x=x[mask]
    z=z[mask]
    q=q[mask]
    x,mask1=nan_inf(x)
    y=y[mask1]
    z=z[mask1]
    q=q[mask1]
    z,mask2=nan_inf(z)
    x=x[mask2]
    y=y[mask2]
    q=q[mask2]
    q,mask3=nan_inf(q)
    x=x[mask3]
    y=y[mask3]
    z=z[mask3]
    k=0
    while k<4:
        y,mask4=sigma_clip(y)
        x=x[mask4]
        z=z[mask4]
        q=q[mask4]
        x,mask5=sigma_clip(x)
        y=y[mask5]
        z=z[mask5]
        q=q[mask5]
        z,mask6=sigma_clip(z)
        x=x[mask6]
        y=y[mask6]
        q=q[mask6]
        q,mask7=sigma_clip(q)
        x=x[mask7]
        y=y[mask7]
        z=z[mask7]
        k+=1
    for a in range(len(x)):
        if x[a]*y[a]<=0 or z[a]*q[a]<=0:
            pass
        else:
            BPTx.append(np.log10(x[a]/y[a]))
            BPTy.append(np.log10(z[a]/q[a]))
    diz['log([OIII]/HBeta) ']=np.array(BPTx)
    diz['log([NII]/HAlpha) ']=np.array(BPTy)
    return diz
 
####writing a function for each theoretical relation 
def relCM(dizionario,arr1):
    dizionario['mass ']=arr1
    x=dizionario['mass ']
    u_rx=random(x)
    r=-0.495+0.25*u_rx
    dizionario['y_relation']=r
    dizionario['x_relation']=u_rx
    return dizionario

def relSFR(dizionario,arr1,arr2):
    dizionario['mass ']=arr1
    dizionario['starformation rate ']=arr2
    x=dizionario['mass ']
    SFRx=random(x)
    r=-8.64+0.76*SFRx
    dizionario['y_relation']=r   
    dizionario['x_relation']=SFRx
    return dizionario

def relBPT(dizionario,arr1):
    x=dizionario[arr1]
    a=random(x)
    a=a[a>0]
    b=np.linspace(1,np.max(x),len(x)-len(a))    
    for i in b:
        a=np.append(a,i)
    r=0.61/(np.log10(a)-0.05)+1.3
    dizionario['y_relation']=r
    dizionario['x_relation']=np.log10(a)
    return dizionario

####writing a function for the colour-mass and the SFR-mass diagram
def scatter3(dizionario,keyx,keyy,arrz,t_rx,t_ry,name):
    k=0
    h=0
    x=dizionario[keyx]
    y=dizionario[keyy]
    z=arrz
    t_rx=dizionario[t_rx]
    t_ry=dizionario[t_ry]
    y,mask=nan_inf(y)
    x=x[mask]
    z=z[mask]
    x,mask1=nan_inf(x)
    y=y[mask1]
    z=z[mask1]
    z,mask2=nan_inf(z)
    x=x[mask2]
    y=y[mask2]
    while h<4:
        y,mask3=sigma_clip(y)
        x=x[mask3]
        z=z[mask3]
        x,mask4=sigma_clip(x)
        y=y[mask4]
        z=z[mask4]
        z,mask5=sigma_clip(z)
        x=x[mask5]
        y=y[mask5]
        h+=1    
    plt.figure(1,figsize=(5,5))
    sc=plt.scatter(x,y,marker='o',c=z,edgecolor='k',zorder=2,label='Data')
    plt.plot(t_rx,t_ry,color='r',linestyle='-')
    plt.xlim(8,12)
    cb=plt.colorbar(sc)
    cb.set_label('Redshift')
    plt.ylabel(keyy)
    plt.xlabel(keyx)   
    os.chdir(dir_scatter)
    if os.path.exists(name+'.png') is False:
        while k<1:
            print('Plotting the scatter for '+name+'...')
            k+=1
        plt.savefig(name,dpi=150)
        plt.close()
    else:
        plt.close()



####writing a function to plot the BPT diagram
def scatterBPT(dizionario,keyx,keyy,t_rx,t_ry,name):
    k=0
    h=0
    x=dizionario[keyx]
    y=dizionario[keyy]
    t_rx=dizionario[t_rx]
    t_ry=dizionario[t_ry]    
    y,mask=nan_inf(y)
    x=x[mask]
    x,mask1=nan_inf(x)
    y=y[mask1]
    while h<4:
        y,mask3=sigma_clip(y)
        x=x[mask3]
        x,mask4=sigma_clip(x)
        y=y[mask4]
        h+=1
    plt.figure(1,figsize=(5,5))
    sc=plt.scatter(x,y,marker='.',color='k',zorder=2,label='Data')
    plt.plot(t_rx,t_ry,color='r',linestyle='-')
    plt.ylim(-2,2)
    plt.xlim(-1.5,0.5)
    plt.ylabel(keyy)
    plt.xlabel(keyx)
    os.chdir(dir_scatter)
    if os.path.exists(name+'.png') is False:
        while k<1:
            print('Plotting the scatter for '+name+'...')
            k+=1
        plt.savefig(name,dpi=150)
        plt.close()
    else:
        plt.close() 

####writing a function in order to slice the arrays according to the theoretical relation
def slice_arr(diz,arr1,arr2,array_relation,nome_sopra_rel,nome_sotto_rel):
    r=diz[array_relation]
    x=diz[arr1]
    y=diz[arr2]
    l1=[]
    l2=[]
    l3=[]
    l4=[]
    for i in range(len(x)):
        if y[i] < r[i]:          
            l1.append(x[i])
            l2.append(y[i])
            i+=1
        else:
            l3.append(x[i])
            l4.append(y[i])
            i+=1
    diz[nome_sopra_rel+str(arr1)]=np.array(l1)
    diz[nome_sopra_rel+str(arr2)]=np.array(l2)
    diz[nome_sotto_rel+str(arr1)]=np.array(l3)
    diz[nome_sotto_rel+str(arr2)]=np.array(l4)
    return diz
               
#plotting histograms for color-mass,str-mass,bpt
def hist3(dizionario,i,k,filename):
    dizionario1= copy.deepcopy(dizionario)
    dizionario1[i],mask=nan_inf(dizionario1[i])
    dizionario1[k],mask1=nan_inf(dizionario1[k])
    f=0
    h=0
    while h<4:
        dizionario1[i],mask2=sigma_clip(dizionario1[i])
        dizionario1[k],mask3=sigma_clip(dizionario1[k])
        h+=1
    plt.figure(1, figsize=(10,5))
    gs=gridspec.GridSpec(1, 2)
        
    ax1=plt.subplot(gs[0])
    counts1, bins1, ign=ax1.hist(dizionario1[i],density = True, bins=100,histtype='bar',align='mid', edgecolor='k',label=i)
    xm1,mod1=Gauss(bins1,np.mean(dizionario1[i]),np.std(dizionario1[i]))
    ax1.plot(xm1,mod1,c='green',ls='--',lw=2,label='gaussian')
    ax1.axvline(np.mean(dizionario1[i]), color='r',label='mean')
    ax1.axvspan(np.mean(dizionario1[i])-np.std(dizionario1[i]),np.mean(dizionario1[i])+np.std(dizionario1[i]),color='0.3',alpha=0.25)
    ax1.axvline(np.median(dizionario1[i]),color='k',label='median')
    ax1.set_xlabel(i)
    plt.legend(loc='upper left',fontsize='x-small')

    ax2=plt.subplot(gs[1])
    counts2, bins2, ign=ax2.hist(dizionario1[k],density = True, bins=100,histtype='bar',align='mid', edgecolor='k',label=k)
    xm2,mod2=Gauss(bins2,np.mean(dizionario1[k]),np.std(dizionario1[k]))
    ax2.plot(xm2,mod2,c='green',ls='--',lw=2,label='gaussian')
    ax2.axvline(np.mean(dizionario1[k]), color='r',label='mean')
    ax2.axvspan(np.mean(dizionario1[k])-np.std(dizionario1[k]),np.mean(dizionario1[k])+np.std(dizionario1[k]),color='0.3',alpha=0.25)
    ax2.axvline(np.median(dizionario1[k]),color='k',label='median')
    ax2.set_xlabel(k)
    plt.legend(loc='upper left',fontsize='x-small')
    os.chdir(dir_hist)
    if os.path.exists(filename+'.png') is False :
        while f<1:
            print('Plotting the histograms for '+filename+ '...')
            f+=1
        plt.savefig(filename,dpi=150)
        plt.close()
    else:
        plt.close()

color_mass_names=['Red Galaxies mass ', 'Red Galaxies colour', 'Blue Galaxies mass ', 'Blue Galaxies colour']
SFR_mass_names=['Starforming Galaxies mass ', 'Starforming Galaxies starformation rate ', 'Passive Galaxies mass ', 'Passive Galaxies starformation rate ']
BPT_names=['Upper Galaxies log([OIII]/HBeta) ', 'Upper Galaxies log([NII]/HAlpha) ', 'Lower Galaxies log([OIII]/HBeta) ', 'Lower Galaxies log([NII]/HAlpha) ']
def file2(dizionario,filename,arr1,arr2,arr3,arr4):
    os.chdir(dir_out)
    if os.path.exists(filename+'.fits') is True:
        print('File '+filename+' already existing ')
    else:
        print('Saving '+filename+' to output')
        col1=pf.Column(name=arr1,format='D',array=dizionario[arr1])
        col2=pf.Column(name=arr2,format='D',array=dizionario[arr2])
        col3=pf.Column(name=arr3,format='D',array=dizionario[arr3])
        col4=pf.Column(name=arr4,format='D',array=dizionario[arr4])
        cols=pf.ColDefs([col1,col2,col3])
        tbhdu=pf.BinTableHDU.from_columns(cols)
        tbhdu.writeto(filename+'.fits')


#creating working directories
os.chdir(os.getcwd())
dir_home=os.getcwd()
dir_data=dir_home+"/data/"
dir_out=dir_home+"/output/"
dir_plots=dir_home+'/plots/'
dir_hist=dir_plots+'/histograms/'
dir_scatter=dir_plots+'/scatters/'
dir_parent=dir_plots+'/parent sample/'
    
if  os.path.exists(dir_data) is False :
    os.mkdir(dir_data)
else:
    print("Directory dir_data already existing. ")
if os.path.exists(dir_out) is False:
    os.mkdir(dir_out)
else:
    print("Directory dir_out already existing. ")
if  os.path.exists(dir_plots) is False :
    os.mkdir(dir_plots)
else:
    print("Directory dir_plots already existing. ")
if os.path.exists(dir_hist) is False:
    os.mkdir(dir_hist)
else:
    print('Directory histograms already existing. ')
if os.path.exists(dir_scatter) is False:
    os.mkdir(dir_scatter)
else:
    print('Directory scatters already existing. ')
if os.path.exists(dir_parent) is False:
    os.mkdir(dir_parent)
else:
    print('Directory parent sample already existing. ')

#moving data_SDSS_Info
os.chdir(dir_data)
# Specify the absolute path to the file
file_path = r"C:\Users\carlo\Documents\tesine triennale\ProjectPython\Project\data_SDSS_Info.fit"
if os.path.exists(file_path) is False:
    print('yea')
    file_name = 'data_SDSS_Info.fit'
    cmd= f"move {file_name} {dir_home}"
    os.system(cmd)
else:
    pass
    

####opening the file
hdul=pf.open(file_path)
hdul.info()
print(hdul[1].header)
cols=hdul[1].columns
cols.info()
values=hdul[1].data
nameid=values['ID']
my_nameid=values['ID']==27
#print(values[my_nameid])




#####functions#########
parent(parent_samp,data_names)
sub(parent_samp,My_ID,data_names)
redshift(My_ID['Redshift '])
hist(parent_samp,data_names,dir_parent)
hist(My_ID,data_names,dir_hist)
diz(My_ID,data_names)
file(My_ID,'My_ID','Mean ','Median ','Std ')
diz(parent_samp,data_names)
file(parent_samp,'parent_sample','Mean ','Median ','Std ')
for i in data_names2:
    scatter(My_ID,'Redshift ',i)
correlation('Redshift ',data_names2,My_ID,subarray)
diz(subarray,subarray.keys())
file(subarray,'Subarrays','Mean ','Median ','Std ')
col(color_mass,'colour',My_ID['Apparent Magnitude_u '],My_ID['Apparent Magnitude_r '])
bptxy(BPT,My_ID['[NII] '],My_ID['Flux HAlpha '],My_ID['[OIII] '],My_ID['Flux HBeta '])
relCM(color_mass,My_ID['Mass '])
relBPT(BPT,'log([OIII]/HBeta) ')
relSFR(SFR_mass,My_ID['Mass '],My_ID['Starformation Rate '])
scatter3(color_mass,'mass ','colour',My_ID['Redshift '],'x_relation','y_relation','color-mass ')
scatter3(SFR_mass,'mass ','starformation rate ',My_ID['Redshift '],'x_relation','y_relation','SFR-mass')
scatterBPT(BPT,'log([OIII]/HBeta) ','log([NII]/HAlpha) ','x_relation','y_relation','BPT')
slice_arr(color_mass,'mass ','colour','y_relation','Red Galaxies ','Blue Galaxies ')
slice_arr(SFR_mass,'mass ','starformation rate ','y_relation','Starforming Galaxies ','Passive Galaxies ')
slice_arr(BPT,'log([OIII]/HBeta) ','log([NII]/HAlpha) ','y_relation','Upper Galaxies ','Lower Galaxies ')
hist3(color_mass,'Blue Galaxies mass ','Red Galaxies mass ','Blue-Red Galaxies mass ')
hist3(color_mass,'Blue Galaxies colour','Red Galaxies colour','Blue-Red Galaxies colour ')
hist3(SFR_mass,'Passive Galaxies mass ','Starforming Galaxies mass ','Passive-Starforming Galaxies Mass ')
hist3(SFR_mass,'Passive Galaxies starformation rate ','Starforming Galaxies starformation rate ','Passive-Starforming Galaxies Staformation Rate ')
hist3(BPT,'Upper Galaxies log([OIII]/HBeta) ','Lower Galaxies log([OIII]/HBeta) ','Uperr-Lower Galaxies log([OIII]:HBeta) ')
hist3(BPT,'Upper Galaxies log([NII]/HAlpha) ','Lower Galaxies log([NII]/HAlpha) ','Upper-Lower Galaxies log([NII]:HAlpha) ')
diz(color_mass,color_mass_names)
diz(SFR_mass,SFR_mass_names)
diz(BPT,BPT_names)
file2(color_mass,'Color-mass diagram','Mean ','Median ','Std ','y_relation')
file2(SFR_mass,'SFR-mass diagram','Mean ','Median ','Std ','y_relation')
file2(BPT,'BPT-diagram','Mean ','Median ','Std ','y_relation')














