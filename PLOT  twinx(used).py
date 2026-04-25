


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


fig, ax = plt.subplots(nrows=2, ncols=2,figsize=(8,8))

# Plot y(x) on left axi
###################part 1 1
df=pd.read_csv("00.csv")
k=list(set(df["k"]))
k.sort()
# data
J=[]
J_ste=[]
K=[]
K_ste=[]
bs=[]
bs_ste=[]

for item in k:
        boresh=df[df["k"]==item]
        J.append(np.mean(boresh["J"]))
        J_ste.append(np.std(boresh["J"])/np.sqrt(len(boresh)))
        K.append(np.mean(boresh["K"]))
        K_ste.append(np.std(boresh["K"])/np.sqrt(len(boresh)))
        bs.append(np.mean(boresh["std"]))
        bs_ste.append(np.std(boresh["std"])/np.sqrt(len(boresh)))
J=np.array(J)
K=np.array(K)
bs=np.array(bs)
J_ste=np.array(J_ste)
K_ste=np.array(K_ste)
bs_ste=np.array(bs_ste)
J=J/max(J)
K=K/max(K)
J_ste=J_ste/max(J_ste)
K_ste=J_ste/max(K_ste)

########********************########
ax1=ax[0][0]
line1=ax1.plot(k, J,"-o",label="J(A)",markersize=7,linewidth=2,color="#8B3A62")
line2=ax1.plot(k, K,"-s",label="K(A)",markersize=7,linewidth=2,color="#8B4513")
#ax1.set_xlabel('k',fontweight=1000,fontsize=14)
ax1.set_ylabel('K,J',fontweight=1000,fontsize=10)
ax1.tick_params(which="both",direction="in")

# Create right axis for z(x)
ax2 = ax1.twinx()
line3=ax2.errorbar(k, bs, bs_ste,label="Avg bub. size",fmt="-^",capsize=5,markersize=7,color="#008080",linewidth=2)
#ax2.set_ylabel("Average bubble size",fontweight=1000,fontsize=10)
ax2.tick_params(which="both",direction="in")

lines=[line1[0],line2[0],line3[0]]

labels = [line.get_label() for line in lines[:-1]]+["Avg bub. size"]
plt.legend(lines,labels)

for a in [ax1,ax2]:
    for label in a.get_xticklabels() + a.get_yticklabels():
        label.set_fontweight(1000)
        label.set_fontsize(10)  # Optional: adjust size too

plt.title('(a)')


###################part 1 2
df=pd.read_csv("25.csv")
k=list(set(df["k"]))
k.sort()
# data
J=[]
J_ste=[]
K=[]
K_ste=[]
bs=[]
bs_ste=[]

for item in k:
        boresh=df[df["k"]==item]
        J.append(np.mean(boresh["J"]))
        J_ste.append(np.std(boresh["J"])/np.sqrt(len(boresh)))
        K.append(np.mean(boresh["K"]))
        K_ste.append(np.std(boresh["K"])/np.sqrt(len(boresh)))
        bs.append(np.mean(boresh["std"]))
        bs_ste.append(np.std(boresh["std"])/np.sqrt(len(boresh)))
J=np.array(J)
K=np.array(K)
bs=np.array(bs)
J_ste=np.array(J_ste)
K_ste=np.array(K_ste)
bs_ste=np.array(bs_ste)
J=J/max(J)
K=K/max(K)
J_ste=J_ste/max(J_ste)
K_ste=J_ste/max(K_ste)

#############**************##############
ax1=ax[0][1]
line1=ax1.plot(k, J,"-o",label="J(A)",markersize=7,linewidth=2,color="#8B3A62")
line2=ax1.plot(k, K,"-s",label="K(A)",markersize=7,linewidth=2,color="#8B4513")
#ax1.set_xlabel('k',fontweight=1000,fontsize=14)
#ax1.set_ylabel('K,J',fontweight=1000,fontsize=10)
ax1.tick_params(which="both",direction="in")

# Create right axis for z(x)
ax2 = ax1.twinx()
line3=ax2.errorbar(k, bs, bs_ste,label="Avg bub. size",fmt="-^",capsize=5,markersize=7,color="#008080",linewidth=2)
ax2.set_ylabel("Average bubble size",fontweight=1000,fontsize=10)
ax2.tick_params(which="both",direction="in")

lines=[line1[0],line2[0],line3[0]]

labels = [line.get_label() for line in lines[:-1]]+["Avg bub. size"]

for a in [ax1,ax2]:
    for label in a.get_xticklabels() + a.get_yticklabels():
        label.set_fontweight(1000)
        label.set_fontsize(10)  # Optional: adjust size too

plt.title('(b)')


###################part 2 1
df=pd.read_csv("50.csv")
k=list(set(df["k"]))
k.sort()
# data
J=[]
J_ste=[]
K=[]
K_ste=[]
bs=[]
bs_ste=[]

for item in k:
        boresh=df[df["k"]==item]
        J.append(np.mean(boresh["J"]))
        J_ste.append(np.std(boresh["J"])/np.sqrt(len(boresh)))
        K.append(np.mean(boresh["K"]))
        K_ste.append(np.std(boresh["K"])/np.sqrt(len(boresh)))
        bs.append(np.mean(boresh["std"]))
        bs_ste.append(np.std(boresh["std"])/np.sqrt(len(boresh)))
J=np.array(J)
K=np.array(K)
bs=np.array(bs)
J_ste=np.array(J_ste)
K_ste=np.array(K_ste)
bs_ste=np.array(bs_ste)
J=J/max(J)
K=K/max(K)
J_ste=J_ste/max(J_ste)
K_ste=J_ste/max(K_ste)

#############***************##############
ax1=ax[1][0]
line1=ax1.plot(k, J,"-o",label="J(A)",markersize=7,linewidth=2,color="#8B3A62")
line2=ax1.plot(k, K,"-s",label="K(A)",markersize=7,linewidth=2,color="#8B4513")
ax1.set_xlabel('k',fontweight=1000,fontsize=14)
ax1.set_ylabel('K,J',fontweight=1000,fontsize=10)
ax1.tick_params(which="both",direction="in")

# Create right axis for z(x)
ax2 = ax1.twinx()
line3=ax2.errorbar(k, bs, bs_ste,label="Avg bub. size",fmt="-^",capsize=5,markersize=7,color="#008080",linewidth=2)
#ax2.set_ylabel("Average bubble size",fontweight=1000,fontsize=10)
ax2.tick_params(which="both",direction="in")

lines=[line1[0],line2[0],line3[0]]

labels = [line.get_label() for line in lines[:-1]]+["Avg bub. size"]

for a in [ax1,ax2]:
    for label in a.get_xticklabels() + a.get_yticklabels():
        label.set_fontweight(1000)
        label.set_fontsize(10)  # Optional: adjust size too

plt.title('(c)')



###################part 2 2
df=pd.read_csv("75.csv")
k=list(set(df["k"]))
k.sort()
# data
J=[]
J_ste=[]
K=[]
K_ste=[]
bs=[]
bs_ste=[]

for item in k:
        boresh=df[df["k"]==item]
        J.append(np.mean(boresh["J"]))
        J_ste.append(np.std(boresh["J"])/np.sqrt(len(boresh)))
        K.append(np.mean(boresh["K"]))
        K_ste.append(np.std(boresh["K"])/np.sqrt(len(boresh)))
        bs.append(np.mean(boresh["std"]))
        bs_ste.append(np.std(boresh["std"])/np.sqrt(len(boresh)))
J=np.array(J)
K=np.array(K)
bs=np.array(bs)
J_ste=np.array(J_ste)
K_ste=np.array(K_ste)
bs_ste=np.array(bs_ste)
J=J/max(J)
K=K/max(K)
J_ste=J_ste/max(J_ste)
K_ste=J_ste/max(K_ste)

#############***************############
ax1=ax[1][1]
line1=ax1.plot(k, J,"-o",label="J(A)",markersize=7,linewidth=2,color="#8B3A62")
line2=ax1.plot(k, K,"-s",label="K(A)",markersize=7,linewidth=2,color="#8B4513")
ax1.set_xlabel('k',fontweight=1000,fontsize=14)
#ax1.set_ylabel('K,J',fontweight=1000,fontsize=10)
ax1.tick_params(which="both",direction="in")

# Create right axis for z(x)
ax2 = ax1.twinx()
line3=ax2.errorbar(k, bs, bs_ste,label="Avg bub. size",fmt="-^",capsize=5,markersize=7,color="#008080",linewidth=2)
ax2.set_ylabel("Average bubble size",fontweight=1000,fontsize=10)
ax2.tick_params(which="both",direction="in")

lines=[line1[0],line2[0],line3[0]]

labels = [line.get_label() for line in lines[:-1]]+["Avg bub. size"]

for a in [ax1,ax2]:
    for label in a.get_xticklabels() + a.get_yticklabels():
        label.set_fontweight(1000)
        label.set_fontsize(10)  # Optional: adjust size too

plt.title('(d)')
fig.tight_layout()
plt.show()
