#matrix[row][column]
#matrix[i][j] interaction of i with the node j
#matrix[asar pazir][asar gozar] interaction of i with the node j
print("hello again")
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import datetime
def rand():
	return 0.0*(random.random()-.5)

"""
In this varsion I added the kriess constant calculator fumction
as a new property clculator fpr network object
"""
J=[]
J_ste=[]
K=[]
K_ste=[]
bs=[]
bs_ste=[]
mbs=[]

error_mbs=[]
import pandas as pd
df=pd.read_csv("data_full_std10e7_new-v14_test.csv")
p=list(set(df["p"]))
p.sort()
for item in p:
        boresh=df[df["p"]==item]
        J.append(np.mean(boresh["J"]))
        J_ste.append(np.std(boresh["J"])/np.sqrt(len(boresh)))
        K.append(np.mean(boresh["K"]))
        K_ste.append(np.std(boresh["K"])/np.sqrt(len(boresh)))
        bs.append(np.mean(boresh["std"]))
        bs_ste.append(np.std(boresh["std"])/np.sqrt(len(boresh)))


fig, ax = plt.subplots(nrows=1, ncols=3,figsize=(12,4))
ax[0].plot(p,K,"-s",label="K(A)",markersize=7,linewidth=2,color="#8B4513")
ax[0].plot(p,J,"-o",label="J(A)",markersize=7,linewidth=2,color="#8B3A62")
ax[0].set_xlabel("r",fontweight=1000,fontsize=14)
ax[0].set_ylabel("K,J",fontweight=1000,fontsize=10)
ax[0].legend()
ax[1].errorbar(p,bs,bs_ste,fmt="-o",capsize=5,markersize=7,color="#8B3A62",linewidth=2)
ax[1].set_xlabel("r",fontweight=1000,fontsize=14)
ax[1].set_ylabel("Average bubble size",fontweight=1000,fontsize=10)


ax[2].errorbar(J,bs,bs_ste,fmt="-o",capsize=5,markersize=7,color="#8B3A62",linewidth=2)
plt.tick_params(right=True,top=True)
plt.tick_params(which="both",direction="in")
ax[2].set_xlabel("J",fontweight=1000,fontsize=14)
ax[2].set_ylabel("Average bubble size",fontweight=1000,fontsize=10)
#plt.xticks(fontweight=1000)
#plt.yticks(fontweight=1000)
for a in ax:
    for label in a.get_xticklabels() + a.get_yticklabels():
        label.set_fontweight(1000)
        label.set_fontsize(10)  # Optional: adjust size too
ax[0].set_title("(a)")
ax[1].set_title("(b)")
ax[2].set_title("(c)")
plt.tight_layout() 
plt.show()

"""
for edg_prb in p:
	print(edg_prb)
	new_net=network(dt=.05,N_nodes=40,edge_weigth=1.,self_interaction=-.1,search_path="dpath_pos2dpath_neg",edge_prob=edg_prb,need_j=True,need_k=True)
	new_net.mean_bubble_height_calc(n_bubles=800,show_bubble_heights=False)
	J.append(new_net.J)
	K.append(new_net.K)
	mbs.append(new_net.mean_bubble_height[0])
	error_mbs.append(new_net.mean_bubble_height[1])

df=pd.DataFrame({"K":K,"J":J,"mbs":mbs,"error":error_mbs})
df.to_csv("linedata3.csv")
"""



"""
deltat=.05
EW=1.
for k1 in range(100):
	for edg_prb in np.arange(.5,1.000001,.1):
		start = datetime.datetime.now()
		print("making a network with p = "+str(edg_prb))
		#new_net=network(dt=.000002,N_nodes=400,edge_weigth=2400,self_interaction=-1000,search_path="dhaerarcial2haerarcial",edge_prob=edg_prb)
		new_net=network(dt=deltat,N_nodes=500,edge_weigth=EW,self_interaction=-1,search_path="dhaerarcial2haerarcial",edge_prob=edg_prb)
		#new_net.plot_network()
		print(max(new_net.eigenvalues))
		print("network is made")
		print("saving: delta s2 p n flambda J ")
		#saving properries
		#saving the high
		print("calculating the std")
		new_net.calc_std_collective_opinion(N=1000000)
		print("all needed data are gathered!")
		print("1-N: ",new_net.N_nodes)
		print("2-P_edges: ",new_net.edge_prob)
		print("3-J: ",new_net.J)
		print("4-DELTA: ",new_net.DELTA)
		print("5-fast lamda: ",new_net.lambda_fast)
		print("6-S2: ",new_net.m2_S)
		print("7-std:",new_net.std_collective_opinion)
		now = datetime.datetime.now()
		print("this datapoint took:",now-start)
		print("another round is completed!")
		print("saving data in ecxel")
		with open("data_full_std10e7_new-v14_test.csv","a") as file:
			file.write(str(np.real(new_net.J))+","+str(EW)+","+str(np.real(new_net.K))+","+str(new_net.m2_S) +","+str(new_net.N_nodes
	) +","+str(edg_prb)+","+str(np.real(new_net.std_collective_opinion))+"\n")
		print("********************************************","p = " ,str(edg_prb))
		print("********************************************","j = ",new_net.J)
		print("********************************************","std = ",new_net.std_collective_opinion)
		print("*************************************************************************************")
		print("*************************************************************************************")


"""
