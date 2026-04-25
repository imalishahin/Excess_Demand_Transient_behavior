#matrix[row][column]
#matrix[i][j] interaction of i with the node j
#matrix[asar pazir][asar gozar] interaction of i with the node j
print("hello again")
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import datetime

"""
in this version we want to see add a serch path avout adding negetivly wigted
edges to heirarkies
"""
class network():
	def __init__(self,dt=.01,N_nodes=1000,noise_amp=1.0,search_path="dstar2completedag",edge_weigth=1.0,edge_prob=0.0,self_interaction=-1.0,need_j=True,need_k=False):
		self.dt=dt
		self.edge_prob=edge_prob
		self.N_nodes=N_nodes
		self.self_interaction=self_interaction
		self.edge_weigth=edge_weigth
		self.noise_amp=noise_amp
		#print("this network is " , (1-edge_prob)*100 ,"percent similar to a star")
		#####making the adjacency matrix
		
		self.state=self.N_nodes*[0.0]
		self.collective_opinion=0.0
		self.eigenvalues=(1,2)
		whilecounter=0
		while max(np.real(self.eigenvalues))>0: 
			self.make_the_network(search_path)
			self.dynamic_matrix=self.adjacency_matrix.copy()		
			for i in range(self.N_nodes):
				self.dynamic_matrix[i][i]=-sum(np.abs(self.adjacency_matrix[i]))+self.self_interaction
		#print(self.dynamic_matrix)
			self.eigenvalues, self.eigenvectors = np.linalg.eig(self.dynamic_matrix)
			whilecounter+=1
			if whilecounter>100:
				self.eigenvalues=None
				break
		self.eigenvectors=np.transpose(self.eigenvectors)
		#sorting eigs
		#idx=np.argsort(self.eigenvalues)
		#self.eigenvalues=self.eigenvalues[idx]
		#self.eigenvectors = self.eigenvectors[:, idx]
		#######end sorting eigs
		#######adding super positions
		non_repeat_eigs=set(self.eigenvalues)
		#print("eigen vectors\n   ",np.round(self.eigenvectors,3))
		#print("eigen vectors\n   ",np.round(self.eigenvalues,3))
		for item in non_repeat_eigs:
			idxs=np.where(self.eigenvalues==item)
			if len(idxs[0])==1:
				continue
			eigvec_superpos=np.sum(self.eigenvectors[idxs], axis=0)
			norm = np.linalg.norm(eigvec_superpos)
			eigvec_superpos=eigvec_superpos/norm
			self.eigenvectors = np.append(self.eigenvectors,  np.array([eigvec_superpos]), axis=0)
			self.eigenvalues=np.append(self.eigenvalues, item)

		#print("eigen vectors\n   ",np.round(self.eigenvectors,3))
		#print("eigen vectors\n   ",np.round(self.eigenvalues,3))
		#######end adding super positions
		self.lambda_fast=min(self.eigenvalues)
		self.lambda_vec=(self.eigenvalues/self.lambda_fast)
		self.m2_lambdas=np.dot(self.lambda_vec,self.lambda_vec)
		self.m1_lambdas=sum(self.lambda_vec)
		self.var_lambdas=(self.m2_lambdas)/self.N_nodes-(self.m1_lambdas/self.N_nodes*self.m1_lambdas/self.N_nodes)
	

		#print(np.transpose(self.eigenvectors))
		U=np.array(self.N_nodes*[1.0])/np.sqrt(self.N_nodes)
		cos_theta=[]
		for i in range(len(self.eigenvectors)):
			cos_theta.append(np.dot(self.eigenvectors[i],U))
		self.cos_theta=np.array(cos_theta)
		#print(self.cos_theta)
		self.DELTA=self.N_nodes*np.dot(self.cos_theta,self.cos_theta)
		
		self.L_vector=self.cos_theta*self.lambda_vec
		self.m2_L=np.dot(self.L_vector,self.L_vector)
		self.m1_L=sum(self.L_vector)
		self.var_L=self.m2_L-self.m1_L*self.m1_L/self.N_nodes
		self.S_vector=[]
		for i in range(self.N_nodes):
			S_=sum(self.adjacency_matrix[:,i])-sum(self.adjacency_matrix[i,:])+self.self_interaction
			self.S_vector.append(S_)
		self.S_vector=np.array(self.S_vector)
		self.m2_S=np.dot(self.S_vector,self.S_vector)
		self.min_arc=np.round(np.arccos(np.sqrt((self.N_nodes-1)/self.N_nodes)),3)-.001

		if need_j:
			self.calculate_J()
		if need_k:
			self.calculate_K()
		else:
			self.K="-"
	def make_the_network(self,search_path):
		self.adjacency_matrix = np.array(self.N_nodes*[self.N_nodes*[0.0]])
		if search_path=="star2complete":
			for i in range(1,self.N_nodes):
				self.adjacency_matrix[i][0]=self.edge_weigth
			for i in range(1,self.N_nodes):
				for j in range(self.N_nodes):
					if (i != j) and (random.random()<=self.edge_prob):
						self.adjacency_matrix[j][i]=self.edge_weigth
		elif search_path=="dstar2completedag":
			for i in range(1,self.N_nodes):
				self.adjacency_matrix[i][0]=self.edge_weigth
			for i in range(1,self.N_nodes):
				for j in range(i+1,self.N_nodes):
					if (random.random()<=self.edge_prob):
						self.adjacency_matrix[j][i]=self.edge_weigth
		elif search_path=="dstar2star":
			for i in range(1,self.N_nodes):
				self.adjacency_matrix[i][0]=self.edge_weigth
			for i in range(1,self.N_nodes):
				if (random.random()<=self.edge_prob):
					self.adjacency_matrix[0][i]=self.edge_weigth
		elif search_path=="dstar2star_nedg":
			for i in range(1,self.N_nodes):
				self.adjacency_matrix[i][0]=self.edge_weigth
			for i in range(1,1+round((self.N_nodes-1)*self.edge_prob)):
				self.adjacency_matrix[0][i]=self.edge_weigth
		elif search_path=="dstar2star_nedg_2c":
			#maybe 2c means 2 centers i dont really remeber why i made it.
			for i in range(2,self.N_nodes):
				self.adjacency_matrix[i][0]=self.edge_weigth
				self.adjacency_matrix[i][1]=self.edge_weigth
			for i in range(2,int(self.N_nodes*self.edge_prob)):
				self.adjacency_matrix[0][i]=self.edge_weigth
				self.adjacency_matrix[1][i]=self.edge_weigth
		elif search_path=="dpath2completedag":
			for i in range(0,self.N_nodes-1):
				self.adjacency_matrix[i+1][i]=self.edge_weigth
			for i in range(1,self.N_nodes):
				for j in range(i+1,self.N_nodes):
					if (random.random()<=self.edge_prob):
						self.adjacency_matrix[j][i]=self.edge_weigth
		elif search_path=="dpath2path":
			for i in range(0,self.N_nodes-1):
				self.adjacency_matrix[i+1][i]=self.edge_weigth
			for i in range(0,round(self.edge_prob*self.N_nodes)-1):
				self.adjacency_matrix[self.N_nodes-2-i][self.N_nodes-1-i]=self.edge_weigth
		elif search_path=="dhaerarcial2haerarcial":
			rank_dic={}
			rank_dic[0]=0
			for i in range(1,self.N_nodes):
				choosing_list=[]
				for j in range(0,i):
					score=int(sum(self.adjacency_matrix[:,j]))+1
					choosing_list=choosing_list+score*[j]
				source=random.choice(choosing_list)
				self.adjacency_matrix[i][source]=self.edge_weigth
				rank_dic[i]=rank_dic[source]+1
			###recirocity probability
				a=100
				b=0
				theta=self.edge_prob#.21
				l=rank_dic[i]
				p=theta/(1+np.exp(-a*(l-b)))#-(theta/(1+np.exp(-a*(1-b))))
				#print(l,"-----",p)
				if random.random()<p:
					self.adjacency_matrix[source][i]=self.edge_weigth
		
		elif search_path=="dpath_pos2dpath_neg":
			for i in range(0,self.N_nodes-1):
				self.adjacency_matrix[i+1][i]=self.edge_weigth
				#self.adjacency_matrix[i][i+1]=self.edge_weigth/5*0
			for i in range(0,round(self.edge_prob*self.N_nodes)-1):
				self.adjacency_matrix[i+1][i]=-self.adjacency_matrix[i+1][i]
				#self.adjacency_matrix[i][i+1]=-self.adjacency_matrix[i][i+1]
		elif search_path=="dhaerarcial2negetivedhaerarcial":
			rank_dic={}
			rank_dic[0]=0
			for i in range(1,self.N_nodes):
				choosing_list=[]
				for j in range(0,i):
					score=int(sum(np.abs(self.adjacency_matrix[:,j])))+1
					choosing_list=choosing_list+score*[j]
				source=random.choice(choosing_list)
				if random.random()<self.edge_prob:
					randomvaribale=-1
				else:
					randomvaribale=1
				self.adjacency_matrix[i][source]=randomvaribale*self.edge_weigth
				rank_dic[i]=rank_dic[source]+1

		
		elif search_path=="dpath_pos2dpath_neg":
			for i in range(0,self.N_nodes-1):
				self.adjacency_matrix[i+1][i]=self.edge_weigth
				#self.adjacency_matrix[i][i+1]=self.edge_weigth/5*0
			for i in range(0,round(self.edge_prob*self.N_nodes)-1):
				self.adjacency_matrix[i+1][i]=-self.adjacency_matrix[i+1][i]
				#self.adjacency_matrix[i][i+1]=-self.adjacency_matrix[i][i+1]
		
###
	def calculate_K(self,square_side=10):
	    #this function searchs for the maximum on a square of size  square_side   
	    N=self.N_nodes
	    i=0+1j
	    val=0
	    for x in np.arange(.0001,square_side,.05):
		    
		    y=0
		    z=x+y*i
		    nval=(x*np.linalg.norm(np.linalg.inv(z*np.identity(N)-self.dynamic_matrix),ord=2))
		    if nval>val:
			    val=nval
		    #print(x,val)	    	
				    
	    self.K=val					
###

	def calculate_J(self,ring_segments_number=30):
		M_max=0
		Neigens=len(self.eigenvalues)
		for icj in range(Neigens):
			sqrtN_times_costhetai=sum(self.eigenvectors[icj])
			if icj%100==0:
				print("the circle we are searching in: ",icj)
			for jcj in range(icj+1,Neigens):
				if ((self.eigenvalues[jcj]-self.eigenvalues[icj]))==0:
					#print("repeating eigenvector skiped")
					continue
				alpha_ij=np.arccos(np.dot(self.eigenvectors[icj],self.eigenvectors[jcj]))
				
				if ((abs(alpha_ij) < self.min_arc) or (abs(alpha_ij-np.pi) < self.min_arc)) :
					#print("repiting eigenvector",alpha_ij)
					continue
				sqrtN_times_costhetaj=sum(self.eigenvectors[jcj])
				lmdj_minus_lmdi_vr_lmj=(self.eigenvalues[jcj]-self.eigenvalues[icj])/self.eigenvalues[jcj]
				lmdi_vr_lmdj_minus_lmi=self.eigenvalues[icj]/(self.eigenvalues[jcj]-self.eigenvalues[icj])
				costi_li_vr_costj_lj=sqrtN_times_costhetai*self.eigenvalues[icj]/sqrtN_times_costhetaj/self.eigenvalues[jcj]
				
				for alpha in np.arange(0,np.pi,np.pi/ring_segments_number):
					c_i=np.cos(alpha)-np.sin(alpha)/np.tan(alpha_ij)
					#c_j=np.cos(alpha-alpha_ij)+np.sin(alpha-alpha_ij)/np.tan(alpha_ij)
					t_star=np.log(costi_li_vr_costj_lj*np.sin(alpha-alpha_ij)/np.sin(alpha))/(self.eigenvalues[jcj]-self.eigenvalues[icj])
					if np.isnan(t_star) or t_star<0:
						continue
					M=sqrtN_times_costhetai*lmdj_minus_lmdi_vr_lmj*c_i*(costi_li_vr_costj_lj*np.sin(alpha-alpha_ij)/np.sin(alpha))**lmdi_vr_lmdj_minus_lmi
					if M>500:
						print("catch one!")
						print("J value: " , M)
						print("T star: ", t_star)
						print("sqrtN_times_costhetai: ",sqrtN_times_costhetai)
						print("lmdj_minus_lmdi_vr_lmj: ",lmdj_minus_lmdi_vr_lmj)
						print(" c_i: ",c_i)
						print("costi_li_vr_costj_lj : ",costi_li_vr_costj_lj)
						print("if zero np.sin(alpha) : ",np.sin(alpha))
						print("lmdi_vr_lmdj_minus_lmi : ",lmdi_vr_lmdj_minus_lmi)
						print(" alpha_ij: ",alpha_ij)
						print("their idxs icj: ", icj)
						print("their idxs: jcj", jcj)
						print(" self.eigenvectors[icj]: ",np.round(self.eigenvalues[icj],6))
						print("self.eigenvectors[jcj]: ",np.round(self.eigenvalues[jcj],6))
						print(" self.eigenvectors[icj]: ",np.round(self.eigenvectors[icj],6))
						print("self.eigenvectors[jcj]: ",np.round(self.eigenvectors[jcj],6))
						print("eigen vectors:\n", np.round(self.eigenvectors,3))
						print("np.tan(alpha_ij) : ",np.tan(alpha_ij))
						input("press enter to continue ")
					if abs(M)>M_max:
						M_max=abs(M)                                        
					#print(icj,jcj,M)
					#term1and5=sum(self.eigenvectors[icj])
					#term2=
					#print(term2)
		self.J=M_max
		#plt.plot(ms)
		#plt.show()
		#print("i can work sir")	
		#input(":/")
	def reset_state(self):
		self.state=self.N_nodes*[0.0]
		self.collective_opinion=0.0
	def update_state_onetimestep(self):
		self.state+=np.matmul(self.dynamic_matrix,self.state)*self.dt
		self.collective_opinion=sum(self.state)
	def update_state_onetimestep_langevian(self):
		self.state+=((np.matmul(self.dynamic_matrix,self.state)*self.dt)+(np.sqrt(self.dt)*self.noise_amp*np.random.normal(size=self.N_nodes)))
		self.collective_opinion=sum(self.state)
	def plot_state_for_N_ts(self,leader_initial_state=1.0,N=1000,initial_cond="central"):#plot the path of the state for N time steps
		if initial_cond=="central":
			self.state[0]=leader_initial_state
		elif initial_cond=="random":
			self.state=random_array = np.random.uniform(-1, 1, self.N_nodes)
		path=[]
		path.append(self.state.copy())
		for t in range(N):
			self.update_state_onetimestep()
			path.append(self.state.copy())
		path=np.array(path)
		plt.plot(path)
		plt.show()
		self.reset_state()
	def plot_state_for_N_ts_langevian(self,leader_initial_state=0.0,N=1000,initial_cond="central"):#plot the path of the state for N time steps
		if initial_cond=="central":
			self.state[0]=leader_initial_state
		elif initial_cond=="random":
			self.state=random_array = np.random.uniform(-1, 1, self.N_nodes)
		path=[]
		path.append(self.state.copy())
		for t in range(N):
			self.update_state_onetimestep_langevian()
			path.append(self.state.copy())
		path=np.array(path)
		plt.plot(path)
		plt.show()
		self.reset_state()
	def plot_collective_opinion_for_N_ts(self,leader_initial_state=1.0,N=1000,initial_cond="central"):#plot the path of the state for N time steps
		if initial_cond=="central":
			self.state[0]=leader_initial_state
		elif initial_cond=="random":
			self.state=random_array = np.random.uniform(-1, 1, self.N_nodes)
		path=[1.0]
		for t in range(N-1):
			self.update_state_onetimestep()
			path.append(self.collective_opinion)
		time=self.dt*np.array((range(N)))
		plt.plot(time,path)
		plt.show()
		self.reset_state()
	def plot_collective_opinion_for_N_ts_langevian(self,leader_initial_state=0.0,N=1000,initial_cond="central"):#plot the path of the state for N time steps
		if initial_cond=="central":
			self.state[0]=leader_initial_state
		elif initial_cond=="random":
			self.state=random_array = np.random.uniform(-1, 1, self.N_nodes)
		path=[0.0]
		for t in range(N-1):
			self.update_state_onetimestep_langevian()
			path.append(self.collective_opinion)
		time=self.dt*np.array((range(N)))
		plt.plot(time,path)
		plt.show()
		self.reset_state()
	def calc_max_collective_opinion(self,leader_initial_state=1.0,N=1000,initial_cond="central"):#plot the path of the state for N time steps
		if initial_cond=="central":
			self.state[0]=leader_initial_state
		elif initial_cond=="random":
			self.state=random_array = np.random.uniform(-1, 1, self.N_nodes)
			self.state[0]=leader_initial_state
		max_collective_opinion=abs(sum(self.state))
		for t in range(N):
			self.update_state_onetimestep()
			if abs(self.collective_opinion)>max_collective_opinion:
				max_collective_opinion=abs(self.collective_opinion)
		self.max_collective_opinion=max_collective_opinion
		self.reset_state()
	def calc_std_collective_opinion(self,leader_initial_state=0.0,N=100000,initial_cond="central"):#plot the path of the state for N time steps
		if initial_cond=="central":
			self.state[0]=leader_initial_state
		elif initial_cond=="random":
			self.state=random_array = np.random.uniform(-1, 1, self.N_nodes)
		sum_col_opi_pow2=leader_initial_state**2
		for t in range(N):
			self.update_state_onetimestep_langevian()
			sum_col_opi_pow2+=(self.collective_opinion**2)
			if abs(self.collective_opinion)>100000:
				sum_col_opi_pow2=np.inf
				break
		self.std_collective_opinion=np.sqrt(sum_col_opi_pow2/(N+1))
		self.reset_state()
	def mean_bubble_height_calc_2see(self,n_bubles=10,leader_initial_state=0.0,N=1000000,initial_cond="central"):
		if initial_cond=="central":
			self.state[0]=leader_initial_state
			self.state= np.random.uniform(-1, 1, self.N_nodes)
		elif initial_cond=="random":
			self.state= np.random.uniform(-1, 1, self.N_nodes)
		co=[]
		co.append(self.collective_opinion)
		self.update_state_onetimestep_langevian()
		sign_bubble_start=np.sign(self.collective_opinion)
		co.append(self.collective_opinion)
		bubble_counter=0
		bubble_heights=[]
		bubble_height=0
		for t in range(N):
			self.update_state_onetimestep_langevian()
			co.append(self.collective_opinion)

			if (sign_bubble_start*np.sign(self.collective_opinion)<=0):
				bubble_counter+=1
				sign_bubble_start*=(-1)
				bubble_heights.append(bubble_height)
				bubble_height=float(abs(self.collective_opinion))
				if (bubble_counter>=n_bubles):
					break
			if abs(self.collective_opinion)>bubble_height:
				bubble_height=float(abs(self.collective_opinion))
		print(bubble_heights)
		plt.plot(co,"-")
		plt.plot(len(co)*[0])
		plt.show()
		self.reset_state()
	def mean_bubble_height_calc(self,n_bubles=200,show_bubble_heights=False,leader_initial_state=0.0,N=1000000,initial_cond="central"):
		if initial_cond=="central":
			self.state[0]=leader_initial_state
			self.state= np.random.uniform(-1, 1, self.N_nodes)
		elif initial_cond=="random":
			self.state= np.random.uniform(-1, 1, self.N_nodes)
		self.update_state_onetimestep_langevian()
		sign_bubble_start=np.sign(self.collective_opinion)
		bubble_counter=0
		bubble_heights=[]
		bubble_height=0
		for t in range(N):
			self.update_state_onetimestep_langevian()

			if (sign_bubble_start*np.sign(self.collective_opinion)<=0):
				bubble_counter+=1
				sign_bubble_start*=(-1)
				bubble_heights.append(bubble_height)
				bubble_height=(abs(self.collective_opinion))
				if (bubble_counter>=n_bubles):
					break
			if abs(self.collective_opinion)>bubble_height:
				bubble_height=(abs(self.collective_opinion))
		if show_bubble_heights:
			plt.plot(bubble_heights,".")
			plt.show()

		self.mean_bubble_height=(np.mean(bubble_heights) , np.std(bubble_heights)/np.sqrt(n_bubles))
	
	def plot_network(self):
		# Create a directed graph from the adjacency matrix
		G = nx.from_numpy_array(np.transpose(self.adjacency_matrix), create_using=nx.DiGraph)
		if nx.is_planar(G):
			pos = nx.planar_layout(G)
			pos[0]=[0,0]
			for node in list(G.nodes):
				pos[node]=[pos[node][0],-nx.shortest_path_length(G, source=0, target=node)-1]

			
			nx.draw(G, pos, with_labels=True,node_size=70, node_color='lightblue', font_size=10, font_weight='bold')
			plt.show()
		else:
			# Draw the directed graph
			pos = nx.circular_layout(G)  # positions for all nodes
			pos[0]=(0,0)
			for i in range(1,self.N_nodes):
				pos[i]=(np.cos((i-1)*6.28/(self.N_nodes-1)),np.sin((i-1)*6.28/(self.N_nodes-1)))
			nx.draw(G, pos, with_labels=True, arrows=True, node_size=70, node_color='lightblue', font_size=10, font_weight='bold')
			plt.axis('equal') 
			# Display the plot
			#plt.title("Directed Graph from Adjacency Matrix")
			plt.show()

deltat=.004

self_int=-.1

for k1 in range(100):
	for EW in np.arange(0,5.0001,.25):
		start = datetime.datetime.now()
		print("making a network with EW = "+str(EW))
		#new_net=network(dt=.000002,N_nodes=400,edge_weigth=2400,self_interaction=-1000,search_path="dhaerarcial2haerarcial",edge_prob=edg_prb)
		new_net=network(dt=deltat,N_nodes=100,edge_weigth=EW,self_interaction=self_int,need_k=True,search_path="dhaerarcial2haerarcial",edge_prob=0)
		#new_net.plot_network()
		#print(np.round(new_net.dynamic_matrix),1)
		print("3-J: ",new_net.J)
		print("4-kreiss constant: ",new_net.K)
		print(max(new_net.eigenvalues))
		print("network is made")
		print("saving: delta s2 p n flambda J ")
		#saving properries
		#saving the high
		print("calculating the std")
		new_net.mean_bubble_height_calc()
		print("all needed data are gathered!")
		print("1-N: ",new_net.N_nodes)
		print("2-P_edges: ",new_net.edge_prob)
		print("3-J: ",new_net.J)
		#print("4-kreiss constant: ",new_net.K)
		print("5-fast lamda: ",new_net.lambda_fast)
		print("6-S2: ",new_net.m2_S)
		print("7-std:",new_net.mean_bubble_height)
		now = datetime.datetime.now()
		print("this datapoint took:",now-start)
		print("another round is completed!")
		print("saving data in ecxel")
		with open("data_full_std10e7_new-v14_test.csv","a") as file:
			file.write(str(np.real(new_net.J))+","+str(EW)+","+str(np.real(new_net.K))+","+str(new_net.m2_S) +","+str(new_net.N_nodes
	) +","+str(new_net.edge_prob)+","+str(np.real(new_net.mean_bubble_height[0]))+"\n")
		print("********************************************","EW = " ,str(EW))
		print("********************************************","j = ",new_net.J)
		print("********************************************","std = ",new_net.mean_bubble_height[0])
		print("*************************************************************************************")
		print("*************************************************************************************")
		

