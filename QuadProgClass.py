class QuadProgClass:
	def __init__(self,orien,pos_v,ang_v,w_t,v_t,epsilon,inc_pos_v,inc_ang_v,er,ep,upper_dq_bounds):
		self.stop=False
		I3 = np.eye(3)
		self.ex = I3[:,0]
		self.ey = I3[:,1]
		self.ez = I3[:,2]
		
		h1 = ez
		h2 = ey
		h3 = ey
		h4 = ex
		h5 = ey
		h6 = ex
		self.P = np.array([[0,0,0], [0.32, 0, 0.78], [0, 0, 1.075], [0, 0, 0.2], [1.142, 0, 0], [0.2, 0, 0], [0,0,0]]).T
		self.q = np.zeros((6, 1))
		self.H = np.array([h1, h2, h3, h4, h5, h6]).T
		self.ttype = np.zeros((1, 6))
		""" """
		self.n = 6
		
		dq_bounds = np.array([[100,110], [90,90], [90,90], [170,190], [120,140], [190,235]]).T
		self.dq_bounds = dq_bounds*np.pi/180
		self.q = np.array([0,0,0,0,np.pi/2,0]).reshape(6, 1)
		self.dq = np.zeros((int(self.n),1))
		dq = np.zeros((int(self.n),1))
		self.R,self.pos = self.__fwdkin(self)
		self.orien=orien
		self.pos_v=pos_v
		self.ang_v=ang_v
		self.w_t=w_t
		self.v_t=v_t
		self.epsilon=epsilon
		self.inc_pos_v=inc_pos_v
		self.inc_ang_v=inc_ang_v
		self.er=er
		self.ep=ep
		self.upper_dq_bounds=upper_dq_bounds
		
	# find closest rotation matrix 
	# A=A*inv(sqrt(A'*A))   
	def Closest_Rotation(R):
		R_n = np.dot(R, inv(sqrtm(np.dot(R.T, R))))
    
		return R_n
		
	__Closest_Rotation=Closest_Rotation

	
	# ROT Rotate along an axis h by q in radius

	def rot(h, q):
		h=h/norm(h)
		R = np.eye(3) + np.sin(q)*self.__hat(h) + (1 - np.cos(q))*np.dot(self.__hat(h), self.__hat(h))
		
		return R
	
	__rot=rot
	
	def hat(h):
		h_hat = np.array([[0, -h[2], h[1]], [h[2], 0, -h[0]], [-h[1], h[0], 0]])
		
		return h_hat
	
	__hat=hat
	
	def fwdkin(self):
		R=np.eye(3)
		p=np.zeros((3,1))
		
		for i in range(self.n):        
			h_i = self.H[0:3,i].reshape(3, 1)
			Ri = np.eye(3)
			
			if self.ttype[0][i] == 0: 
				#rev
				pi = self.P[0:3,i].reshape(3, 1)
				p = p+np.dot(R, pi)
				Ri = self.__rot(h_i,self.q[i])
				R = np.dot(R, Ri)
				R = self.__Closest_Rotation()
			elif self.ttype[i] == 1: 
				#pris
				pi = (self.P[:,i]+self.q[i]*h_i).reshape(3, 1)
				p = p+np.dot(R, pi)
			else: 
				#default pris
				pi = (self.P[:,i]+self.q[i]*h_i).reshape(3, 1)
				p = p+np.dot(R, pi)
	  
		#End Effector T
		p=p+np.dot(R, P[0:3,self.n].reshape(3, 1))
		
		return R, p
	
	__fwdkin=fwdkin
	
	def fwdkin_alljoints(self):
		R=np.eye(3)
		p=np.zeros((3,1))
		RR = np.zeros((3,3,self.n+1))
		pp = np.zeros((3,self.n+1))
		
		for i in range(self.n):
			h_i = self.H[0:3,i]
		   
			if self.ttype[0][i] == 0:
			#rev
				pi = self.P[0:3,i].reshape(3, 1)
				p = p+np.dot(R,pi)
				Ri = self.__rot(h_i,self.q[i])
				R = np.dot(R,Ri)
				R = self.__Closest_Rotation()
			elif self.ttype[i] == 1: 
			#pris
				pi = (self.P[:,i]+self.q[i]*h_i).reshape(3, 1)
				p = p+np.dot(R,pi)
			else: 
			# default pris
				pi = (self.P[:,i]+self.q[i]*h_i).reshape(3, 1)
				p = p+np.dot(R,pi)
			
			pp[:,[i]] = p
			RR[:,:,i] = R
		
		# end effector T
		p=p+np.dot(R, self.P[0:3,self.n].reshape(3, 1))
		pp[:,[n]] = p
		RR[:,:,n] = R
		# parameters for qp
		self.pos=pp[:, -1]
		orien_tmp = Quaternion(matrix=RR[:, :, -1])
		self.orien=np.array([orien_tmp[0], orien_tmp[1], orien_tmp[2], orien_tmp[3]]).reshape(1, 4)
		
		return pp, RR	
	__fwdkin_alljoints=fwdkin_alljoints

	def Joint2Collision(self,Closest_Pt,pp):
		link_dist = []

		for i in range(5):
			link = pp[:,i+1]-pp[:,i]
			link = link/norm(link)
			pp2c = Closest_Pt - pp[:,i]
			
			link_dist.append(norm(pp2c - abs(np.dot(pp2c, link))*link))

		J2C_Joint = link_dist.index(min(link_dist)) + 1
		if(J2C_Joint==1):
			J2C_Joint=2
			
		return J2C_Joint
		
	def getJacobian(self):
		num_joints = len(self.q)

		P_0_i = np.zeros((3,num_joints+1))
		R_0_i = np.zeros((3,3,num_joints+1))


		P_0_i,R_0_i=self.__fwdkin_alljoints()
		
		P_0_T = P_0_i[:,num_joints]

		J = np.zeros((6,num_joints))
		
		for i in range(num_joints):
			if self.ttype[0][i] == 0:
				J[:,i] = np.hstack((np.dot(R_0_i[:,:,i],H[:,i]), np.dot(self.__hat(np.dot(R_0_i[:,:,i], self.H[:,i])), (P_0_T - P_0_i[:,i]))))
		""" """
		
		return J
	# return jacobian of the closest point on panel  
	def getJacobian3(self,Closest_Pt, J2C_Joint):
		num_joints = len(self.q)

		P_0_i = np.zeros((3,num_joints+1))
		R_0_i = np.zeros((3,3,num_joints+1))


		P_0_i,R_0_i=self.__fwdkin_alljoints()
		"""  """
		
		P_0_T = Closest_Pt

		J = np.zeros((6,num_joints))
		
		for i in range(num_joints):
			if self.ttype[0][i] == 0:
				J[:,i] = np.hstack((np.dot(R_0_i[:,:,i],self.H[:,i]), np.dot(self.__hat(np.dot(R_0_i[:,:,i], self.H[:,i])), (P_0_T - P_0_i[:,i]))))
		
		return J

	# return jacobian of the closest point on robot        
	def getJacobian2(self,Closest_Pt,J2C_Joint):

		num_joints = len(self.q)

		P_0_i,R_0_i = self.__fwdkin_alljoints(self)

		P_0_T = P_0_i[:,num_joints]

		J = np.zeros((6,num_joints))

		for i in range(num_joints):
			if self.ttype[0][i] == 0:
				J[:,i] = np.hstack((np.dot(R_0_i[:,:,i], self.H[:,i]), np.dot(self.__hat(np.dot(R_0_i[:,:,i], self.H[:,i])), (P_0_T - P_0_i[:,i]))))

		J[:,J2C_Joint:7] = 0
		link_c = P_0_i[:,J2C_Joint]-P_0_i[:,J2C_Joint-1]
		link_c = link_c/norm(link_c)
		
		P_0_tmp = P_0_i[:,J2C_Joint-1]+ abs(np.dot(Closest_Pt-P_0_i[:,J2C_Joint-1],link_c))*link_c
		
		return J,P_0_tmp
	
	def quat2axang(self):

		s = norm(self.q[0][1:4])
		if s >= 10*np.finfo(np.float32).eps:
			vector = self.q[0][1:4]/s
			theta = 2*np.arctan2(s,self.q[0][0])
		else:
			vector = np.array([0,0,1])
			theta = 0
		axang = np.hstack((vector,theta))
		
		return axang
		
	def getqp_H(self,J, vr):
		n = len(self.dq)
		H1 = np.dot(np.hstack((J,np.zeros((6,2)))).T,np.hstack((J,np.zeros((6,2)))))
		
		tmp = np.vstack((np.hstack((np.hstack((np.zeros((3,n)),vr)),np.zeros((3,1)))),np.hstack((np.hstack((np.zeros((3,n)),np.zeros((3,1)))),self.pos_v)))) 
		H2 = np.dot(tmp.T,tmp)

		H3 = -2*np.dot(np.hstack((J,np.zeros((6,2)))).T, tmp)
		H3 = (H3+H3.T)/2;
		
		tmp2 = np.vstack((np.array([0,0,0,0,0,0,np.sqrt(self.er),0]),np.array([0,0,0,0,0,0,0,np.sqrt(self.ep)])))
		H4 = np.dot(tmp2.T, tmp2)

		H = 2*(H1+H2+H3+H4)

		return H
	
	def getqp_f(self):
		f = -2*np.array([0,0,0,0,0,0,self.er,self.ep]).reshape(8, 1)
    
		return f
		
	# quaternion multiply
	def quatmultiply(q1, q0):
		w0, x0, y0, z0 = q0[0][0], q0[0][1], q0[0][2], q0[0][3]
		w1, x1, y1, z1 = q1[0][0], q1[0][1], q1[0][2], q1[0][3]
		
		return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
						 x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
						 -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
						 x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64).reshape(1, 4)
					  
	__quatmultiply=quatmultiply
	
	def func_xbox(self,button):
		# vx, vy, vz
		if (button[1] == 1):
			self.pos_v = self.pos_v + np.matrix([button[0]*self.inc_pos_v,0,0]).T
		if (button[2] == 1):
			self.pos_v = self.pos_v + np.matrix([0,button[0]*self.inc_pos_v,0]).T
		if (button[3] == 1):
			self.pos_v = self.pos_v + np.matrix([0,0,button[0]*self.inc_pos_v]).T
			
		# wx, wy, wz       
		if (button[4] == 1):
			self.ang_v = self.__quatmultiply(np.array([np.cos(self.inc_ang_v/2), button[0]*np.sin(self.inc_ang_v/2), 0, 0]).reshape(1, 4), self.ang_v)
		if (button[5] == 1):
			self.ang_v = self.__quatmultiply(np.array([np.cos(self.inc_ang_v/2), 0, button[0]*np.sin(self.inc_ang_v/2), 0]).reshape(1, 4), self.ang_v)
		if (button[6] == 1):
			self.ang_v = self.__quatmultiply(np.array([np.cos(self.inc_ang_v/2), 0, 0, button[0]*np.sin(self.inc_ang_v/2)]).reshape(1, 4), self.ang_v)
		