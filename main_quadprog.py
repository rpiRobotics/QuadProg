# -*- coding: utf-8 -*-

import numpy as np
from numpy.linalg import inv
from scipy.linalg import logm, norm, sqrtm
import pygame
from ControlParams import *
from OpenRAVE_rr_server import *
import rpi_abb_irc5
import time
import timeit
from pyquaternion import Quaternion
import quadprog


""" """

def inequality_bound(h,c,eta,epsilon,e):
    sigma = np.zeros((h.shape))
    h2 = h - eta
    sigma[np.array(h2 >= epsilon)] = -np.tan(c*np.pi/2)
    sigma[np.array(h2 >= 0) & np.array(h2 < epsilon)] = -np.tan(c*np.pi/2/epsilon*h2[np.array(h2 >= 0) & np.array(h2 < epsilon)])
    sigma[np.array(h >= 0) & np.array(h2 < 0)] = -e*h2[np.array(h >= 0) & np.array(h2 < 0)]/eta
    sigma[np.array(h < 0)] = e
    
    return sigma

   

def main():
    #Init the joystick
    pygame.init()
    pygame.joystick.init()

    joy=pygame.joystick.Joystick(0)
    joy.init()
    clock=pygame.time.Clock()

    egm=rpi_abb_irc5.EGM()

    OpenRAVE_obj = OpenRAVEObject()

    # Initialize Robot Parameters    
    #ex,ey,ez,n,P,q_ver,H,ttype,dq_bounds = robotParams()
	
    # Initialize Control Parameters
    # initial joint angles
    
    """ """
    #q = np.zeros((6, 1))
    q = np.array([0,0,0,0,np.pi/2,0]).reshape(6, 1)

    
    orien = Quaternion(matrix=R)
    orien = np.array([orien[0], orien[1], orien[2], orien[3]]).reshape(1, 4)

    pos_v = np.zeros((3, 1))
    ang_v = np.array([1,0,0,0])
    
    # joint limits
    lower_limit = np.transpose(np.array([-170*np.pi/180, -65*np.pi/180, -np.pi, -300*np.pi/180, -120*np.pi/180, -2*np.pi]))
    upper_limit = np.transpose(np.array([170*np.pi/180, 85*np.pi/180, 70*np.pi/180, 300*np.pi/180, 120*np.pi/180, 2*np.pi]))
	
    # inequality constraints
    h = np.zeros((15, 1))
    sigma = np.zeros((13, 1))
    dhdq = np.vstack((np.hstack((np.eye(6), np.zeros((6, 1)), np.zeros((6, 1)))), np.hstack((-np.eye(6), np.zeros((6, 1)), np.zeros((6, 1)))), np.zeros((1, 8))))

    # velocities
    w_t = np.zeros((3, 1))
    v_t = np.zeros((3, 1))
    
    # keyboard controls
    # define position and angle step
    inc_pos_v = 0.01 # m/s
    inc_ang_v = 0.5*np.pi/180 # rad/s

    # optimization params
    er = 0.05
    ep = 0.05
    epsilon = 0 # legacy param for newton iters
    
    # parameters for inequality constraints
    c = 0.5
    eta = 0.1
    epsilon_in = 0.15
    E = 0.005
    
    Ke = 1
    
    # create a handle of these parameters for interactive modifications
    quadprogobj = QuadProgClass(orien,pos_v,ang_v.reshape(1, 4),w_t,v_t,epsilon,inc_pos_v,inc_ang_v,er,ep,0)

    dt = 0
	tic=0
    
    while not quadprogobj.stop:

        # Loop reading the joysticks
        for event in pygame.event.get():
            pass

        #used to get time increment it takes for loop to execute
		toc = timeit.default_timer()
		dt = toc - tic
		print dt
        tic = timeit.default_timer()

        
		quadprogobj.q = quadprogobj.q + quadprogobj.dq*dt*0.1
		res, state = egm.receive_from_robot(0.01)

		if res:
			a = np.array(state.joint_angles)
			a = a * 180 / np.pi
			print "Joints: " + str(a)
			egm.send_to_robot([float(x)*180/np.pi for x in quadprogobj.q])
			print "Target Joints: " + str([float(x)*180/np.pi for x in quadprogobj.q])

		pp,RR = quadprogobj.fwdkin_alljoints()
		
		
		stop, Closest_Pt, Closest_Pt_env = OpenRAVE_obj.CollisionReport(quadprogobj.q[0],quadprogobj.q[1],quadprogobj.q[2],quadprogobj.q[3],quadprogobj.q[4],quadprogobj.q[5])
		
		# check self-collision
		if (stop):
			print 'robot is about to self-collide, robot stopped.'
			quadprogobj.pos_v = np.array([0,0,0]).reshape(3, 1)
			quadprogobj.ang_v = np.array([1,0,0,0]).reshape(1, 4)
		
		J2C_Joint = quadprogobj.Joint2Collision(Closest_Pt, pp)
		
		J_eef = quadprogobj.getJacobian()
		
		v_tmp = Closest_Pt-quadprogobj.pos
		
		v_tmp2 = (pp[:, -1] - pp[:, -3]) 
		p_norm2 = norm(v_tmp2)
		v_tmp2 = v_tmp2/p_norm2
		
		# determine if the closest point is on the panel
		if (norm(v_tmp) < 1.5 and np.arccos(np.inner(v_tmp, v_tmp2)/norm(v_tmp))*180/np.pi < 95):
			print '---the closest point is on the panel---'
			J2C_Joint = 6
			J = quadprogobj.getJacobian3(Closest_Pt,J2C_Joint)
			#J,p_0_tmp = getJacobian2(obj.params['controls']['q'], obj.params['defi']['ttype'], obj.params['defi']['H'], obj.params['defi']['P'], obj.params['defi']['n'],Closest_Pt,J2C_Joint)
			
		#if (J2C_Joint < 4):
		else:
			J,p_0_tmp = quadprogobj.getJacobian2(Closest_Pt,J2C_Joint)
		
		#else:            
		 #   J = getJacobian3(obj.params['controls']['q'], obj.params['defi']['ttype'], obj.params['defi']['H'], obj.params['defi']['P'], obj.params['defi']['n'], Closest_Pt,J2C_Joint)
		
		x = joy.get_axis(0)
		if (abs(x) < .2):
			x = 0
		else:
			x = (abs(x) - .2) / .8 * cmp(x, 0)

		# control of linear velocity
		b1 = joy.get_button(0)
		b2 = joy.get_button(1)
		b3 = joy.get_button(2)

		# control of angular velocity
		b4 = joy.get_button(3)
		b5 = joy.get_button(4)
		b6 = joy.get_button(5)
	
		# emergency stop
		b9 = joy.get_button(8)
	
		if (b9 == 1):
			print 'robot stopped'
			quadprogobj.pos_v = np.array([0,0,0]).reshape(3, 1)
			quadprogobj.ang_v = np.array([1,0,0,0]).reshape(1, 4)
  
		button = [x, b1, b2, b3, b4, b5, b6]
		quadprogobj.func_xbox(button)
		
		# update joint velocities
		axang = quadprogobj.quat2axang()

		# desired rotational velocity
		vr = axang[3]*axang[0:3]
		
		# desired linear velocity
		V_desired = quadprogobj.pos_v
					
		Q = quadprogobj.getqp_H(J_eef, vr.reshape(3, 1)) 
		
		# make sure Q is symmetric
		Q = 0.5*(Q + Q.T)
		
		f = quadprogobj.getqp_f()
		f = f.reshape((8, ))

		#Q = matrix(Q, tc='d')
		#f = matrix(f, tc='d')
		
		# bounds for qp
		if quadprogobj.upper_dq_bounds:
			bound = quadprogobj.dq_bounds[1, :]
		else:
			bound = quadprogobj.dq_bounds[0, :]

		LB = np.vstack((-0.1*bound.reshape(6, 1),0,0))
		UB = np.vstack((0.1*bound.reshape(6, 1),1,1))
		LB = matrix(LB, tc = 'd')
		UB = matrix(UB, tc = 'd')
				
		# inequality constrains A and b
		h[0:6] = quadprogobj.q - lower_limit.reshape(6, 1)
		h[6:12] = upper_limit.reshape(6, 1) - quadprogobj.q
		
		dx = Closest_Pt_env[0] - Closest_Pt[0]
		dy = Closest_Pt_env[1] - Closest_Pt[1]
		dz = Closest_Pt_env[2] - Closest_Pt[2]
		
		""" """
		dist = np.sqrt(dx**2 + dy**2 + dz**2)
		
		# derivative of dist w.r.t time
		der = np.array([dx*(dx**2 + dy**2 + dz**2)**(-0.5), dy*(dx**2 + dy**2 + dz**2)**(-0.5), dz*(dx**2 + dy**2 + dz**2)**(-0.5)])

		""" """
		h[12] = dist - 0.05
		""" """ """ """
		#dhdq[12, 0:6] = np.dot(-der.reshape(1, 3), J_eef2[3:6,:])
		dhdq[12, 0:6] = np.dot(-der[None, :], J[3:6,:])
		
		sigma[0:12] =inequality_bound(h[0:12], c, eta, epsilon_in, E)
		sigma[12] = inequality_bound(h[12], c, eta, epsilon_in, E)           
		
		A = dhdq
		b = sigma
		
		A = np.vstack((A, np.eye(8), -np.eye(8)))
		b = np.vstack((b, LB, -UB))
		b = b.reshape((29, ))

		# solve the quadprog problem
		dq_sln = quadprog.solve_qp(Q, -f, A.T, b)[0]
	   
		b7 = joy.get_button(6)
		b8 = joy.get_button(7)
				   
		# set desired eef position and orientation
		if (b7 == 1 and b8 == 1):
			x_des = pp[:, -1]
			R_des = RR[:,:,-1]
			print 'desired position set'
			print 'desired position set'
		
		# it seems that the equality constraints works better when close to the obstacle
		# equality constraints for maintaining end-effector position (pure rotation)    
		if (b7 == 1):
			print 'pure rotational movement'
			A_eq_pos = np.hstack((J_eef[3:6,:], np.zeros((3, 2))))
			b_eq_pos = 0.05*Ke*(x_des - pp[:, -1])
			
			# stack equality constrains on top of the inequality constraints
			A = np.vstack((A_eq_pos, A))
			b = np.concatenate((b_eq_pos.reshape((3, 1)), b.reshape((29, 1))), axis=0)
			b = b.reshape((32, ))
			
			# the last argument specify the number of equality constraints
			dq_sln = quadprog.solve_qp(Q, -f, A.T, b, A_eq_pos.shape[0])[0]
			
			A = np.delete(A, [0, 1, 2], axis=0)
			b = np.delete(b, [0, 1, 2])
			
		# equality constraints for maintaining end-effector orientation (pure translation)
		if (b8 == 1):
			print 'pure translational movement'
			A_eq = np.hstack((J_eef[0:3,:], np.zeros((3, 2))))            
			w_skew = logm(np.dot(RR[:,:,-1],R_des.T))
			w = np.array([w_skew[2, 1], w_skew[0, 2], w_skew[1, 0]])
			b_eq = -0.05*Ke*w
			
			# stack equality constrains on top of the inequality constraints
			A = np.vstack((A_eq, A))
			b = np.concatenate((b_eq.reshape((3, 1)), b.reshape((29, 1))), axis=0)
			b = b.reshape((32, ))
			
			# the last argument specify the number of equality constraints
			dq_sln = quadprog.solve_qp(Q, -f, A.T, b, A_eq.shape[0])[0]
			
			A = np.delete(A, [0, 1, 2], axis=0)
			b = np.delete(b, [0, 1, 2])
			
		if len(dq_sln) < quadprogobj.n:
			quadprogobj.dq = np.zeros((6,1))
			V_scaled = 0
			print 'No Solution'
		else:
			quadprogobj.dq = dq_sln[0: int(quadprogobj.n)]
			quadprogobj.dq = quadprogobj.dq.reshape((6, 1))
			V_scaled = dq_sln[-1]*V_desired
			vr_scaled = dq_sln[-2]*vr.reshape(3,1)
		
		V_linear = np.dot(J_eef[3:6,:], quadprogobj.dq)
		V_rot = np.dot(J_eef[0:3,:], quadprogobj.dq)
					   
		print '------Scaled desired linear velocity------'
		print V_scaled
				   
		print '------Real linear velocity by solving quadratic programming------'
		print V_linear
		
		print '------Scaled desired angular velocity------'
		print vr_scaled
		
		print '------Real angular velocity by solving quadratic programming------'
		print V_rot

    pygame.quit()



        
if __name__ == '__main__':
    main()
    
