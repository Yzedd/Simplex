# ======================================================================== #
#                                                                          #
#       Revised Simplex and Revised Dual Simplex Method (Two-Phase)        #
#                            Zhan Yang, 2019                               #
#                                                                          #
# ======================================================================== #

import numpy as np

# =================================== #
# Initialize the optimization problem #
# =================================== #

"""
design for
min z(x) = ct*x
s.t
Ax <= b
x >= 0
"""
#"""
opti = input("Make sure it's a minimization problem => ")

c = input("Objective coefficients => ")
c = np.array([float(n) for n in c.split()])

b = input("Constraints => ")
b = np.array([float(n) for n in b.split()])

equ = input("Signs (only '<' or '=') => ")
equ = np.array([str(n) for n in equ.split()])

A = input("Constraint coefficients => ")
A = [float(n) for n in A.split()]
A = np.array(A).reshape((len(c), len(b)), order='F')
print(A.T)
#"""

# ============================= #
# Analyze the Constraint Matrix #
# ============================= #

def check(A,b,c):
    
    a2 = np.zeros((1, len(A[0,:])), dtype=np.float64)
    a3 = np.zeros((1, len(A[0,:])), dtype=np.float64)
    B = np.zeros((len(b), len(b)), dtype=np.float64)
    for i in range(len(B[0,:])):
        B[i,i] = 1
    
    # Basic variables
    x = list(np.zeros(len(b), dtype=np.int))
    
    for i in range(len(A[:,0])):
        for j in range(len(A[0,:])):
            
            if A[i,j] == 1:
                a2[0,:] = A[i,:]
                a2[0,j] = 0                
                if (a2 == a3).all():
                    x[j] = i+1
    
    for i in range (len(x)):
        
        if x[i] == 0:
            A = np.vstack((A, B[:,i]))
            x[i] = len(A[:,0])
    
    A = np.transpose(A)
    d = np.zeros(len(A[0,:]) - len(c))
    c = np.hstack((c,d))
    B,D,a = basic(A,x)
    
    return A,B,D,c,x

# ============================= #
# Construct Basic and Non-basic #
# ============================= #
    
def basic(A,x):
    
    B = np.zeros((len(A[:,0]), len(A[:,0])), dtype=np.float64)    
    D = np.zeros((len(A[:,0]), len(A[0,:])-len(A[:,0])), dtype=np.float64)
    a = list(np.zeros(len(A[0,:]), dtype=np.int))

    for i in range (len(a)):
        
        a[i] = i+1
    
    for i in range (len(x)):
        
        a.remove(x[i])
        B[:,i] = A[:,x[i]-1]
    
    for i in range(len(a)):
                  
        D[:,i] = A[:,a[i]-1]

    return B,D,a

# ========================== #
# Check Relative Cost Vector #
# ========================== #

def cost_vector(B,D,c,x,a):
    
    cb = np.zeros(len(B[:,0]))
    cd = np.zeros(len(D[0,:]))
    
    for i in range (len(x)):
        
        cb[i] = c[x[i]-1]
    
    for i in range (len(a)):
        
        cd[i] = c[a[i]-1]
         
    else:
        pi = np.dot(cb, np.linalg.inv(B))
        cd_ = cd - pi.dot(D)            
    
    return cd_

# ========= #
# Phase One #
# ========= #

def phase1(A,b,c,equ):
       
    # Artificial variables
    ar = list(np.zeros(len(equ), dtype=int))
    A1 = np.vstack((A.T, c)).T
    equ = np.hstack((equ, '='))
    b1 = np.hstack((b, 0))
    cr = np.zeros(len(c), dtype=np.float64)

    B = np.zeros((len(equ), len(equ)), dtype=np.float64)
    for i in range(len(B[0,:])):
        B[i,i] = 1
        
    for i in range (len(equ)):
        
        if equ[i] == '=':
            a1 = np.delete(A1, i, axis=1)
            
            w = 0
            for j in range (len(a1[:,0])):
                if (a1[j,:] == 0).all():
                    w = 1
            if w == 0:
                A1 = np.vstack((A1, B[:,i]))
                cr = np.hstack((cr, 1))
    
    cr[len(cr)-1] = 0
    j = 0
    for i in range (len(cr)):
        if cr[i] == 1:
            ar[j] = i+1
            j += 1                   
    
    A1,B,D,cr,x = check(A1,b1,cr)
    B,D,a = basic(A1,x)
    cd_ = cost_vector(B,D,cr,x,a)
            
    iteration = 0

    if np.linalg.det(B) == 0:
        print("No feasible solution")
        
    else:
        b_ = np.dot(np.linalg.inv(B), b1)

        if min(cd_) >= 0:
    
            x = x
            np.asarray
        else:
            while list(set(x) & set(ar)) != []:
        
                for i in range (len(cd_)):
                    
                    if cd_[i] == min(cd_):
                        if cd_[i] < 0:
                            enter = a[i]
                
                a_ = np.dot(np.linalg.inv(B), A1[:,enter-1])
            
                if max(a_) <= 0:
                    print("The problem is unbounded.")
                    x = [0]
                    break
                    
                if iteration >= 50:
                    x = [0]
                    print("Cycling")
                    break
        
                else:
                    iteration = iteration + 1
                    print("Phase One Iteration", iteration)
                    
                    nz = np.zeros(len(a_), dtype = int)
                    b_ = np.delete(b_, len(b_)-1)
                    a_ = np.delete(a_, len(a_)-1)
                    for i in range (len(a_)):
                        if a_[i] > 0:
                            nz[i] = i + 1
                            nz2 = np.nonzero(nz)
                            nz2 = nz2[0]
                            a2 = np.zeros(len(nz2))
                            b2 = np.zeros(len(nz2))
                                    
                    for i in range (len(a2)):
                                        
                        a2[i] = a_[nz[nz2[i]] - 1]
                        b2[i] = b_[nz[nz2[i]] - 1]
                        
                    for i in range (len(a_)):
                                        
                        if a_[i] > 0:
                            if b_[i]/a_[i] == min(b2/a2):
                                leave = i
                        
                    x[leave] = enter                    
                    B,D,a = basic(A1,x)
                    b_ = np.dot(np.linalg.inv(B), b1)
                    cd_ = cost_vector(B,D,cr,x,a)
    
    A = np.transpose(A)
    x.remove(x[len(x)-1])
    B,D,a = basic(A,x)
    
    return A,B,D,c,x

# ====================== #
# Revised Simplex Method #
# ====================== #
    
def revised_simplex(a_,b_,cd_,a,x):

    for i in range (len(cd_)):
                    
        if cd_[i] == min(cd_):
            if cd_[i] < 0:
                enter = a[i]    
                
    nz = np.zeros(len(a_), dtype = int)
    for i in range (len(a_)):
        if a_[i] > 0:
            nz[i] = i + 1
            nz2 = np.nonzero(nz)
            nz2 = nz2[0]
            a2 = np.zeros(len(nz2))
            b2 = np.zeros(len(nz2))
                    
    for i in range (len(a2)):
                        
        a2[i] = a_[nz[nz2[i]] - 1]
        b2[i] = b_[nz[nz2[i]] - 1]
        
    for i in range (len(a_)):
                        
        if a_[i] > 0:
            if b_[i]/a_[i] == min(b2/a2):
                leave = i
        
    x[leave] = enter
        
    return x

# =========================== #
# Revised Dual Simplex Method #
# =========================== #
    
def revised_dual(B,D,b,c,a,x):
    
    b_ = np.dot(np.linalg.inv(B), b)
    iteration = 0
    
    while min(b_) < 0:
        iteration = iteration + 1
        print("Iteration", iteration)
        
        bp = np.dot(np.linalg.inv(B), D)
        bp1 = np.zeros(len(bp[0]))
        
        if (bp >= 0).all():
            print("No feasible solution")
            x = [0]
            break

        else:
            for i in range (len(b_)):
                        
                if b_[i] == min(b_):
                    if b_[i] < 0:
                        bp1 = bp[i]
            
            cb = np.zeros(len(B[:,0]))
            cd = np.zeros(len(D[0,:]))
        
            for i in range (len(x)):
            
                cb[i] = c[x[i]-1]
        
            for i in range (len(a)):
            
                cd[i] = c[a[i]-1]
                
            zc = np.dot(cb.dot(np.linalg.inv(B)), D) - cd
            
            ng = np.zeros(len(bp1), dtype = int)
            for i in range (len(bp1)):
                if bp1[i] < 0:
                    ng[i] = i + 1
                    ng2 = np.nonzero(ng)
                    ng2 = ng2[0]
                    zc2 = np.zeros(len(ng2))
                    bp2 = np.zeros(len(ng2))
    
            for i in range (len(bp2)):
                            
                zc2[i] = zc[ng[ng2[i]] - 1]
                bp2[i] = bp1[ng[ng2[i]] - 1]
            
            for i in range (len(bp1)):
                            
                if bp1[i] < 0:
                    if abs(zc[i]/bp1[i]) == min(abs(zc2/bp2)):
                        enter = a[i]
            
            for i in range (len(b_)):
                        
                if b_[i] == min(b_):
                    if b_[i] < 0:
                        leave = i
            
            x[leave] = enter
            if len(x) > len(set(x)):
                print("No feasible solution.")
                x = [0]
                break
            B,D,a = basic(A,x)
            b_ = np.dot(np.linalg.inv(B), b)
        
    return B,D,a,x

# =============================== #
# Solver for Optimization Problem #
# =============================== #
    
def solver(A,B,D,a,b,x):
    
    iteration = 0

    if np.linalg.det(B) == 0:
        print("No feasible solution")
        
    else:
        b_ = np.dot(np.linalg.inv(B), b)
                  
        if min(b_) < 0:
                    
            # =========================== #
            # Revised Dual Simplex Method #
            # =========================== #
                        
            B,D,a,x = revised_dual(B,D,b,c,a,x)
                    
        else:
            
            # ====================== #
            # Revised Simplex Method #
            # ====================== #
            cd_ = cost_vector(B,D,c,x,a)
            
            if min(cd_) >= 0:
    
                x = x
            
            else:
                while min(cd_) < 0:
        
                    for i in range (len(cd_)):
                    
                        if cd_[i] == min(cd_):
                            if cd_[i] < 0:
                                enter = a[i]
                
                    a_ = np.dot(np.linalg.inv(B), A[:,enter-1])
            
                    if max(a_) <= 0:
                        print("The problem is unbounded.")
                        x = [0]
                        break
                    
                    if iteration >= 50:
                         x = [0]
                         print("Cycling")
                         break
        
                    else:
                        iteration = iteration + 1
                        print("Iteration", iteration)
                  
                        x = revised_simplex(a_,b_,cd_,a,x)
                        B,D,a = basic(A,x)
                        b_ = np.dot(np.linalg.inv(B), b)
                        cd_ = cost_vector(B,D,c,x,a)

    return B,D,a,x
    
# ================ #
# Optimal Solution #
# ================ #

def optimal(x,c,B,b):
    
    if x == [0]:
        opti_x = []
        z = []
        
    else:    
        opti_x = np.zeros(len(c))
        cb = np.zeros(len(x))
        cv = np.dot(np.linalg.inv(B), b)
        
        for i in range (len(x)):
            
            opti_x[x[i]-1] = cv[i]
            cb[i] = c[x[i]-1]
            
        z = np.dot(cb, cv)
    
    return opti_x,z    

# ============================== #
# Solve the Optimization Problem #
# ============================== #
    
"""
# Test
A = np.transpose(np.asarray([[5,-4,13,-2,1],[1,-1,5,-1,1]], dtype=np.float64))
c = np.asarray([3,-1,-7,3,1], dtype=np.float64)
b = np.asarray([20,8], dtype=np.float64)
equ = np.asarray(['=','='], dtype=np.str)
"""

if (equ == '=').any():
    A,B,D,c,x = phase1(A,b,c,equ)
else:
    A,B,D,c,x = check(A,b,c)

B,D,a = basic(A,x)

B,D,a,x = solver(A,B,D,a,b,x)
opti_x,z = optimal(x,c,B,b)
    
print("Optimal solution: x =", opti_x, ", z =", z)