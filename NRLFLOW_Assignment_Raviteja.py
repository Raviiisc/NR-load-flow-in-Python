import numpy as np

import cmath
import math

# program for reading the input data and copying into busdata array


def bus_data(filename): 
     file = open(filename,"r") #Read the input file
     list_data = list(file.readlines())
     list_data = list_data[1:]

     branch_dline = []
     ref = 0
     for d in list_data:
        if "BUS DATA FOLLOWS" in d:
            ref = 1
            continue

        if ref == 1:
            if "-999" in d:
                break

 
        branch_dline.append(d)
  
  


     line_data = []
     
     bspace = ""
     
     index = 1
     
     for branchData in branch_dline:

        bspace = branchData.strip() # removes extra space
        
        temp_list = bspace.split() #converts to array

        L0 = [index]

        L2 = temp_list[6:8]

        L3 = temp_list[9:13]
        
        L4 = temp_list[14:18]
        
        L5 = [temp_list[19]]
        

        Lsum = L0 + L2 + L3 + L4 + L5  #Combining the lists

       

        line_data_list = list(map(float, Lsum))

        line_data.append(line_data_list) #creating 2D array
  
  
        index = index + 1
        
     return line_data
  
  
data_list = bus_data("yBusData.txt")


busdata = np.array(data_list)

busno=np.array(busdata[:,0])
btype=np.array(busdata[:,1])

finalvoltage=np.array(busdata[:,2])
Pl=np.array(busdata[:,3])
Ql=np.array(busdata[:,4])
Pg=np.array(busdata[:,5])
Qg=np.array(busdata[:,6])
Qmax=np.array(busdata[:,8])
Qmin=np.array(busdata[:,9])
Theta=np.array(busdata[:,11])



#Main program#


# Newton-Raphson Load Flow Analysis in rectangular form

from NRYBUS_Assignment_Raviteja import*   #calling Y bus from NRYBUS_Assignment_Raviteja
                             
busd = busdata;           # Calling busdatas
BMva = 100;                # Base MVA taken as 100
bus = busdata[:,0];        #Bus Number
typeofbus = busd[:,1];     # Type of Bus 3-Slack, 2-PV, 0-PQ

V = busd[:,2];              # Specified Voltage..
          
delta =busd[:,11]; 


        
Pg = busd[:,5]/BMva;        # PGi
Qg = busd[:,6]/BMva;        # QGi
Pl = busd[:,3]/BMva;        # PLi
Ql = busd[:,4]/BMva;        # QLi
Qmin = busd[:,9]/BMva;      # Minimum Reactive Power Limit
Qmax = busd[:,8]/BMva;      #Maximum Reactive Power Limit

P=np.zeros([int(nbus),1])
Q=np.zeros([int(nbus),1])
Psp=np.zeros([int(nbus),1])
Qsp=np.zeros([int(nbus),1])
P = np.subtract(Pg,Pl);     # Pi = PGi - PLi..
Q =np.subtract(Qg ,Ql);     # Qi = QGi - QLi..
Psp = P;                    # P Specified
Qsp = Q;                    #Q specified
                              
G = Y.real;     #conductance

B = Y.imag;    #susceptance



pv =[] #initialize in list form
pq=[]                        
for i in range(0,int(nbus)):
  if np.bitwise_or(typeofbus[i]==3,typeofbus[i]==2):
    pv.append(i) 
    
  else:                      
    pq.append(i)
                         
                          
pvn =len(pv)           # No. of PV buses..
pqn =len(pq)           # No. of PQ buses..



nbus=nbus.astype(int)

Iter = 1;

while(Iter<=3):   # Iteration starting
    
# Only three iterations are taken in order to converge the solution. 
# It is notoced that after third iteration, the values diverge.  
    P =np.zeros([int(nbus),1]);
    Q =np.zeros([int(nbus),1]);
  
    # Calculate P and Q  in rectangular form
    
    for i in range(0,int(nbus)):
       for k in range(0,int(nbus)): 
           
           #Converting into rectangular form
           
            g=((G[i][k])*np.cos(delta[i]-delta[k]) + (B[i][k])*np.sin(delta[i]-delta[k])) # Real part
            
            P[i] = P[i] + (V[i])*(V[k])*g
            
            e=(G[i,k]*np.sin(delta[i]-delta[k]) - (B[i][k])*np.cos(delta[i]-delta[k]))  #Imaginary part
            
            Q[i] = Q[i] + (V[i])*(V[k])*e
    
    
# Check Q-limit violations

    if(np.logical_and(Iter <= 3,Iter >2)):    
      for n in range(1,nbus):
           if (typeofbus[n-1] == 2):
               QG = Q[n]+Ql[n]
               if(QG < Qmin[n]):
                   V[n] = V[n] + 0.01;
               if(QG > Qmax[n]):
                   V[n] = V[n] - 0.01;

    

   #Change from specified value
   
   
    dPa =np.zeros([int(pvn+pqn),1])
    dQa =np.zeros([int(pvn+pqn),1])

    
    for i in range(0,int(nbus)):
        dPa[i][0]=Psp[i]-P[i]
        
    for i in range(0,int(nbus)):
        dQa[i][0]=Qsp[i]-Q[i]
        

                             
   
    k = 0;    
    
    dQ = np.zeros([pqn,1]);
    for i in range(0,int(nbus)): 
        if (typeofbus[i] == 3):
            
            dQ[k][0] = dQa[i][0];
            k = k+1;              
    
    dP=np.zeros([int(pvn+pqn-1),1])
    for i in range(1,int(nbus)):
        dP[i-1][0]=dPa[i][0]
        
    # Mismatch vector
    
    M=np.zeros([int(2*pqn+pvn-1),1])
    M=np.concatenate((dP,dQ),axis=0)
 
    # J1 matrix
    
    J1 = np.zeros([nbus-1,nbus-1]);
    for i in range(0,(nbus)-1):       
        m = i+1;  #edited m=m+1
        for k in range(0,int(nbus)-1):    
            n = k+1;
            if(n == m):
                for n in range(0,int(nbus)):
                    
                    g1=(-G[m][n]*np.sin(delta[m]-delta[n]))
                    
                    e1=(B[m][n]*np.cos(delta[m]-delta[n]))
                    
                    J1[i][k] = J1[i][k] + V[m]* V[n]*(g1+e1)
                
                J1[i][k] = J1[i][k] - V[m]*V[m]*B[m][m]
            else:
               g1=(G[m][n]*np.sin(delta[m]-delta[n]))
               
               e1=(B[m][n]*np.cos(delta[m]-delta[n]))
               
               J1[i][k] = V[m]* V[n]*(g1-e1)
            
    

    
    # J2 matrix
    
    J2 =np.zeros([nbus-1,pqn])
    for i in range(0,int(nbus)-1):
        m = i+1;
        for k in range(0,pqn): 
            n = pq[k];
            if (n == m):
                for n in range(0,nbus):
                   g2=(G[m][n]*np.cos(delta[m]-delta[n]))
                   
                   e2=(B[m][n]*np.sin(delta[m]-delta[n]))
                   
                   J2[i][k] = J2[i][k] + V[n]*(g2+e2)
                
                J2[i][k] = J2[i][k] + V[m]*G[m][m]
            else:
                g2=(G[m][n]*np.cos(delta[m]-delta[n]))
                
                e2=(B[m][n]*np.sin(delta[m]-delta[n]))
                
                J2[i][k] = V[m]*(g2+e2)
            
    
    
    # J3 matrix
    
    J3 =np.zeros([pqn,nbus-1])
    for i in range(0,pqn): 
        m = pq[i];
        for k in range(0,int(nbus)-1): 
            n = k+1;
            if(n == m):
                for n in range(0,int(nbus)-1):
                    
                    g3=(G[m][n]*np.cos(delta[m]-delta[n]))
                    
                    e3=(B[m][n]*np.sin(delta[m]-delta[n]))
                    
                    J3[i][k] = J3[i][k] + V[m]* V[n]*(g3+e3)
                
                J3[i][k] = J3[i][k] - V[m]*V[m]*G[m][m]
                
            else:
                g3=(-G[m][n]*np.cos(delta[m]-delta[n])) 
                
                e3=(B[m][n]*np.sin(delta[m]-delta[n]))
                
                J3[i][k] = V[m]* V[n]*(g3-e3)
                
            
    
    # J4 matrix
    
    J4 =np.zeros([pqn,pqn])
    for i in range(0,pqn): 
        m = pq[i]
        for k in range(0,pqn): 
            n = pq[k]
            if(n == m):
                for n in range(0,nbus):
                    b4=(G[m][n]*np.sin(delta[m]-delta[n]))
                    
                    e4=(B[m][n]*np.cos(delta[m]-delta[n]))
                    
                    J4[i][k] = J4[i][k] + V[n]*(b4-e4)
                    
                J4[i][k] = J4[i][k] - V[m]*B[m][m]
                
            else:
                
                b4=(G[m][n]*np.sin(delta[m]-delta[n])) 
                
                e4=(B[m][n]*np.cos(delta[m]-delta[n]))
                
                J4[i][k] = V[m]*(b4-e4)
            
        
    

    ST=2*int(pqn)+int(pvn)-1
   
    Jt=np.concatenate((J1,J2),axis=1)
    
    Ju=np.concatenate((J3,J4),axis=1)
    
    J=np.concatenate((Jt,Ju),axis=0)       #Jacobian matrix
    
    
    X=np.zeros([int(2*pqn+pvn-1),1])    #Initializing the difference vector
    
  
# Calculating the inverse of the matrix 
  
    
    Jinv=np.zeros([int(2*pqn+pqn-1),int(2*pqn+pvn-1)])
    Jinv=np.linalg.inv(J)
   
 # calculating the difference vector X
   
    for i in range(0,int(2*pqn+pvn-1)):
        for j in range(0,int(2*pqn+pvn-1)):
            X[i]=X[i]+Jinv[i][j]*M[j]
            
    
   #correction vectors    
   
    dTh=np.zeros([int(pvn+pqn-1),1])
    for i in range(0,int(nbus)-1):
        dTh[i][0]=X[i-1][0]
     
    dV=np.zeros([int(pqn),1])
    for i in range(0,int(pqn)):
        dV[i][0]=X[i+int(pqn)][0]

    
    
    # Updating State Vectors

    for i in range(0,nbus-1):
        delta[i]=(np.sum(dTh[i][0],delta[i].astype(int)))
       
    k = 0;
    for i in range(0,int(nbus)):
        if (typeofbus[i] == 3):
                                    
            V[i]=abs(dV[k][0]+V[i])
            
            k = k+1;
    
    
    Iter = Iter + 1;
    
#End of while loop       
        

for i in range(0, int(nbus)):  # Radians to degrees
    delta[i]=np.rad2deg(delta[i])
    
print(V)
print(delta)


# Store the output in a file with name "Nroutput"
   
file1 = open("Nroutput", "w")  
content1 = str(V) 
content2=str(delta)
file1.write(content1)
file1.write("     ")
file1.write(content2) 
file1.close()
    
# End of the program

