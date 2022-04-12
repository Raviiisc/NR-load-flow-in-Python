import numpy as np
import cmath
import math

# Program to read the data from the input file and copy into the linedata array

def getline_data(filename):

    file = open(filename, "r")
    list_data = list(file.readlines())
    list_data = list_data[2:]
    branch_dline = []
    ref = 0
    for data in list_data:
        if "BRANCH DATA FOLLOWS" in data:
            ref = 1
            continue

        if ref == 1:
            if "-999" in data:
                break

            branch_dline.append(data)



    line_data = []
    b = ""
    index = 1
    for branchData in branch_dline:

        b = branchData.strip() #removes extra space
        temp_list = b.split() #converts to array

        L0 = [index]

        L1 = temp_list[0:2]

        L2 = temp_list[6:9]

        L3 = [temp_list[14]]

        Lsum = L0 + L1 + L2 + L3 # Append the datas


        line_data_list = list(map(float, Lsum))

        line_data.append(line_data_list) # 2D array

        index = index + 1



    return line_data



#Susceptance matrix
Shuntdata=[0,0, 0 ,0, 0, 0, 0, 0, cmath.sqrt(-1) *0.190, 0, 0, 0, 0, 0]


data_list = getline_data("yBusData.txt") # Calling the data for Y bus matrix



linedata = np.array(data_list)

#Ybus program
    
f=np.array(linedata[:,1]) #from bus
t=np.array(linedata[:,2]) #to bus
R=np.array(linedata[:,3]) # Resistance
X=np.array(linedata[:,4]) # Reactance
B=np.array(linedata[:,5]) # shunt susceptance
T=np.array(linedata[:,6]) #Tap ratio
q=max(f)

lineno=np.array(linedata[:,0])
maxim=max(lineno)

for x in range(0,int(maxim)): #Changing the zero tap ratio given in IEEE data sheet to 1
     if T[x]==0:
         T[x]=1


Z=np.array(R+complex(0,1)*X)

y = np.reciprocal(Z);

b=1j*B;


nbus = max(max(f),max(t)); # Number of buses




nbranch = len(f); # no of branch


Y = np.zeros([int(nbus),int(nbus)]) #Initialize Y bus 

Y=Y.astype(complex)


f=f.astype(int) #convert to integer
t=t.astype(int)


for k in range(0,nbranch):
       Y[f[k]-1][t[k]-1]=np.subtract(Y[f[k]-1][t[k]-1],np.divide(y[k],T[k]));
       Y[t[k]-1][f[k]-1]=Y[f[k]-1][t[k]-1];



#formation of diagonal element


for m in range(0,int(nbus)):
    for n in range(0,int(nbranch)):
        if f[n]== m:
           Y[m-1][m-1]= Y[m-1][m-1]+np.divide(y[n],T[k]*T[k])+b[n]   
        if t[n]== m:
           Y[m-1][m-1]= Y[m-1][m-1]+y[n]+b[n];
  

#entry of shuntdata in IEEE n Bus system

    
for i in range(0,int(nbus)):
    if Shuntdata[i]=='':
        Shuntdata[i]=0
        
    Y[i][i]= Y[i][i] + Shuntdata[i];
    

#print(Y) 
# Print the Ybus in a file named "Ybus".

fileYbus = open("Ybus", "w")  
content = str(Y) 
fileYbus.write(content)
fileYbus.write("     ")
fileYbus.close()
    















