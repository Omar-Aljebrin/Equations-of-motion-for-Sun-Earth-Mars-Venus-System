
# =============================================================================
# GLOBAL VAR - GLOBAL MODS
# -----------------------------------------------------------------------------
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt



# -----------------------------------------------------------------------------
# END GLOBAL VAR - GLOBAL MODS
# =============================================================================

# =============================================================================
# FUNCTIONS
# -----------------------------------------------------------------------------

def cromer(t0,tf,h,f,rv0):
    
    #prep
    nt=int((tf-t0)/h)  #3number of points in time 
    
    r=np.zeros((3,nt,3))
    v=np.zeros((3,nt,3))
    
    #intailaizing postions
    
    #earth
    # ip t  xyz
    r[0][0][0]=rv0[0][0]
    r[0][0][1]=rv0[0][1]
    r[0][0][2]=rv0[0][2]
    
    #mars
    r[1][0][0]=rv0[1][0]
    r[1][0][1]=rv0[1][1]
    r[1][0][2]=rv0[1][2]
    
    #venus
    r[2][0][0]=rv0[2][0]
    r[2][0][1]=rv0[2][1]
    r[2][0][2]=rv0[2][2]
 
    
     #intalizizing velocity
     
     #earth
    v[0][0][0]=rv0[0][3]
    v[0][0][1]=rv0[0][4]
    v[0][0][2]=rv0[0][5]
     
     #mars
    v[1][0][0]=rv0[1][3]
    v[1][0][1]=rv0[1][4]
    v[1][0][2]=rv0[1][5]
     
     #venus
    v[2][0][0]=rv0[2][3]
    v[2][0][1]=rv0[2][4]
    v[2][0][2]=rv0[2][5]
                    
    # time
    for k in range(nt-1):
        a=f(r[0,k,:],r[1,k,:],r[2,k,:])
        
        # planets
        for i in range(3):
       
            #moves in compenents 
            for j in range(3):
           
           
               
               
               v[i][k+1][j]=v[i][k][j]+h*a[i][j]
               
               r[i][k+1][j]=r[i][k][j]+h*v[i][k+1][j]
               
    return r,v,nt
       


        
def doon(rse,rsm,rsv):
    #prep for sum
    a=np.zeros((3,3))
    G=6.67e-11  #si value of G
  
    #mass of each stellar object
    Ms=1.9891e30   # (kg) mass of the Sun 
    Me=5.97219e24  # (kg) mass of the Earth
    Mm=6.41710e23  # (kg) mass of the Mars
    Mv=4.86747e24  # (kg) mass of the Venus
    
    M=[[Ms,Mm,Mv],[Ms,Me,Mv],[Ms,Mm,Me]]
    
   
    
    #earth
    rem=np.subtract(rsm,rse)  #rsm-rse
    rev=np.subtract(rsv,rse)  #rsv-rse
   
    #mars
    rmv=np.subtract(rsv,rsm) 
    rme=-1*rem
    
    #3venus
    rve=-1*rev
    rvm=-1*rmv
    
    rv=np.array([[rse,rem,rev],[rsm,rme,rmv],[rsv,rvm,rve]])

 
    
    # finding the magnatuide of rs plnates from rv whihc has all the vectors with theri compenst
    rmag=np.zeros((3,3))    
       
    for i in range(3):
        for k in range(3):
            s=0
            for j in range(3):
                s=s+(rv[i][k][j])**2
            s=np.sqrt(s)
            rmag[i][k]=s
        
    
    #finding a 
    for i in range(3):
        for j in range(3):
            s=0 #maybe wronf for this postions
            for k in range(3):
                #bc how rv is made k is in middle but doesnt matter 
                s=s -((M[i][k]*rv[i][k][j])/(rmag[i][k])**3)
            s=G*s
                
         
            a[i][j]=s

            
            
            
         
            

    return a
            
# -----------------------------------------------------------------------------
# END FUNCTIONS
# =============================================================================

# =============================================================================
# MAIN SCRIPT
# -----------------------------------------------------------------------------
#unit conversions
yr_to_sec=3.15e7
km_to_m=1e3

#importing data
ns=3 #number of sterllar objects 
nc=3 #3number on compents  
#directory name path
dname='../../data/final'   
p = Path(dname)

#extracting multiple files an placing them in a list
fns=list(p.glob('*.csv'))

#room to plave all the data imprted 
data=[None]*ns

#room to place only the intail values of postion and velocity for all 3 objects
rv0=[None]*ns


for i in range(ns):
    data[i]=np.genfromtxt(fns[i],delimiter=',',skip_header=2, usecols=(range(1,7)))

#just taking the first row bc that is where the intail values are 
for i in range(ns):
    rv0[i]=data[i][0]


# moves in the planets data 
for i in range(ns):
    #moves in compents 
    for j in range(nc):
        #foor potion
        rv0[i][j]=(rv0[i][j])*km_to_m   #m
        
        #3for velocity
        rv0[i][j+3]=(rv0[i][j+3])*(km_to_m)   #(m/s)


t0=0 *yr_to_sec  #seec
tf=8 *yr_to_sec #sec 
h=.001 *yr_to_sec #sec 

r,v,nt=cromer(t0,tf,h,doon,rv0)
t=np.linspace(t0,tf,num=nt)

#-------------------------PLOT------------------------------------------------------
#prep for plots for time series plots

lnp=['x(m)','y(m)','z(m)']   #3line name for postion
lnv=['vx(m/s)','vy(m/s)','vz(m/s)']   #3line name for velocity
Tn=['Earth(Euler_Cromer)','Mars(Euler_Cromer)']  #title name
fign=[['E_cr_postions','M_cr_postions'],['E_cr_velocity','M_cr_velocity']]   #figure names

#trojectory plots

#earth and mars x,y
plt.figure(0)
plt.plot(r[0,:,1],r[0,:,0],label='Earth')
plt.plot(r[1,:,1],r[1,:,0], label='Mars')
plt.xlabel('y(m)',fontsize=16)
plt.ylabel('x(m)',fontsize=16)
plt.title('Earth,Mars-Euler_Cromer',fontsize=16)
plt.grid()
plt.legend(fontsize=16)
plt.savefig('EMcrxy_plot.png', dpi=250)

#earth and mars z,x
plt.figure(1)
plt.plot(r[0,:,2],r[0,:,0],label='Earth')
plt.plot(r[1,:,2],r[1,:,0],label='Mars')
plt.xlabel('z(m)',fontsize=16)
plt.ylabel('x(m)',fontsize=16)
plt.title('Earth,Mars-Euler_Cromer',fontsize=16)
plt.grid()
plt.legend(fontsize=16)
plt.savefig('EMcrxz', dpi=250)

#earth and mars z,y
plt.figure(2)
plt.plot(r[0,:,2],r[0,:,1],label='Earth')
plt.plot(r[1,:,2],r[1,:,1],label='Mars')
plt.xlabel('z(m)',fontsize=16)
plt.ylabel('y(m)',fontsize=16)
plt.title('Earth,Mars-Euler_Cromer',fontsize=16)
plt.grid()
plt.legend(fontsize=16)
plt.savefig('EMcryz', dpi=250)


# time series plots

#postion
for i in range(2):
    plt.figure(i+3)
    for j in range(3):
        plt.plot(t,r[i,:,j], label='{}'  .format(lnp[j]))
        plt.title('{}'  .format(Tn[i]) ,fontsize=16)
    plt.xlabel("t(s)",fontsize=16) 
    plt.ylabel('x,y,z (m)',fontsize=16)
    plt.legend(fontsize=16)
    plt.grid()
    plt.savefig('{}_plot.png'.format(fign[0][i]), dpi=250)

#velocity
for i in range(2):
    plt.figure(i+6)
    for j in range(3):
        plt.plot(t,v[i,:,j],label='{}'  .format(lnv[j]))
        plt.title('{}'  .format(Tn[i]) ,fontsize=16)

    plt.xlabel("t(s)",fontsize=16) 
    plt.ylabel('vx,vy,vz (m/s)',fontsize=16)
    plt.legend(fontsize=16)
    plt.grid()
    plt.savefig('{}_plot.png'.format(fign[1][i]), dpi=250)
# -----------------------------------------------------------------------------
# END MAIN SCRIPT
# =============================================================================
