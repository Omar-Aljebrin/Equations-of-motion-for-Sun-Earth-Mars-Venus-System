
# =============================================================================
# GLOBAL VAR - GLOBAL MODS
# -----------------------------------------------------------------------------
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# END GLOBAL VAR - GLOBAL MODS
# =============================================================================

# =============================================================================
# FUNCTIONS
# -----------------------------------------------------------------------------

def dUdt(t,U,G,Ms,Me,Mm,Mv):
    no=3  #number of objects
    nc=3  #number of compents
    a=np.zeros((no,nc))

    #input for loop
   
  
    
    M=[[Ms,Mm,Mv],[Ms,Me,Mv],[Ms,Mm,Me]]
    
    #nee work respect to U remmber how you definded U0 it was not like that 
    #inputting U values
    
    #U's postion need reworks 
    rse=np.array([U[0],U[1],U[2]])  #from e to s
    rsm=np.array([U[6],U[7],U[8]])  #3from v to s
    rsv=np.array([U[12],U[13],U[14]])
    
    #velocity
    
    #earth
    vxe=U[3]
    vye=U[4]
    vze=U[5]
    #mars
    vxm=U[9]
    vym=U[10]
    vzm=U[11]
    #venus
    vxv=U[15]
    vyv=U[16]
    vzv=U[17]
    
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
       
    for i in range(no):
        for k in range(no):
            s=0
            for j in range(nc):
                s=s+(rv[i][k][j])**2
            s=np.sqrt(s)
            rmag[i][k]=s
        
    
    #finding a 
    for i in range(no):
        for j in range(nc):
            s=0 
            for k in range(no):
                #bc how rv is made k is in middle but doesnt matter 
                s=s -((M[i][k]*rv[i][k][j])/(rmag[i][k])**3)
            s=G*s
                
         
            a[i][j]=s
            
#-------------------Done finding a--------------------------------------------
    
#constructing dUdt
    dUdt=[
        #earth
        vxe,vye,vze,a[0][0],a[0][1],a[0][2],
        #mars
        vxm,vym,vzm,a[1][0],a[1][1],a[1][2],
        #venus
        vxv,vyv,vzv,a[2][0],a[2][1],a[2][2]
          ]
    
    
    return np.array(dUdt)
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
    #using 'usecol to not inculde time col
    
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

#---------------------------Done importing data---------------------------------

#other given values

t0=0 *yr_to_sec  #seec
tf=12 *yr_to_sec #sec 
h=.001 *yr_to_sec #sec 

nt=int((tf-t0)/h)+1 #number of points for t bc ivp needs to know homany m[points]
#prep for solve_ivp
tspan=[t0,tf]
t=np.linspace(t0,tf,num=nt)
U0 = np.concatenate(rv0)
#rs0 is intail condestion but they are blocked not in one elist in serirs


#the agrs for the funcations 

G=6.67e-11  #si value of G
  
#mass of each stellar object 


Ms=1.9891e30   # (kg) mass of the Sun 
Me=5.97219e24  # (kg) mass of the Earth
Mm=6.41710e23  # (kg) mass of the Mars
Mv=4.86747e24  # (kg) mass of the Venus




#using ivp to appromxate

sol = solve_ivp(dUdt,tspan,U0,
                method='Radau',
                t_eval=t,  # times at which to store 
                            # the computed solution, 
                            # must be sorted and 
                            # lie within t_span.
                args=(G,Ms,Me,Mm,Mv))  #maybe orr of arg matter?

# Extracting ivp values
t= sol.t

#earth

# postion

#earth
rse=[sol.y[0,:],sol.y[1,:],sol.y[2,:]]
#mars
rsm=[sol.y[6,:],sol.y[7,:],sol.y[8,:]]
#venus
rsv=[sol.y[12,:],sol.y[13,:],sol.y[14,:]]

r=[rse,rsm,rsv]
#velocity

#earth
ve=[sol.y[3,:],sol.y[4,:],sol.y[5,:]]
#mars
vm=[sol.y[9,:],sol.y[10,:],sol.y[11,:]]
#venus
vv=[sol.y[15,:],sol.y[16,:],sol.y[17,:]]

v=[ve,vm,vv]

#-------------------------PLOT------------------------------------------------------
#prep for plots for time series plots

lnp=['x(m)','y(m)','z(m)']   #3line name for postion
lnv=['vx(m/s)','vy(m/s)','vz(m/s)']   #3line name for velocity
Tn=['Earth(ivp)','Mars(ivp)']  #title name
fign=[['E_ivp_postions','M_ivp_postions'],['E_ivp_velocity','M_ivp_velocity']]   #figure names

#trojectory plots

#earth and mars x,y
plt.figure(0)
plt.plot(r[0][1],r[0][0],label='Earth')
plt.plot(r[1][1],r[1][0], label='Mars')
plt.xlabel('y(m)',fontsize=16)
plt.ylabel('x(m)',fontsize=16)
plt.title('Earth,Mars-ivp',fontsize=16)
plt.grid()
plt.legend(fontsize=16)
plt.savefig('EMivpxy_plot.png', dpi=250)

#earth and mars z,x
plt.figure(1)
plt.plot(r[0][2],r[0][0],label='Earth')
plt.plot(r[1][2],r[1][0],label='Mars')
plt.xlabel('z(m)',fontsize=16)
plt.ylabel('x(m)',fontsize=16)
plt.title('Earth,Mars-ivp',fontsize=16)
plt.grid()
plt.legend(fontsize=16)
plt.savefig('EMivpxz', dpi=250)

#earth and mars z,y
plt.figure(2)
plt.plot(r[0][2],r[0][1],label='Earth')
plt.plot(r[1][2],r[1][1],label='Mars')
plt.xlabel('z(m)',fontsize=16)
plt.ylabel('y(m)',fontsize=16)
plt.title('Earth,Mars-ivp',fontsize=16)
plt.grid()
plt.legend(fontsize=16)
plt.savefig('EMvvyz', dpi=250)


# time series plots

#postion
for i in range(2):
    plt.figure(i+3)
    for j in range(3):
        plt.plot(t,r[i][j], label='{}'  .format(lnp[j]))
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
        plt.plot(t,v[i][j],label='{}'  .format(lnv[j]))
        plt.title('{}'  .format(Tn[i]) ,fontsize=16)

    plt.xlabel("t(s)",fontsize=16) 
    plt.ylabel('vx,vy,vz (m/s)',fontsize=16)
    plt.legend(fontsize=16)
    plt.grid()
    plt.savefig('{}_plot.png'.format(fign[1][i]), dpi=250)


#plotting Earth Venus dance
plt.figure(20)
ax = plt.gca()
plt.plot([rse[0],rsv[0]],[rse[1],rsv[1]],
              alpha=0.005)
             
ax.set_aspect(aspect=1)
plt.title('Earth-Venus Dance')
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.grid()
plt.savefig('Earth-Venus-Dance_plot.png', dpi=250)
plt.show()

# -----------------------------------------------------------------------------
# END MAIN SCRIPT5
# =============================================================================
