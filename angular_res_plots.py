#Program to Plot Angular Resolution as a function of angle for different separation distances, using the results of the monte carlo simulation
import numpy as np
import matplotlib.pyplot as plt

v = 0.994*3e8
time_res_1 = 2.774e-9 #time resolution of scintillator 1, determined in time res section of report
time_res_2 = 2.683e-9
u_t = np.sqrt(time_res_1**2+time_res_1**2) #this is just propagation from time resolution for each scintillator
def uncertainty(theta, v, u_t, u_d, d):
    u_costheta = np.cos(theta)*(np.sqrt((v*(u_t)/(d*np.cos(theta)))**2+(u_d/d)**2))
    u_theta = u_costheta/(np.sin(theta))
    return np.abs(u_theta)

distances = [2.082, 2.865, 3.625, 7.5]
u_distances = [0.205,0.204, 0.205,0.204]

#print(uncertainty(0,v, u_t,u_distances[1],distances[1])) #must exclude this point as gives inf

for i in range(len(distances)):
    x = np.linspace(0.001,np.pi-0.001, 1000)
    print()
    y = uncertainty(x, v, u_t,u_distances[i],distances[i])
    label1 = "d = " + str(distances[i])+ "+/-" + str(u_distances[i]) + "m"
    plt.plot(x*180/np.pi, y*180/np.pi, label = label1)
plt.ylim(0,180)
plt.xlabel('\u03F4 (degrees)')
plt.ylabel('u(\u03F4) (degrees)')
plt.legend()
plt.show()

plt.figure()
for i in range(len(distances)):
    x = np.linspace(0.001,np.pi-0.001, 1000)
    y = ((distances[i]*np.cos(x))/v)*1e9 #Here have multipled by 1e9 to get into nanoseconds
    label1 = "d = " + str(distances[i])+ "+/-" + str(u_distances[i]) + "m"
    plt.plot(x*180/np.pi, y*180/np.pi, label = label1)

plt.xlabel('\u03F4')
plt.ylabel('Expected Time Difference (ns)')
plt.legend()
plt.show()