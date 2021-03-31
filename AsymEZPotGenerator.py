import numpy as np
import scipy.misc as scpy
#import sympy
import math
import matplotlib.pyplot as plt
import collections

# currently works for CB (carbon beta) or CG (carbon gamma)
term_atoms = "CG"
print("Generating bin values for %s"%term_atoms)

def sigmoidal(z, deltaE0, zmid, n, asym = False):
    if asym == False:
        if z > 0:
            z = -z
    return( deltaE0/( 1 + (z/zmid)**n ),( ((-1*deltaE0*n)*(z/zmid)**n) / ( z* (((z/zmid)**n + 1)**2) ) ) )

def gaussian(z, deltaEmin, zmin, sigma2):
    return( deltaEmin* np.exp( -((z-zmin)**2)/(2*sigma2)), -deltaEmin*(2*z - 2*zmin)*np.exp(-(z - zmin)**2/(2*sigma2))/(2*sigma2))

params = collections.OrderedDict()
with open("2012SamishData/Samish%sparams.txt"%(term_atoms), "r") as inData:
    for line in inData:
        line = line.strip().split()
        params[line[0]] = line[1:]


z_coords = np.arange(-30,31,1)
count = 0
with open("AsymEZ_%s.txt"%(term_atoms), "w+") as outData, open("AsymEZ_%s_Derivs.txt"%(term_atoms), "w+") as derivs:
    for key in params.keys():
        print(key)
        count += 1
        Paq = float(params[key][0])
        deltaEz = np.zeros(len(z_coords))
        derivEz = np.zeros(len(z_coords))
        if key[-1] == "1":
            E0 = float(params[key][7])
            Zmid = float(params[key][8])
            N = float(params[key][9])
            for z in range(len(z_coords)):
                result = sigmoidal(z_coords[z], E0, Zmid, N)
                deltaEz[z] = result[0]
                derivEz[z] = result[1]
        elif key[-1] == "2":
            E0 = float(params[key][7])
            Zmid = float(params[key][8])
            N = float(params[key][9])
            for z in range(len(z_coords)):
                if z_coords[z] > float(params[key][12]) and z_coords[z] < float(params[key][13]):
                    result = sigmoidal(z_coords[z], E0, Zmid, N)
                    deltaEz[z] = result[0]
                    derivEz[z] = result[1]
                elif z_coords[z] <= float(params[key][12]):
                    A = float(params[key][1])
                    mu = float(params[key][2])
                    sigma2 = float(params[key][3])
                    result = gaussian(z_coords[z], A, mu, sigma2)
                    deltaEz[z] = result[0]
                    derivEz[z] = result[1]
                else:
                    A = float(params[key][4])
                    mu = float(params[key][5])
                    sigma2 = float(params[key][6])
                    result = gaussian(z_coords[z], A, mu, sigma2)
                    deltaEz[z] = result[0]
                    derivEz[z] = result[1]
        elif key[-1] == "3":
            #continue
            E0 = float(params[key][7])
            Zmid = float(params[key][8])
            N = float(params[key][9])
            for z in range(len(z_coords)):
                if z_coords[z] > float(params[key][12]):
                    result = sigmoidal(z_coords[z], E0, Zmid, N)
                    deltaEz[z] = result[0]
                    derivEz[z] = result[1]
                elif z_coords[z] <= float(params[key][12]):
                    A = float(params[key][1])
                    mu = float(params[key][2])
                    sigma2 = float(params[key][3])
                    result = gaussian(z_coords[z], A, mu, sigma2)
                    deltaEz[z] = result[0]
                    derivEz[z] = result[1]
        elif key[-1] == "4":
            #continue
            E0 = float(params[key][7])
            Zmid = float(params[key][8])
            N = float(params[key][9])
            for z in range(len(z_coords)):
                if z_coords[z] < float(params[key][13]):
                    result = sigmoidal(z_coords[z], E0, Zmid, N)
                    deltaEz[z] = result[0]
                    derivEz[z] = result[1]
                elif z_coords[z] >= float(params[key][13]):
                    A = float(params[key][4])
                    mu = float(params[key][5])
                    sigma2 = float(params[key][6])
                    result = gaussian(z_coords[z], A, mu, sigma2)
                    deltaEz[z] = result[0]
                    derivEz[z] = result[1]
        elif key[-1] == "5":
            #continue
            E0 = float(params[key][7])
            Zmid = float(params[key][8])
            N = float(params[key][9])
            for z in range(int(len(z_coords)/2)):
                result = sigmoidal(z_coords[z], E0, Zmid, N, asym = True)
                deltaEz[z] = result[0]
                derivEz[z] = result[1]
            Zmid = float(params[key][10])
            N = float(params[key][11])
            for z in range(int(len(z_coords)/2), len(z_coords)):
                result = sigmoidal(z_coords[z], E0, Zmid, N, asym = True)
                deltaEz[z] = result[0]
                derivEz[z] = result[1]
        else:
            #continue
            A = float(params[key][1])
            mu = float(params[key][2])
            sigma2 = float(params[key][3])
            for z in range(int(len(z_coords)/2)):
                result = gaussian(z_coords[z], A, mu, sigma2)
                deltaEz[z] = result[0]
                derivEz[z] = result[1]
            A = float(params[key][4])
            mu = float(params[key][5])
            sigma2 = float(params[key][6])
            for z in range(int(len(z_coords)/2), len(z_coords)):
                result = gaussian(z_coords[z], A, mu, sigma2)
                deltaEz[z] = result[0]
                derivEz[z] = result[1]
        #left end discontinuities to smooth
        if key in ["ASP2", "GLU3", "LYS2", "ARG3"]:
            right_end = int(np.where(z_coords == int(params[key][14]))[0])
            left_end = int(np.where(z_coords == int(params[key][12]))[0])
            slope = (deltaEz[right_end] - deltaEz[left_end]) / (int(params[key][14]) - int(params[key][12]))
            intercept = deltaEz[right_end] - slope*int(params[key][14])
            for z in range(left_end+1, right_end):
                deltaEz[z] = slope * z_coords[z] + intercept
                derivEz[z] = slope
        
        if key in ["HIS4", "ASP2", "ASN4"]:
            right_end = int(np.where(z_coords == int(params[key][13]))[0])
            left_end = int(np.where(z_coords == int(params[key][15]))[0])
            slope = (deltaEz[right_end] - deltaEz[left_end]) / (int(params[key][13]) - int(params[key][15]))
            intercept = deltaEz[right_end] - slope*int(params[key][13])
            for z in range(left_end+1, right_end):
                deltaEz[z] = slope * z_coords[z] + intercept
                derivEz[z] = slope
 
        if np.array_equal(deltaEz, np.zeros(len(z_coords))):
            continue
        else:
            #outData.write(key[:-1] + " " + " ".join(map(str, deltaEz)) + "\n")
            #derivs.write(key[:-1] + " " + " ".join(map(str, derivEz)) + "\n")
            #outData.write(key[:-1] + "\t" + "\t".join(map(str, deltaEz)) + "\n")
            #derivs.write(key[:-1] + "\t" + "\t".join(map(str, derivEz)) + "\n")
            outData.write(key[:-1])
            derivs.write(key[:-1])
            for num in deltaEz:
                outData.write("\t%.4f" % num)
            for num in derivEz:
                derivs.write("\t%.4f" % num)
            outData.write("\n")
            derivs.write("\n")
            #print(sigEz)
            #continue
            plt.figure(count)
            plt.scatter(z_coords, deltaEz)
            plt.title(key[:-1])
            plt.ylim(-2,2)
            plt.xlim(-30,30)
            plt.savefig("EnergyGraphs/2012EzPotential_%s_%s.png"%(key[:-1], term_atoms))


