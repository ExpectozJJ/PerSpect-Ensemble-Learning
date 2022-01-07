import csv
import sys
import os
from subprocess import call
import glob
from collections import defaultdict
import numpy as np
import math
import gudhi as gd
from scipy.sparse import *
from scipy import *

def faces(simplices):
    faceset = set()
    for simplex in simplices:
        numnodes = len(simplex)
        for r in range(numnodes, 0, -1):
            for face in combinations(simplex, r):
                faceset.add(tuple(sorted(face)))
    return faceset

def n_faces(face_set, n):
    return filter(lambda face: len(face)==n+1, face_set)

def boundary_operator(face_set, i):
    source_simplices = list(n_faces(face_set, i))
    target_simplices = list(n_faces(face_set, i-1))
    #print(source_simplices, target_simplices)

    if len(target_simplices)==0:
        S = dok_matrix((1, len(source_simplices)), dtype=np.float64)
        S[0, 0:len(source_simplices)] = 1
    else:
        source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}

        S = dok_matrix((len(target_simplices), len(source_simplices)), dtype=np.float64)
        for source_simplex in source_simplices:
            for a in range(len(source_simplex)):
                target_simplex = source_simplex[:a]+source_simplex[(a+1):]
                i = target_simplices_dict[target_simplex]
                j = source_simplices_dict[source_simplex]
                S[i, j] = -1 if a % 2==1 else 1
    
    return S

def compute_stat(dist):
    dist = np.array(dist)
    if len(dist) == 0:
        feat = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] # If the input is empty, output all-zero stat
    else:
        feat = []
        feat.append(np.count_nonzero(dist<1e-3))          # persistent multiplicity

        dist = dist[dist>=1e-3]
        if len(dist) == 0:
            feat = [feat[0],0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        else:
            feat.append(np.min(dist))                       # persistent min (Fiedler Value - Algebraic Connectivity)
            feat.append(np.max(dist))                       # persistent max
            feat.append(np.mean(dist))                      # persistent mean
            feat.append(np.std(dist))                       # persistent std

            feat.append(np.sum(np.abs(dist)))               # persistent Laplacian Graph Energy 
            #feat.append(len(dist))                          # persistent number of non-zero eigenvalues
            s, t, u, v, w, x, y, z = 0, 0, [], 0, 0, 0, 0, 0
            for l in dist:
                z += l*l*l*l*l
                y += l*l*l*l
                x += l*l*l
                w += l*l
                u.append(abs(l-feat[3])/len(dist))
                v += 1/l
                t += l**(-2)
                s += math.log(l)
            
            u = np.array(u)
            feat.append(np.sum(u))                           # persistent Laplacian Generalised Mean Graph Energy

            feat.append(w)                                   # persistent spectral 2nd moment
            #feat.append(x)                                   # persistent spectral 3rd moment
            #feat.append(y)                                   # persistent spectral 4th moment
            #feat.append(z)                                   # persistent spectral 5th moment
            feat.append(t)                                   # persistent zeta(2) of laplacian
            feat.append((len(dist)+1)*v)                     # persistent quasi-Wiener Index
            feat.append(s-math.log(len(dist)+1))             # persistent spanning tree number

    return feat

def read_eig(args):
    foldername = "./"+args[0] +'_'+args[1]+'_'+args[2]+'_'+args[3]+'_'+args[4]
    output = []

    print("{}".format(foldername+"/binding_PH"))
    os.chdir(foldername+"/binding_PH")
    file1 = open("wild_binding_eig.txt")
    file2 = open("mut_binding_eig.txt")
    c1, c2 = file1.readlines(), file2.readlines()
    wb = defaultdict(defaultdict)
    if len(c1) > 0:
        for i in range(len(c1)):
            tmp = c1[i].split()
            if tmp[0]+tmp[1] not in wb.keys():
                wb[tmp[0]+tmp[1]] = defaultdict(list)
            wb[tmp[0]+tmp[1]][float(tmp[-2])].append(float(tmp[-1]))
        
        for key, val in wb.items():
            for k1, subval in val.items():
                result = compute_stat(subval)
                #print(key, k1, result[0], np.shape(result))
                output.append(result)
    else:
        for ele in ["CC", "CO", "CN", "NC", "NN", "NO", "OC", "ON", "OO"]:
            for i in range(48):
                output.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

    mb = defaultdict(defaultdict)
    if len(c2) > 0:
        for i in range(len(c2)):
            tmp = c2[i].split()
            if tmp[0]+tmp[1] not in mb.keys():
                mb[tmp[0]+tmp[1]] = defaultdict(list)
            mb[tmp[0]+tmp[1]][float(tmp[-2])].append(float(tmp[-1]))
    
        for key, val in mb.items():
            for k1, subval in val.items():
                result = compute_stat(subval)
                #print(key, k1, result[0], np.shape(result))
                output.append(result)
    else:
        for ele in ["CC", "CO", "CN", "NC", "NN", "NO", "OC", "ON", "OO"]:
            for i in range(48):
                output.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    
    os.chdir('..')
    os.chdir('..')

    print("{}".format(foldername+"/mutation_PH"))
    os.chdir(foldername+"/mutation_PH")
    file1 = open("wild_mutation_eig.txt")
    file2 = open("mut_mutation_eig.txt")
    c1, c2 = file1.readlines(), file2.readlines()
    wb = defaultdict(defaultdict)
    if len(c1) > 0:
        for i in range(len(c1)):
            tmp = c1[i].split()
            if tmp[0]+tmp[1] not in wb.keys():
                wb[tmp[0]+tmp[1]] = defaultdict(list)
            wb[tmp[0]+tmp[1]][float(tmp[-2])].append(float(tmp[-1]))
        
        for key, val in wb.items():
            for k1, subval in val.items():
                result = compute_stat(subval)
                #print(key, k1, result[0], np.shape(result))
                output.append(result)
    else:
        for ele in ["CC", "CO", "CN", "NC", "NN", "NO", "OC", "ON", "OO"]:
            for i in range(48):
                output.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

    mb = defaultdict(defaultdict)
    if len(c2) > 0:
        for i in range(len(c2)):
            tmp = c2[i].split()
            if tmp[0]+tmp[1] not in mb.keys():
                mb[tmp[0]+tmp[1]] = defaultdict(list)
            mb[tmp[0]+tmp[1]][float(tmp[-2])].append(float(tmp[-1]))
        
        for key, val in mb.items():
            for k1, subval in val.items():
                result = compute_stat(subval)
                #print(key, k1, result[0], np.shape(result))
                output.append(result)

    else:
        for ele in ["CC", "CO", "CN", "NC", "NN", "NO", "OC", "ON", "OO"]:
            for i in range(48):
                output.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

    os.chdir('..')
    os.chdir('..')

    return output

def call_script(args):
    os.system("./run_topnettree.sh {}".format(args))

def hodge(args):
    foldername = "./"+args[0] +'_'+args[1]+'_'+args[2]+'_'+args[3]+'_'+args[4]
    os.system("echo {}".format(foldername))
    os.chdir(foldername)
    #os.chdir('..')
    #os.system("python2 Gen_pointclouds.py {} {} {}".format(args[0],args[1],args[2]))
    #os.system("rm ./binding_PH/*.txt")
    #os.system("rm ./mutation_PH/*.txt")
    os.system("rm gen_hodge.m")
    os.system("rm computeVRComplex.m")
    os.system("rm gen_hodge.py")
    os.chdir('..')

    os.system("cp ../../computeVRComplex.m " + foldername+"/binding_PH")
    os.system("cp ../../binding_hodge.m " + foldername+"/binding_PH")
    os.system("echo {}".format(foldername+"/binding_PH"))
    os.chdir(foldername+"/binding_PH")
    os.system("matlab -nodisplay -r 'run binding_hodge.m'")
    os.chdir('..')
    os.chdir('..')

    os.system("cp ../../computeVRComplex.m " + foldername+"/mutation_PH")
    os.system("cp ../../mutation_hodge.m " + foldername+"/mutation_PH")
    os.system("echo {}".format(foldername+"/mutation_PH"))
    os.chdir(foldername+"/mutation_PH")
    os.system("matlab -nodisplay -r 'run mutation_hodge.m'")

    #os.system("python gen_hodge.py")
    os.chdir('..')
    os.chdir('..')

def alpha(args):
    foldername = "./"+args[0] +'_'+args[1]+'_'+args[2]+'_'+args[3]+'_'+args[4]
    os.system("echo {}".format(foldername))

    pm1, pm2 = [], []
    os.chdir(foldername+"/binding_PH")
    for typ in ['mut', 'wild']:
        A, B = defaultdict(list), defaultdict(list)
        for ele1 in ['C', 'N', 'O']:
            with open("{}_binding_C_{}.pts".format(typ, ele1)) as csvfile:
                data = np.array(list(csv.reader(csvfile, delimiter=' ', quotechar='|')), dtype=float)
                csvfile.close()

            if len(data)>0:
                A[ele1] = data[data[:, 0]==0][:, 1:]

            with open("{}_binding_{}_C.pts".format(typ, ele1)) as csvfile:
                data = np.array(list(csv.reader(csvfile, delimiter=' ', quotechar='|')), dtype=float)
                csvfile.close()

            if len(data)>0:
                B[ele1] = data[data[:, 0]==1][:, 1:]
            
        for e_set in [['C'], ['N'], ['O'], ['C', 'N'], ['C', 'O'], ['N', 'O'], ['C', 'N', 'O']]:
            Acloud = A[e_set[0]]
            Bcloud = B[e_set[0]]
            for i in range(1, len(e_set)):
                Acloud = np.concatenate((Acloud, A[e_set[i]]), axis=0)
                Bcloud = np.concatenate((Bcloud, B[e_set[i]]), axis=0)
                
            #print(typ, e_set, len(Acloud))
            #print(typ, e_set, len(Bcloud))
                
            alpha_complex = gd.AlphaComplex(Acloud)
            st = alpha_complex.create_simplex_tree()
            dgmsalpha = st.persistence()
            betti0, betti1, betti2 = [], [], []
            for r in dgmsalpha:
                if r[0] == 0:
                    betti0.append([r[1][0], r[1][1]])
                elif r[0] == 1:
                    betti1.append([r[1][0], r[1][1]])
                elif r[0] == 2:
                    betti2.append([r[1][0], r[1][1]])

            # Using circumradius, we take sqrt of F and multiply by 2  
            betti0 = np.array(np.sqrt(betti0)*2)
            betti1 = np.array(np.sqrt(betti1)*2)
            betti2 = np.array(np.sqrt(betti2)*2)
            betti = [betti0, betti1, betti2]

            betti0 = sorted(betti[0], key=lambda x: x[0])
            betti0 = np.flip(betti0, axis=0)
            betti1 = sorted(betti[1], key=lambda x: x[0])
            betti1 = np.flip(betti1, axis=0)
            betti2 = sorted(betti[2], key=lambda x: x[0])
            betti2 = np.flip(betti2, axis=0)

            pbn1, pbn2 = np.zeros(80), np.zeros(80)
            for i in range(10, 90):
                for j in range(len(betti1)):
                    if betti1[j][0] <= 0.1*(i) and betti1[j][1] >= 0.1*(i+1):
                        pbn1[i-10] += 1
                for j in range(len(betti2)):
                    if betti2[j][0] <= 0.1*(i) and betti2[j][1] >= 0.1*(i+1):
                        pbn2[i-10] += 1

            pm1.append(pbn1)
            pm2.append(pbn2)
            
            alpha_complex = gd.AlphaComplex(Bcloud)
            st = alpha_complex.create_simplex_tree()
            dgmsalpha = st.persistence()
            betti0, betti1, betti2 = [], [], []
            for r in dgmsalpha:
                if r[0] == 0:
                    betti0.append([r[1][0], r[1][1]])
                elif r[0] == 1:
                    betti1.append([r[1][0], r[1][1]])
                elif r[0] == 2:
                    betti2.append([r[1][0], r[1][1]])

            # Using circumradius, we take sqrt of F and multiply by 2  
            betti0 = np.array(np.sqrt(betti0)*2)
            betti1 = np.array(np.sqrt(betti1)*2)
            betti2 = np.array(np.sqrt(betti2)*2)
            betti = [betti0, betti1, betti2]

            betti0 = sorted(betti[0], key=lambda x: x[0])
            betti0 = np.flip(betti0, axis=0)
            betti1 = sorted(betti[1], key=lambda x: x[0])
            betti1 = np.flip(betti1, axis=0)
            betti2 = sorted(betti[2], key=lambda x: x[0])
            betti2 = np.flip(betti2, axis=0)

            pbn1, pbn2 = np.zeros(80), np.zeros(80)
            for i in range(10, 90):
                for j in range(len(betti1)):
                    if betti1[j][0] <= 0.1*(i) and betti1[j][1] >= 0.1*(i+1):
                        pbn1[i-10] += 1
                for j in range(len(betti2)):
                    if betti2[j][0] <= 0.1*(i) and betti2[j][1] >= 0.1*(i+1):
                        pbn2[i-10] += 1

            pm1.append(pbn1)
            pm2.append(pbn2)
    
    os.chdir("..")
    os.chdir("..")

    os.chdir(foldername+"/mutation_PH")
    for typ in ['mut', 'wild']:
        A, B = defaultdict(list), defaultdict(list)
        for ele1 in ['C', 'N', 'O']:
            with open("{}_mutation_C_{}.pts".format(typ, ele1)) as csvfile:
                data = np.array(list(csv.reader(csvfile, delimiter=' ', quotechar='|')), dtype=float)
                csvfile.close()

            if len(data)>0:
                A[ele1] = data[data[:, 0]==0][:, 1:]

            with open("{}_mutation_{}_C.pts".format(typ, ele1)) as csvfile:
                data = np.array(list(csv.reader(csvfile, delimiter=' ', quotechar='|')), dtype=float)
                csvfile.close()

            if len(data)>0:
                B[ele1] = data[data[:, 0]==1][:, 1:]
            
        for e_set in [['C'], ['N'], ['O'], ['C', 'N'], ['C', 'O'], ['N', 'O'], ['C', 'N', 'O']]:
            Acloud = A[e_set[0]]
            Bcloud = B[e_set[0]]
            for i in range(1, len(e_set)):
                Acloud = np.concatenate((Acloud, A[e_set[i]]), axis=0)
                Bcloud = np.concatenate((Bcloud, B[e_set[i]]), axis=0)
                
            #print(typ, e_set, len(Acloud))
            #print(typ, e_set, len(Bcloud))
                
            alpha_complex = gd.AlphaComplex(Acloud)
            st = alpha_complex.create_simplex_tree()
            dgmsalpha = st.persistence()
            betti0, betti1, betti2 = [], [], []
            for r in dgmsalpha:
                if r[0] == 0:
                    betti0.append([r[1][0], r[1][1]])
                elif r[0] == 1:
                    betti1.append([r[1][0], r[1][1]])
                elif r[0] == 2:
                    betti2.append([r[1][0], r[1][1]])

            # Using circumradius, we take sqrt of F and multiply by 2  
            betti0 = np.array(np.sqrt(betti0)*2)
            betti1 = np.array(np.sqrt(betti1)*2)
            betti2 = np.array(np.sqrt(betti2)*2)
            betti = [betti0, betti1, betti2]

            betti0 = sorted(betti[0], key=lambda x: x[0])
            betti0 = np.flip(betti0, axis=0)
            betti1 = sorted(betti[1], key=lambda x: x[0])
            betti1 = np.flip(betti1, axis=0)
            betti2 = sorted(betti[2], key=lambda x: x[0])
            betti2 = np.flip(betti2, axis=0)

            pbn1, pbn2 = np.zeros(80), np.zeros(80)
            for i in range(10, 90):
                for j in range(len(betti1)):
                    if betti1[j][0] <= 0.1*(i) and betti1[j][1] >= 0.1*(i+1):
                        pbn1[i-10] += 1
                for j in range(len(betti2)):
                    if betti2[j][0] <= 0.1*(i) and betti2[j][1] >= 0.1*(i+1):
                        pbn2[i-10] += 1

            pm1.append(pbn1)
            pm2.append(pbn2)
            
            alpha_complex = gd.AlphaComplex(Bcloud)
            st = alpha_complex.create_simplex_tree()
            dgmsalpha = st.persistence()
            betti0, betti1, betti2 = [], [], []
            for r in dgmsalpha:
                if r[0] == 0:
                    betti0.append([r[1][0], r[1][1]])
                elif r[0] == 1:
                    betti1.append([r[1][0], r[1][1]])
                elif r[0] == 2:
                    betti2.append([r[1][0], r[1][1]])

            # Using circumradius, we take sqrt of F and multiply by 2  
            betti0 = np.array(np.sqrt(betti0)*2)
            betti1 = np.array(np.sqrt(betti1)*2)
            betti2 = np.array(np.sqrt(betti2)*2)
            betti = [betti0, betti1, betti2]

            betti0 = sorted(betti[0], key=lambda x: x[0])
            betti0 = np.flip(betti0, axis=0)
            betti1 = sorted(betti[1], key=lambda x: x[0])
            betti1 = np.flip(betti1, axis=0)
            betti2 = sorted(betti[2], key=lambda x: x[0])
            betti2 = np.flip(betti2, axis=0)

            pbn1, pbn2 = np.zeros(80), np.zeros(80)
            for i in range(10, 90):
                for j in range(len(betti1)):
                    if betti1[j][0] <= 0.1*(i) and betti1[j][1] >= 0.1*(i+1):
                        pbn1[i-10] += 1
                for j in range(len(betti2)):
                    if betti2[j][0] <= 0.1*(i) and betti2[j][1] >= 0.1*(i+1):
                        pbn2[i-10] += 1

            pm1.append(pbn1)
            pm2.append(pbn2)
    
    os.chdir("..")
    os.chdir("..")
    
    return np.concatenate((pm1, pm2), axis=0)

#j = int(sys.argv[1])
start, end = int(sys.argv[1]), int(sys.argv[2])

with open('../../skempi.csv') as csvfile:
    mat = []
    data = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in data:
        mat.append(', '.join(row))
csvfile.close()

#print(mat)
#print(len(mat))

args = []
cnt = 0
for i in range(1, len(mat)):

    #if i-1 == j:
    if i >= start and i < end:
        tmp = mat[i].split(",")
        pdb = tmp[0]+".pdb"
        chainid = tmp[1]
        wildtype = tmp[2]
        mutanttype = tmp[4]
        resid = tmp[3]
        #if tmp[0] == "2I9B":
        args.append([tmp[0], chainid, wildtype, resid, mutanttype])
        #call_script(pdb+" "+chainid+" "+wildtype+" "+mutanttype+" "+resid+" 0")
        #print(args[-1])
        cnt += 1

os.system("echo Count: {}".format(cnt))


feature = []
os.chdir("./src/feature")

for i in range(len(args)):
    #hodge(args[i])
    op = read_eig(args[i])
    feature.append(op)
    print(np.shape(feature))

os.chdir("..")
os.chdir("..")
np.save("X_skempi_l0.npy", np.transpose(np.reshape(feature,(len(feature), 36, 48, 15)), axes=(0,2,1,3)))

"""
os.chdir("./src/feature")
alpha_feat = []
for i in range(len(args)):
    alpha_feat.append(alpha(args[i]))
    print(np.shape(alpha_feat))

os.chdir("..")
os.chdir("..")
np.save("X_skempi_alpha_l1_l2.npy", np.transpose(alpha_feat, axes=(0, 2, 1)))
"""
