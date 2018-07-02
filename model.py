# Dependencies
from __future__ import division
import math
import numpy as np
from pylab import plot, show, xlabel, ylabel
import random
from Bio import motifs
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import string

# DATA SET : 18 SAMPLES OF -35 AND -10 SEQUENCES WITH RESPECTIVE STRENGTH

#filename1 = raw_input("Enter /path/to/filename with promoter sequences of regulon db (csv of -35 hexamer)\t")
filename1 = 'Datasets/fimo_35_ff_ed11.csv'
#filename2 = raw_input("Enter /path/to/filename with promoter sequences of regulon db (csv of -10 hexamer)\t")
filename2 = 'Datasets/fimo_10_ff_ed11.csv'

infile1 = open(filename1, 'r')
infile2 = open(filename2, 'r')

instances_regdb = list()
instances2_regdb = list()
#outputResult = list()

for line in infile1.readlines():
   if line == '\n':
      break
   lines = line.split(',')
   lines[1] = lines[1][:-1]
#   print lines
   instances_regdb.append(Seq(lines[1].upper()))
#   outputResult.append([float(lines[3])])

for line in infile2.readlines():
   if line == '\n':
      break
   lines = line.split(',')
   lines[1] = lines[1][:-1]
#   print lines
   instances2_regdb.append(Seq(lines[1].upper()))
   
infile1.close()
infile2.close()
##instances, instances2 and outputResult -- all are instantiated for use in rest of the program

# -35 SEQUENCE
instances = [
Seq("TTGACG"),
Seq("TTTACA"),
Seq("TTGACA"),
Seq("TTGACA"),
Seq("TTTACG"),
Seq("TTTACG"),
Seq("TTTACG"),
Seq("CTGACA"),
Seq("TTTACA"),
Seq("TTTACG"),
Seq("TTGACG"),
Seq("CTGATA"),
Seq("CTGATG"),
Seq("TTTATG"),
Seq("TTTATA"),
Seq("TTGACA"),
Seq("TTGACA"),
Seq("TTGACG")
]
# -10 SEQUENCE
instances2 =[
Seq("TACAGT"),
Seq("TATTAT"),
Seq("TACTGT"),
Seq("TATTGT"),
Seq("TACTAT"),
Seq("TATAGT"),
Seq("TATTAT"),
Seq("TATAAT"),
Seq("GACTGT"),
Seq("TACAAT"),
Seq("TATAGT"),
Seq("GATTAT"),
Seq("GATTAT"),
Seq("TACAAT"),
Seq("TACAAT"),
Seq("GACTAT"),
Seq("GATTGT"),
Seq("TATTGT")
]
# STRENGTH OF -35 AND -10 SEQUENCE
outputResult = [
[1],
[0.7],
[0.86],
[0.72],
[0.24],
[0.47],
[0.36],
[0.51],
[0.04],
[0.33],
[0.58],
[0.01],
[0.01],
[0.1],
[0.15],
[0.16],
[0.06],
[0.56]
]

flag = 0
convert10 = [] #variable for str representation of -10 hexamer for PWM_10 scores generation in regression modelling
convert35 = [] #variable for str representation of -35 hexamer for PWM_35 scores generation in regression modelling
result = []
result2 = []

# GRADIENT DESCENT TTO OBTAIN THE THETA PARAMETER
def gradientDescent(x, y, theta, alpha, m, numIterations):
    J_history = np.zeros(shape=(numIterations, 1))
    xTrans = x.transpose()
    for i in range(0, numIterations):
        hypothesis = np.dot(x, theta)
        loss = hypothesis - y
        cost = np.sum(loss ** 2) / (2 * m)
        #print("Iteration %d | Cost: %f" % (i, cost))
        gradient = np.dot(xTrans, loss) / m
        theta = theta - alpha * gradient
        J_history[i][0] = cost
    return theta, J_history

#choice = raw_input("Do you want to enter pseudocounts or use defaults (Y / N):\t")

choice = raw_input("Would you like to add your newly characterized promoters to the dataset (Y / N):\t")
# print choice

if (choice == 'y' or choice == 'Y'):
    numberOfDatasets = int(raw_input("How many datasets to add?\t"))
    for i in range(0, numberOfDatasets):
        while True:
            sequence35 = raw_input("Enter the -35 Hexamer Sequence:\t").upper()
            if len(sequence35)!=6:
                print "Please enter a hexamer sequence (exactly six nucleotides)"
            else:
                break
        while True:
            sequence10 = raw_input("Enter the -10 Hexamer Sequence:\t").upper()
            if len(sequence10)!=6:
                print "Please enter a hexamer sequence (exactly six nucleotides)"
            else:
                break
        strength   = float(raw_input("Enter the strength"))
        instances.append(Seq(sequence35))
        instances2.append(Seq(sequence10))
        outputResult.append([strength])
else:
    pass

choice = raw_input("Would you like to predict the strength of a promoter region (Y / N):\t")
# print choice

if (choice == 'y' or choice == 'Y'):
    while True:
        sequence35 = raw_input("Enter the -35 Hexamer Sequence:\t").upper()
        if len(sequence35)!=6:
            print "Please enter sequence of exactly six nucleotides"
        else:
            break
    while True:
        sequence10 = raw_input("Enter the -10 Hexamer Sequence:\t").upper()
        if len(sequence10)!=6:
            print "Please enter a sequence of exactly six nucleotides"
        else:
            break
    flag = 1
else:
    pass

# CONVERTING THE -35 SEQUENCE INTO A SUITABLE FORMAT
i = 0
while i < len(instances):
    convert35.append(str(instances[i]))
    i+=1

# GENERATIVE MODELLING -- CONSTRUCTION OF THE PSSM MATRIX FOR THE -35 SEQUENCE BASED ON STRENGTH
m = motifs.create(instances_regdb[:])
#print(m.counts);
pwm = m.counts.normalize(pseudocounts= {'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49}                )
#print(pwm)
pssm = pwm.log_odds()
#print(pssm)

#REGRESSION MODELLING
if flag == 1:
    p,o,l,k,m,n = str(sequence35)
    resultP = pssm[p,0] + pssm[o,1] + pssm[l,2] + pssm[k,3] + pssm[m,4] + pssm[n,5]

def calculateX(a,b,c,d,e,f,x):
    temp1 = pssm[a,0] + pssm[b,1] + pssm[c,2] + pssm[d,3] + pssm[e,4] + pssm[f,5]
    result.append([temp1])


i = 0
while i < len(convert35):
    calculateX(convert35[i][0],convert35[i][1],convert35[i][2],convert35[i][3],convert35[i][4],convert35[i][5],i)
    i +=1

# CONVERTING THE -10 SEQUENCE INTO A SUITABLE FORMAT
i = 0
while i < len(instances2):
    convert10.append(str(instances2[i]))
    i+=1

# GENERATIVE MODELLING -- CONSTRUCTION OF THE PSSM MATRIX FOR THE -35 SEQUENCE
m2 = motifs.create(instances2_regdb)
#print(m2.counts);
pwm2 = m2.counts.normalize(pseudocounts={'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49})
#print(pwm2)
pssm2 = pwm2.log_odds()
#print(pssm2)

#REGRESSION MODELLING
if flag == 1:
    p2,o2,l2,k2,m2,n2 = str(sequence10)
    resultP2 = pssm2[p2,0] + pssm2[o2,1] + pssm2[l2,2] + pssm2[k2,3] + pssm2[m2,4] + pssm2[n2,5]

def calculateX2(a,b,c,d,e,f,x):
    temp1 = pssm2[a,0] + pssm2[b,1] + pssm2[c,2] + pssm2[d,3] + pssm2[e,4] + pssm2[f,5]
    result2.append([temp1])

i = 0
while i < len(convert10):
    calculateX2(convert10[i][0],convert10[i][1],convert10[i][2],convert10[i][3],convert10[i][4],convert10[i][5],i)
    i +=1
# CONSTRUCTION OF MATRIX A - PSSM OF -35 AND -10 SEQUENCE
a = []
i = 0
while i<len(outputResult):
    a.append([1,result[i][0],result2[i][0]])
    i +=1

# CONSTRUCTION OF B MATRIX - STRENGTH OF PROMOTER
b =[]
i = 0
# print len(outputResult)
while i<len(outputResult):
    b += outputResult[i]
    # TO DEAL WITH LOG(0)
    if b[i] < 0.0001:
        b[i] = 0.0001  # Might change cutoff
    b[i] = math.log(b[i])
    i +=1

# outfile2 = open('pwmScores.csv','w')
# for i in range(len(outputResult)):
#     record = str(result[i][0])+','+ str(result2[i][0]) + ',' + str(result[i][0]+result2[i][0]) + ',' +str(outputResult[i][0]) +','+ str(b[i])+'\n'
#     outfile2.write(record)
# outfile2.close()

x = np.asarray(a)
y = np.asarray(b)

# CALLING THE GRADIENT DESCENT FUNCTION TO OBTAIN THE THETA PARAMETER
m, n = np.shape(x)
numIterations= 100000 #c
alpha = 0.015 #c
theta = np.ones(n)
theta, J_history = gradientDescent(x, y, theta, alpha,m,numIterations)

print "Theta", theta

# CONSTRUCTING THE HYPOTHESIS
hx = x.dot(theta)
#print "hx x.dot(theta):\t", hx
diff = hx - y
# print diff

diff_square = diff*diff
# print diff_square

sum = np.sum(diff_square)
# print sum 


# COST FUNCTION
temp = 1/(2*m)
cost = temp*sum

# TO OBTAIN R SQUARE VALUE : CORRELATION COEFFICIENT
meany = np.mean(y)
sumsqmeany = np.sum((y-meany)**2)
sumsqmeanysum = np.sum((y-hx)**2)/sumsqmeany
# print sumsqmeanysum

print ""
print "RESULTS"
print "======="

if flag == 1:
    strength = np.array([1.0, resultP, resultP2 ]).dot(theta)
    print "Predicted Strength:\t\t",
    print(math.exp(strength))
    print "Predicted ln(Strength):\t\t",
    print(strength)

print "R-squared value of model:\t",
R = 1 - sumsqmeanysum
print R
print "Adj. R-squared value of model:\t",
adjR = 1 - (1-R)*(len(instances)-1)/(len(instances)-2-1)
print adjR
print ""

# CONSTRUCTION OF THE MULTIVARIENT LINEAR REGRESSION GRAPH
# THE NUMBER 18 REPRESENTS THE DEFAULT SIZE OF THE DATASET PROVIDED
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
for c, m in [('r','o')]:
    xs = x[0:18,1]
##    print "\t\t\t\t xs Plot for Graph -10 Sequence"
##    print(xs)
##    print ""
    ys = x[0:18,2]
##    print "\t\t\t\t ys Plot for Graph -35 Sequence"
##    print (ys)
##    print ""
    zs = y[0:18]
##    print "\t\t\t\t zs Plot for Graph Strength"
##    print (zs)
##    print ""
    ax.scatter(xs, ys,zs, c=c, marker =m)

md = len(y)
if (md > 18):
    flag2 = 1
    for c,m in [('b','o')]:
        xd = x[18:,1]
        yd = x[18:,2]
        zd = y[18:]
        ax.scatter(xd, yd, zd, c=c, marker =m)

ax.set_xlabel('-10 Hexamer')
ax.set_ylabel('-35 Hexamer')
ax.set_zlabel('Strength of Promoter')
plt.show()
