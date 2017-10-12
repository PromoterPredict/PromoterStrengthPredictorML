# Dependencies
from __future__ import division
import math
import numpy as np
from pylab import plot, show, xlabel, ylabel
import random
import web
from Bio import motifs
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

# DATA SET : 19 SAMPLES OF -35 AND -10 SEQUENCES WITH RESPECTIVE STRENGTH

# -35 SEQUENCE
instances = [
Seq("TTGACG"),
Seq("TTTACA"),
Seq("TTGACA"),
Seq("CTGATA"),
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
Seq("GATTAT"),
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
[0.01],
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

# STRENGTH OF -35 AND -10 SEQUENCE  - COPY OF outputResult
outputResult2 = [
[1],
[0.7],
[0.86],
[0.01],
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
convert10 = []
convert35 = []
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

choice = raw_input("Do you want to enter more Data to the existing Datasets (Y / N)")
print choice

if (choice == 'y' or choice == 'Y'):
    numberOfDatasets = int(raw_input("How many Datasets do you want to Add"))
    for i in range(0, numberOfDatasets):
        sequence35 = raw_input("Enter the -35 Hexamer Sequence")
        sequence10 = raw_input("Enter the -10 Hexamer Sequence")
        strength   = int(raw_input("Enter the strength"))
        instances.append(Seq(sequence35))
        instances2.append(Seq(sequence10))
        outputResult.append([strength])

choice = raw_input("Do you want to Predict the strength of a Promoter (Y / N)")
print choice

if (choice == 'y' or choice == 'Y'):
    sequence35 = raw_input("Enter the -35 Hexamer Sequence")
    sequence10 = raw_input("Enter the -10 Hexamer Sequence")
    instances.append(Seq(sequence35))
    instances2.append(Seq(sequence10))
    instancesP = instances[:]
    instancesP2 = instances2[:]
    flag = 1

# CONVERTING THE -35 SEQUENCE INTO A SUITABLE FORMAT
i = 0
while i < len(instances):
    convert35.append(str(instances[i]))
    i+=1

# CONVERTING THE -10 SEQUENCE INTO A SUITABLE FORMAT
i = 0
while i < len(instances2):
    convert10.append(str(instances2[i]))
    i+=1

# CONSTRUCTION OF THE PSSM MATRIX FOR THE -35 SEQUENCE BASED ON STRENGTH
m = motifs.create(instances)
#print(m.counts);
pwm = m.counts.normalize(pseudocounts= {'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49}                )
#print(pwm)
pssm = pwm.log_odds()
#print(pssm)

if flag == 1:
    mP = motifs.create(instancesP)
    pwmP = mP.counts.normalize(pseudocounts={'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49})
    pssmP = pwmP.log_odds()
    p,o,l,k,m,n = str(sequence35)
    resultP = pssmP[p,0] + pssmP[o,1] + pssmP[l,2] + pssmP[k,3] + pssmP[m,4] + pssmP[n,5]


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

# CONSTRUCTION OF THE PSSM MATRIX FOR THE -35 SEQUENCE BASED ON STRENGTH
m2 = motifs.create(instances2)
#print(m2.counts);
pwm2 = m2.counts.normalize(pseudocounts={'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49})
#print(pwm2)
pssm2 = pwm2.log_odds()
#print(pssm2)

if flag == 1:
    mP2 = motifs.create(instancesP2)
    pwmP2 = mP2.counts.normalize(pseudocounts={'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49})
    pssmP2 = pwmP2.log_odds()
    p2,o2,l2,k2,m2,n2 = str(sequence10)
    resultP2 = pssmP2[p2,0] + pssmP2[o2,1] + pssmP2[l2,2] + pssmP2[k2,3] + pssmP2[m2,4] + pssmP2[n2,5]


def calculateX2(a,b,c,d,e,f,x):
    temp1 = pssm2[a,0] + pssm2[b,1] + pssm2[c,2] + pssm2[d,3] + pssm2[e,4] + pssm2[f,5]
    result2.append([temp1])

i = 0
while i < len(convert10):
    calculateX2(convert10[i][0],convert10[i][1],convert10[i][2],convert10[i][3],convert10[i][4],convert10[i][5],i)
    i +=1

# CONSTRUCTION OF A MATRIX - PSSM OF -35 AND -10 SEQUENCE
a = []
i = 0
while i<len(outputResult):
    a.append([1,result[i][0],result2[i][0]])
    i +=1

print ""
print "\t\t\t\t Matrix A (Input Matrix : -35 and -10 Sequence)"
for x in a:
    print x
    print ""

# CONSTRUCTION OF B MATRIX - STRENGTH OF -35 AND -10 SEQUENCE
b =[]
i = 0
print len(outputResult)
while i<len(outputResult):
    b += outputResult[i]
    # TO DEAL WITH LOG(0)
    if b[i] == 0:
        b[i] = 0.01
    b[i] = math.log(b[i])
    i +=1

print ""
print "\t\t\t\t Matrix B (Strength)"
for x in b:
    print x
    print ""

# FORMATTING MATRIX A AND B TO OBTAIN MATRIX x AND y
print "hi"
print "\t\t\t\t Matrix X (Input Matrix : -35 and -10 Sequence)"
x = np.asarray(a)
print x
print ""
print "\t\t\t\t Matrix Y (Output Matrix : Strength)"
y = np.asarray(b)
print y
print ""
print "lol"

# CALLING THE GRADIENT DESCENT FUNCTION TO OBTAIN THE THETA PARAMETER
m, n = np.shape(x)
numIterations= 100000 #c
alpha = 0.015 #c
theta = np.ones(n)
theta, J_history = gradientDescent(x, y, theta, alpha,m,numIterations)
print "\t\t\t\t Theta : theta"
print(theta)
print ""

# CONSTRUCTING THE HYPOTHESIS
hx = x.dot(theta)
print "\t\t\t\t Hypothesis"
print hx
print ""

# DIFFERENCE BETWEEN THE HYPOTHESIS AND OUTPUT MATRIX Y
print "\t\t\t\t Difference"
diff = hx - y
print diff
print ""

# DIFFERENCE SQUARED
print "\t\t\t\t Difference Square"
diff_square = diff*diff
print diff_square
print ""

# SUM OF DIFFERENCES
print "\t\t\t\t Sum"
sum = np.sum(diff_square)
print sum
print ""

# COST FUNCTION
print "\t\t\t\t Cost function"
temp = 1/(2*m)
cost = temp*sum
print cost
print ""

# TO OBTAIN R SQUARE VALUE : CORRELATION COEFFICIENT
print "\t\t\t\t Cost Function"
meany = np.mean(y)
sumsqmeany = np.sum((y-meany)**2)
sumsqmeanysum = np.sum((y-hx)**2)/sumsqmeany
print sumsqmeanysum
print ""

print "\t\t\t\t R Sqare Value"
R = 1 - sumsqmeanysum
print R
print ""

# CONSTRUCTION OF THE MULTIVARIENT LINEAR REGRESSION GRAPH
# THE NUMBER 19 REPRESENTS THE DEFAULT SIZE OF THE DATASET PROVIDED
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
for c, m in [('r','o')]:
    xs = x[0:19,1]
    print "\t\t\t\t xs Plot for Graph -10 Sequence"
    print(xs)
    print ""
    ys = x[0:19,2]
    print "\t\t\t\t ys Plot for Graph -35 Sequence"
    print (ys)
    print ""
    zs = y[0:19]
    print "\t\t\t\t zs Plot for Graph Strength"
    print (zs)
    print ""
    ax.scatter(xs, ys,zs, c=c, marker =m)

md = len(y)
print "\t\t\t\t Total Number of Elements in the Dataset"
print md
print ""
if (md > 19):
    flag2 = 1
    for c,m in [('b','o')]:
        xd = x[19:,1]
        yd = x[19:,2]
        zd = y[19:]
        ax.scatter(xd, yd, zd, c=c, marker =m)

ax.set_xlabel('-10 Hexamer')
ax.set_ylabel('-35 Hexamer')
ax.set_zlabel('Strength of Promoter')
plt.show()
