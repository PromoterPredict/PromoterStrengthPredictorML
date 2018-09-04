# Dependencies

from __future__ import division
import web
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
import string

#                                          NOTE
#           HOW THE MODEL WORKS FOR PREDICTING AND ALLOWING THE USER TO ADD MORE DATASETS

# OBTAIN ALL THE INPUTS FROM THE CLIENT AS A SINGLE STRING : Eg TGTGTGTGATGA*TGATGATGTGTG#0.5&
# [THIS IS JUST A RANDOM SEQENCE TAKEN FOR THE SAKE OF THIS EXAMPLE]

# EXAMPLE OF STRING SENT FROM SERVER : TGTGTGTGATGA*TGATGATGTGTG#0.5&

# -35 AND -10 SEQUENCE THE USER ENTERS FOR WHICH HE REQUIRES THE PREDICTED STRENGTH (TGTGTGTGATGA*)
# -35 SEQUENCE : TGTGTG
# -10 SEQUENCE : TGATGA

# SEPERATOR *

# -35 AND -10 SEQUENCE THE USER ENTERS AS DYNAMIC INPUT WHICH IS TO BE ADDED TO THE DATA SET (TGATGATGTGTG#)
# -35 SEQUENCE : TGATGA
# -10 SEQUENCE : TGTGTG
# SUPPOSE THE USER ENTERS MORE THAN 1 DYNAMIC INPUT, LETS SAY 2 INPUT SETS (TGATGATGTGTGATGATGGATTGA*)
# SEQUENCE 1 : TGATGATGTGTG
# SEQUENCE 2 : ATGATGGATTGA
# -35 SEQUENCE1 : TGATGA
# -10 SEQUENCE1 : TGTGTG
# -35 SEQUENCE2 : ATGATG
# -10 SEQUENCE2 : GATTGA

# SEPERATOR #

# THE STRENGTH FOR THE RESPECTIVE -35 AND -10 SEQUENCE THE USER ENTERS AS DYNAMIC INPUT WHICH IS TO BE ADDED TO THE DATA SET(0.5&)
# STRENGTH : 0.5
# SUPPOSE THE USER ENTERS MORE THAN 1 DYNAMIC INPUT, LETS SAY 2 INPUT SETS (0.5&0.6&)
# STRENGTH1 : 0.5
# STRENGTH2 : 0.6

# RETURNS RESUTS TO SERVER
def make_text(string):
    return string

# TO CONNECT WITH THE WEBPAGE BEING USED
urls = ('/', 'tutorial')
render = web.template.render('templates/')

app = web.application(urls, globals())

my_form = web.form.Form(
                web.form.Textbox('', class_='textfield', id='textfield'),
                )

class tutorial:
    def GET(self):
        form = my_form()
        return render.tutorial(form, "")

    def POST(self):

        # ALL THE VARIABLES USED
        # CONSIDER EXAMPLE : TGTGTGTGATGA*TGATGATGTGTG#0.5&
        form = my_form()
        form.validates()
        t = form.value['textfield']
        numberOfElements = 0 # 2 x no of datasets added by user
        p = "" # Prediction Sequences -35 and -10 before * Eg TGTGTGTGATGA
        s1 = "" # Prediction Sequence -35 TGTGTG
        s2 = "" # Prediction Sequence -10 TGATGA
        s = "" # User Added dataset -35 and -10 values
        out = "" # User Added dataset outputs with & between each 0.5
        out2 = "" # Copy of out
        userData = [] # User added dataset in array form with -35 and -10 in alternative sequence

        # TEMPORARY VARIABLES - USED TO RETRIVE VAROUS PORTIONS OF THE STRING SENT FROM CLIENT
        sx = 0 # index of *
        sy = 0 # index of #
        j = 6 # index of *
        k = 0
        n = 0
        temp1 = 0
        temp3 = 1.1
        find = 0
        temp2 = ""
        outElements = []
        tempo = []
        convert10 = []
        convert35 = []
        flag = 0 # Used to check if the user has entered the sequence for which prediction is required
        flag2 = 0

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

        # STRENGTH OF -35 AND -10 SEQUENCE  - COPY OF outputResult
        outputResult2 = [
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

        # USED TO SLICE THE STRING TO OBTAIN
        # -35 AND -10 SEQUENCE FOR WHICH PREDICTION IS REQUIRED
        # THE DYNAMIC DATASET ENTERED BY THE USER - IT'S -35 AND -10 SEQUENCE ALONG WITH IT'S STRENGTH
        sy = t.index('#')
        sx = t.index('*')
        p = str(t[0:sx])
        if len(p) == 12:
            s1 = str(t[0:6])
            s2 = str(t[6:12])
            flag = 1
        s = str(t[sx+1:sy])
        out = str(t[sy+1:])
        out2 = str(t[sy+1:])
        while k < len(s) :
            userData.append(str(s[k:j]))
            k += 6
            j += 6
            numberOfElements += 1

        print "\t\t\t\t STRING SPLIT"
        print t
        print ""
        print s1
        print ""
        print s2
        print ""
        print s
        print ""
        print out
        print ""
        print userData
        print ""
        print numberOfElements
        print ""

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

        print "\t\t\t\t -35 SEQUENCE"
        print ""

        # TO APPEND THE -35 SEQUENCE DYNAMICALLY ADDED BY THE USER TO THE EXISTING DATASET
        i = 0
        while i < numberOfElements:
            if((i%2) == 0):
                instances.append(Seq(userData[i]))
            i += 1

        # CREATING A COPY OF INSTANCES ADDED WITH THE -35 SEQUENCE FOR WHICH WE HAVE TO PREDICT THE STRENGTH
        if flag == 1:
            #instances.append(Seq(s1))
            instancesP = instances[:]
            instancesP.append(Seq(s1))

        # CONVERTING THE -35 SEQUENCE INTO A SUITABLE FORMAT
        i = 0
        while i < len(instances):
            convert35.append(str(instances[i]))
            i+=1

        # CONSTRUCTION OF THE PSSM MATRIX FOR THE -35 SEQUENCE BASED ON STRENGTH
        m = motifs.create(instances_regdb)
        #print(m.counts);
        pwm = m.counts.normalize(pseudocounts= {'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49}                )
        #print(pwm)
        pssm = pwm.log_odds()
        #print(pssm)

        #REGRESSION MODELLING
        if flag == 1:
            p,o,l,k,m,n = str(s1)
            resultP = pssm[p,0] + pssm[o,1] + pssm[l,2] + pssm[k,3] + pssm[m,4] + pssm[n,5]

        result = []

        def calculateX(a,b,c,d,e,f,x):
            temp1 = pssm[a,0] + pssm[b,1] + pssm[c,2] + pssm[d,3] + pssm[e,4] + pssm[f,5]
            result.append([temp1])

        i = 0

        while i < len(convert35):
            calculateX(convert35[i][0],convert35[i][1],convert35[i][2],convert35[i][3],convert35[i][4],convert35[i][5],i)
            i +=1

        # EXTRACT THE OUTPUT FROM THE DYNAMICALLY ENTERED USER DATASET AND STORING IT IN OUTPUTRESULTS
        i = 0
        counter = 0
        half = numberOfElements*0.5
        while counter < half:
            find = out.index('&')
            temp2 = str(out[i:find])
            temp3 = float(temp2)
            outputResult.append([temp3])
            out = out[find+1:]
            counter +=1

        print "\t\t\t\t Obtaining input and output values for -35"
        print ""
        print "\t\t\t\t X1 Values -35 Sequence"
        print result
        print ""
        print "\t\t\t\t Y Values Strength"
        print outputResult
        print ""

        print "\t\t\t\t -10"

        # TO APPEND THE -10 SEQUENCE DYNAMICALLY ADDED BY THE USER TO THE EXISTING DATASET
        i = 0
        while i < numberOfElements:
            if((i%2) != 0):
                instances2.append(Seq(userData[i]))
            i += 1

        # CREATING A COPY OF INSTANCES ADDED WITH THE -10 SEQUENCE FOR WHICH WE HAVE TO PREDICT THE STRENGTH
        if flag == 1:
            #instances2.append(Seq(s2))
            instancesP2 = instances2[:]
            instancesP2.append(Seq(s2))

        # CONVERTING THE -10 SEQUENCE INTO A SUITABLE FORMAT
        i = 0
        while i < len(instances2):
            convert10.append(str(instances2[i]))
            i+=1

        # CONSTRUCTION OF THE PSSM MATRIX FOR THE -10 SEQUENCE BASED ON STRENGTH
        m2 = motifs.create(instances2_regdb)
        #print(m2.counts);
        pwm2 = m2.counts.normalize(pseudocounts={'A':0.49, 'C': 0.51, 'G' : 0.51, 'T' : 0.49})
        #print(pwm2)
        pssm2 = pwm2.log_odds()
        #print(pssm2)

        #REGRESSION MODELLING
        if flag == 1:
            p2,o2,l2,k2,m2,n2 = str(s2)
            resultP2 = pssm2[p2,0] + pssm2[o2,1] + pssm2[l2,2] + pssm2[k2,3] + pssm2[m2,4] + pssm2[n2,5]

        result2 = []

        def calculateX2(a,b,c,d,e,f,x):
            temp1 = pssm2[a,0] + pssm2[b,1] + pssm2[c,2] + pssm2[d,3] + pssm2[e,4] + pssm2[f,5]
            result2.append([temp1])

        i = 0
        while i < len(convert10):
            calculateX2(convert10[i][0],convert10[i][1],convert10[i][2],convert10[i][3],convert10[i][4],convert10[i][5],i)
            i +=1

        # EXTRACT THE OUTPUT FROM THE DYNAMICALLY ENTERED USER DATASET AND STORING IT IN OUTPUTRESULTS2 (THIS IS A COPY OF OUTPUTRESULTS)
        i = 0
        counter = 0
        half = numberOfElements*0.5
        while counter < half:
            find = out2.index('&')
            temp2 = str(out2[i:find])
            temp3 = float(temp2)
            outputResult2.append([temp3])
            out2 = out2[find+1:]
            counter +=1

        print "\t\t\t\t Obtaining input and output values for -10"
        print ""
        print "\t\t\t\t X2 Values -10 Sequence"
        print result2
        print ""
        print "\t\t\t\t Y Values Strength"
        print outputResult2

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
        while i<len(outputResult):
              b += outputResult[i]
              # TO DEAL WITH LOG(0)
              if b[i] == 0:
                  b[i] = 0.01 # change in cutoff
              b[i] = math.log(b[i])
              i +=1
        print ""
        print "\t\t\t\t Matrix B (Strength)"
        for x in b:
            print x
            print ""

        # FORMATTING MATRIX A AND B TO OBTAIN MATRIX x AND y
        print "\t\t\t\t Matrix X (Input Matrix : -35 and -10 Sequence)"
        x = np.asarray(a)
        print x
        print ""
        print "\t\t\t\t Matrix Y (Output Matrix : Strength)"
        y = np.asarray(b)
        print y
        print ""

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

        print "\t\t\t\t R Sqare Value"
        meany = np.mean(y)
        sumsqmeany = np.sum((y-meany)**2)
        sumsqmeanysum = np.sum((y-hx)**2)/sumsqmeany
        R = 1 - sumsqmeanysum
        print R
        print ""
        print "\t\t\t\t Adj. R Sqare Value"
        adjR = 1 - (1-R)*(len(instances)-1)/(len(instances)-2-1)
        print adjR
        print ""        

        # CONSTRUCTION OF THE MULTIVARIENT LINEAR REGRESSION GRAPH
        # THE NUMBER 18 REPRESENTS THE DEFAULT SIZE OF THE DATASET PROVIDED
        fig =plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        for c, m in [('r','o')]:
            xs = x[0:18,1]
            print "\t\t\t\t xs Plot for Graph -10 Sequence"
            print(xs)
            print ""
            ys = x[0:18,2]
            print "\t\t\t\t ys Plot for Graph -35 Sequence"
            print (ys)
            print ""
            zs = y[0:18]
            print "\t\t\t\t zs Plot for Graph Strength"
            print (zs)
            print ""
            ax.scatter(xs, ys,zs, c=c, marker =m)

        md = len(y)
        print "\t\t\t\t Total Number of Elements in the Dataset"
        print md
        print ""
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

        # TO SAVE THE GRAPH AS AN IMAGE
        fig1 = plt.gcf()
        # To display image in console : plt.show()
        fig1.savefig('Multivariant.png', format='png', dpi=300)
        data_uri = open('Multivariant.png', 'rb').read().encode('base64').replace('\n', '')
        graph = '<div class="row"><div class="col s12 m12 l12"><h2><center><i class="fa fa-line-chart fa-2x" aria-hidden="true"></i>&nbsp;&nbsp;Regression Model And Predicted Output</center></h2></div></div>'
        img_tag = '<div class="row"><div class="col s12 m7 l7"><img id="reg" src="data:image/png;base64,%s"></div>' % data_uri
        #print theta, J_history
        plot(np.arange(numIterations), J_history)
        xlabel('Iterations')
        ylabel('Cost Function')
        #show()

        # TO SEND THE RESULTS BACK TO THE CLIENT
        if flag == 1:
            strength = np.array([1.0,   resultP, resultP2 ]).dot(theta)
            finalStrength = strength
            print "\t\t\t\t Predicted Strength",
            print math.exp(strength)
            print ""
            print "\t\t\t\t Predicted ln(Strength)",
            print strength
            print ""            
            #print 'Predicted strength of promoter : %s' % (finalStrength)
            #print 'Correlation Measure Sqare Value : %s' % (R)
            pri = '<div class="col s12 m5 l5"><p class="para"><span class="highlight"><i class="fa fa-bullseye fa-2x" aria-hidden="true"></i></span><br><span class="highlight">Predicted ln(Strength)</span><br> %s </p>' % finalStrength
            key1 = '<p class="para"><span class="highlight redDot"><i class="fa fa-circle" aria-hidden="true"></i></span><br><span class="highlight">Existing Dataset</span></p>'
            key2 = '<p class="para"><span class="highlight blueDot"><i class="fa fa-circle" aria-hidden="true"></i></span><br><span class="highlight">User Added Data Points</span> <br> (If Dynamic Model in Use)</p>'
            rsqare = '<p class="para"><span class="highlight"><i class="fa fa-gavel fa-2x" aria-hidden="true"></i></span><br><span class="highlight">Correlation Measure Sqare Value</span> <br> %s </p></div></div>' % R
        else:
            pri = '<div class="col s12 m5 l5"><p></p>'
            key1 = '<p class="para"><span class="highlight redDot"><i class="fa fa-circle" aria-hidden="true"></i></span><br><span class="highlight">Existing Dataset</span></p>'
            key2 = '<p class="para"><span class="highlight blueDot"><i class="fa fa-circle" aria-hidden="true"></i></span><br><span class="highlight">User Added Data Points</span> <br> (If Dynamic Model in Use)</p>'
            rsqare = '<p class="para"><span class="highlight"><i class="fa fa-gavel fa-2x" aria-hidden="true"></i></span><br><span class="highlight">Correlation Measure Sqare Value</span> <br> %s </p></div></div>' % R

        return make_text(graph + img_tag + pri + key1 + key2 + rsqare)

if __name__ == '__main__':
    app.run()

application = app.wsgifunc()
