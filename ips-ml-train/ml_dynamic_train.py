# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:34:45 2019

@author: Dr Clement Etienam
Supevisor :Professor Kody Law
Collaborator: Dr Mark Cianciosa
Train CCR with an Active learning aproach on the Tauth data
-------------------------------------------------------------------------------
This is the Logic
1)We wil first detrmine the optimum number of clusters suitable for the data
2)Do K-means clustering of the input-output pair
3)Train a classifier with the input and the labels from step (2)
4)Sample those points with low probability
5) Generate more points around these uncertain points
6) Get the coresponding simulated (y) data of these generated points
7) Add these points and retrain the machine
-------------------------------------------------------------------------------
"""
#Import the necessary libraries
from __future__ import print_function
print(__doc__)

from sklearn.neural_network import MLPClassifier
import numpy as np
from sklearn.cluster import KMeans

from numpy import linalg as LA

from sklearn.utils import check_random_state

from sklearn.ensemble import RandomForestClassifier
import pickle
from sklearn.metrics import classification_report, confusion_matrix
import pandas as pd
from ml_forwarding import ensembleforwarding
from ml_surrogate import surrogateforward
from ml_machine_model import passiveclement
from sklearn.model_selection import train_test_split

import sys, getopt
opts, args = getopt.getopt(sys.argv[1:], 'c:i:n:q:e:f:x:y:m:h', ['max_clusters=', 'max_iterations=', 'num_ensembles=','max_queried=','epocs=','input_file=','x_data=','y_data=','model_save=','help'])

input_keys=[];
output_keys=[];
commandline_arguments = {'x_data': input_keys, 'y_data': output_keys}

for opt, arg in opts:
    if opt in ('-c','--max_clusters'):
        commandline_arguments['max_clusters'] = np.int(arg)
    if opt in ('-i','--max_iterations'):
        commandline_arguments['max_iterations'] = np.int(arg)
    if opt in ('-n','--num_ensembles'):
        commandline_arguments['num_ensembles'] = np.int(arg)
    if opt in ('-q','--max_queried'):
        commandline_arguments['max_queried'] = np.int(arg)
    if opt in ('-e','--epocs'):
        commandline_arguments['epocs'] = np.int(arg)
    if opt in ('-f','--input_file'):
        commandline_arguments['input_file'] = arg
    if opt in ('-x','--x_data'):
        commandline_arguments['x_data'].append(arg)
    if opt in ('-y','--y_data'):
        commandline_arguments['y_data'].append(arg)
    if opt in ('-m','--model_save'):
        commandline_arguments['model_save'] = arg
    if opt in ('-h','--help'):
        print('usage: ml_dynamic_train.py -[cieh] {option} --[args]={option}')
        print('  -h, --help           Print this help message.')
        print('  -o, --max_clusters   Maximum number of clusters you want to accomodate.')
        print('  -i, --max_iterations Maximum number of iteration you want to run.')
        print('  -e, --num_ensembles  Number of ensemble you want for pertubation.')
        print('  -q, --max_queried    Maximum number of queries to try.')
        print('  -e, --epocs          Number of epocs to use.')
        print('  -f, --input_file     Initial input data file.')
        print('  -x, --x_data         Training input data. Allows multiple entries.')
        print('  -y, --y_data         Training output data. Allows multiple entries.')
        exit(0);

max_queried = 1000
if 'max_queried' in commandline_arguments.keys():
    max_queried = commandline_arguments['max_queried']

max_epocs = 1000
if 'epocs' in commandline_arguments.keys():
    max_epocs = commandline_arguments['epocs']

max_clusters = 0
if 'max_clusters' in commandline_arguments.keys():
    max_clusters = commandline_arguments['max_clusters']
else:
    max_clusters = np.int(input("Enter the maximum number of clusters you want to accomodate: "))

iterr = 1
if 'max_iterations' in commandline_arguments.keys():
    iterr = commandline_arguments['max_iterations']
else:
    iterr = np.int(input("Enter the maximum number of iteration you want to run: "))

Ne = 0
if 'num_ensembles' in commandline_arguments.keys():
    Ne = commandline_arguments['num_ensembles']
else:
    Ne= np.int(input("Enter the number of ensemble you want for pertubation: "))

trainset_size=5
mark = (max_queried/trainset_size) - 1

print('Determine he optimum number of clusers')


def optimalK(data, nrefs, maxClusters):
    """
    Calculates KMeans optimal K using Gap Statistic 
    """
    gaps = np.zeros((max(1, maxClusters - 1),))
    resultsdf = pd.DataFrame({'clusterCount':[], 'gap':[]})
    for gap_index, k in enumerate(range(1, maxClusters)):

# Holder for reference dispersion results
        refDisps = np.zeros(nrefs)

# For n references, generate random sample and perform kmeans getting resulting dispersion of each loop
        for i in range(nrefs):
            
# Create new random reference set
            randomReference = np.random.random_sample(size=data.shape)
            
# Fit to it
            km = KMeans(k)
            km.fit(randomReference)
            
            refDisp = km.inertia_
            refDisps[i] = refDisp

# Fit cluster to original data and create dispersion
        km = KMeans(k)
        km.fit(data)
        
        origDisp = km.inertia_
# Calculate gap statistic
        gap = np.log(np.mean(refDisps)) - np.log(origDisp)

# Assign this loop's gap statistic to gaps
        gaps[gap_index] = gap
        
        resultsdf = resultsdf.append({'clusterCount':k, 'gap':gap}, ignore_index=True)

    return (gaps.argmax() + 1, resultsdf)  # Plus 1 because index of 0 means 1 cluster is optimal, index 2 = 3 clusters are optimal

def run_model(model,X_labeled,y_labeled,X_unlabeled,y_oracle,X_test,y_test):
    from sklearn.metrics import accuracy_score
    # build the model on training data
    random_state = check_random_state(0)
    initial_labeled_samples = 5
    permutation = np.random.choice(trainset_size,initial_labeled_samples,replace=False)
    X_train = X_labeled[permutation]
    y_train = y_labeled[permutation]
    X_train = X_train.reshape((X_train.shape[0], -1))
    queried = initial_labeled_samples
    
    X_valclement = np.copy(X_unlabeled)
    X_val = np.copy(X_unlabeled)
    
    X_val = np.delete(X_val, permutation, axis=0)
    X_valclement = np.delete(X_valclement, permutation, axis=0)
    
    y_val = np.copy(y_oracle)
    y_val = np.delete(y_val, permutation, axis=0)
    
    scaler = MinMaxScaler()
    
    X_train = scaler.fit_transform(X_train)
    X_val   = scaler.transform(X_val)
    X_test  = scaler.transform(X_test)
    X_train = scaler.inverse_transform(X_train)
    X_val   = scaler.inverse_transform(X_val)
    X_test  = scaler.inverse_transform(X_test)
    
    model.fit(X_train, y_train)
    ff_array = np.array([])
    queried_array2 = np.array([])
   
    queried_array = np.empty((numrows, 0))
    queried_clement = np.empty((trainset_size*numcols, 0))
    
    while queried < max_queried:
        probas_val = model.predict_proba(X_val)
        rev = np.sort(probas_val, axis=1)[:, ::-1]
        values = rev[:, 0] - rev[:, 1]
        selection = np.argsort(values)[:initial_labeled_samples]
        uncertain_samples=selection

        X_train = scaler.inverse_transform(X_train)
        X_val   = scaler.inverse_transform(X_val)
        X_test  = scaler.inverse_transform(X_test)
        X_train = scaler.fit_transform(X_train)
        X_val   = scaler.transform(X_val)
        X_test  = scaler.transform(X_test)

        X_train = np.concatenate((X_train, X_val[uncertain_samples]))
        y_train = np.concatenate((y_train, y_val[uncertain_samples]))

        X_val = np.delete(X_val, uncertain_samples, axis=0)
        y_val = np.delete(y_val, uncertain_samples, axis=0)


        X_train = scaler.inverse_transform(X_train)
        X_val   = scaler.inverse_transform(X_val)
        X_test  = scaler.inverse_transform(X_test)
        X_train = scaler.fit_transform(X_train)
        X_val   = scaler.transform(X_val)
        X_test  = scaler.transform(X_test)
        queried += initial_labeled_samples

# make predictions for test data
        model.fit(X_train, y_train)
        pickle.dump(model, open(commandline_arguments['model_save'], 'wb')) #Save it
        labelDA = model.predict(X_test)
        cm = confusion_matrix(y_test, labelDA,
                              labels=model.classes_)
        labelDAA = np.reshape(labelDA,(-1,1))
        queried_array = np.append(queried_array, labelDAA, axis=1)
        uncertainmark = X_valclement[uncertain_samples]
        queried_clement = np.append(queried_clement, np.reshape(uncertainmark,(-1,1),'F'), axis=1)
            
        ff = accuracy_score(y_test, labelDA)*100
        ff_array = np.append(ff_array, ff)
        queried_array2 = np.append(queried_array2, queried)
            
# Uncertainclement are the points we dont know
# label_array = np.append(label_array, labelDA)
        print('The accuracy is',ff)
        print('Finished querying',queried,'points')
        print("Confusion matrix after query",queried)
        print(cm)

        if ff == 100.0:
            break
    
    return labelDAA,ff_array,queried_array,queried_array2,queried_clement


print('Load the Initial Model data')

import json
with open(commandline_arguments['input_file'],'r') as json_ref:
    test = json.load(json_ref)

#  Create input data.
input = np.array(test[input_keys[0]])
for i in range(1,len(input_keys)):
    input = np.concatenate((input, np.array(test[input_keys[i]])), axis=1)

#  Create output data.
output = np.array(test[output_keys[0]])
for i in range(1,len(output_keys)):
    output = np.concatenate((input, np.array(test[output_keys[i]])), axis=1)

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()

# Randomly split the data in half.
permutation = np.random.permutation(input.shape[0])
half = permutation.shape[0]/2

# Split input data.
inputtrain = np.log(input[permutation[:half],:]) #select the first half for training
inputini = inputtrain
inputtest = np.log(input[permutation[half:],:]) #use the remaining data for testing

# Split output data.
outputtrain = np.log(output[permutation[:half],:]) #select the first half for training
outputini = outputtrain
outputtest = (output[permutation[half:],:]) #use the remaining data for testing

numrowstest = len(outputtest)
modelDNN = surrogateforward(inputini,outputini, max_epocs) #The DN surrogate model
for iclement in range (iterr):
    print('Starting Iteration {}'.format(iclement))

    numrows = len(inputtrain) # rows of inout
    numcols = len(input[0]) # columns of input
#Start the loop to dynamically enrich the training set for sparse data coverage
  
# model = MLPClassifier(solver= 'lbfgs',max_iter=5000)
    model = RandomForestClassifier(n_estimators=500)
    matrix = np.concatenate((inputtrain,outputtrain), axis=1)
    matrixtest = np.concatenate((inputtest,outputtest), axis=1)
# use Gap statistics to get optimum number of clustsrs
    nclusters, gapdf = optimalK(matrix, nrefs=5, maxClusters=max_clusters)
    print ('Optimal k is: ', nclusters)

# Now do the K-means clustering
    print('Do the K-means clustering of [X,y] and get the labels')
    kmeans = KMeans(n_clusters=nclusters,max_iter=100).fit(matrix)
    dd = np.reshape(kmeans.labels_,(numrows,1))
    #Xini = inputtrain
    
    print('Do for the test data as well')
    
    X_labeled, X_unlabeled, y_labeled, y_oracle = train_test_split(inputtrain, dd, test_size=0.99)
    print('Enter the active learning dynamics')
    labelDAA,ff_array,queried_array,queried_array2,queried_clement = run_model(model,X_labeled,y_labeled,X_unlabeled,y_oracle,inputtrain,dd)
    unie = np.amax(ff_array)
    ff_array = np.reshape(ff_array,(-1,1))

    label0 = (np.asarray(np.where(ff_array==unie))).T
    finallabel = queried_array[:,label0[0,0]]

    labelDA = queried_array[:,-1]
    print('The highest accuracy is',unie,'at query point',label0[0,0])
    print('generate a distribution about this uncertain points') 
# Routine to run a forward code of the uncertain points
    (xensemble) = ensembleforwarding(queried_clement, Ne)
    print('Save the data')
    np.savetxt('newdata.out', np.reshape(xensemble,(-1,1),'F'), fmt = '%4.6f', newline = '\n')
# Augment this points
# Pause here and use the newdata.out to get new points, (reshape the newdata
# To be np.reshape(newdata(-1,9)))
    xensemble=np.reshape(xensemble,(trainset_size,numcols,mark))
    for i in range(mark):
        inputtrain = np.concatenate((inputtrain, xensemble[:,:,i]))
    
# Get yensemble from the forwarding and reshape doing
# yensemble=np.reshape(yensemble,(-1,1)'F')
# y_keepit = np.concatenate((y_keepit, yensemble))
    np.savetxt('newtrainingset.out', np.exp(inputtrain), fmt = '%4.6f', newline = '\n')

    print('For now use a DNN surrogate model for forwarding')
    
    outputtrain = modelDNN.predict(inputtrain)
    print('Finished Iteration %d'%iclement)
    
print('now some predictions')
(clementanswer) = passiveclement(inputtrain, outputtrain, inputtest, outputtest)
print(' Compute L2 and R2 for the machine')

outputtest = np.reshape(outputtest, (numrowstest, 1))
Lerrorsparse = (LA.norm(outputtest - clementanswer)/LA.norm(outputtest))**0.5
L_2sparse = 1 - Lerrorsparse**2

#Coefficient of determination
outputreq = np.zeros((numrowstest,1))
for i in range(numrowstest):
    outputreq[i,:] = outputtest[i,:] - np.mean(outputtest)


#outputreq=outputreq.T
CoDspa = 1 - (LA.norm(outputtest - clementanswer)/LA.norm(outputreq))
CoDsparse = 1 - (1-CoDspa)**2 ;
print ('R2 of fit using the machine is :', CoDsparse)
print ('L2 of fit using the machine is :', L_2sparse)

print('program executed')
