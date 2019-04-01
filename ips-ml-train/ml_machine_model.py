# -*- coding: utf-8 -*-
"""
Created on Tuesday Feb 05 12:05:47 2019

@author: Dr Clement Etienam
This is the code for learning a machine for FAST-TRAN
We will cluster th data first, use that label from the cluster and learn a
classifier then a regressor
This code is very important for Chi-data
"""
from __future__ import print_function
print(__doc__)

from sklearn.neural_network import MLPClassifier
import numpy as np

from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from datetime import datetime
from numpy import linalg as LA
from sklearn.preprocessing import MinMaxScaler
#------------------Begin Code----------------#
print('Load the Tauth/Fastran data')
#max_clusters = np.int(input("Enter the maximum number of clusters you want to accomodate: ") )
#iterr = np.int(input("Enter the maximum number of iteration you want to run: ") )
#Ne= np.int(input("Enter the number of ensemble you want for pertubation: ") )
import scipy.io as sio
def passiveclement(inputtrain,outputtrain,inputtest,outputtest):
   

    numrowstest=len(outputtest)
    matrix=np.concatenate((inputtrain,outputtrain), axis=1)
    print('Do the K-means clustering with 5 clusters of [X,y] and get the labels')
    kmeans = MiniBatchKMeans(n_clusters=5,random_state=0,batch_size=6,max_iter=10).fit(matrix)
    rowsclement=len(inputtrain)
    colsclement = len(inputtrain[0]) # columns of input
    dd=kmeans.labels_
    dd=dd.T
    dd=np.reshape(dd,(rowsclement,1))
#-------------------#---------------------------------#
    print('Use the labels to train a classifier')
    inputtrainclass=inputtrain
    outputtrainclass=np.reshape(dd,(rowsclement,1))
    ydami=outputtrain
#%outputtest=y(290000+1:end,:);
    def run_model(model):
    # build the model on training data
        model.fit(inputtrainclass, outputtrainclass )

    # make predictions for test data
        labelDA = model.predict(inputtest)
        return labelDA
    print(' Learn the classifer from the predicted labels from Kmeans')
    model = MLPClassifier(solver= 'lbfgs',max_iter=3000)
    print('Predict the classes from the classifier for test data')
    labelDA=run_model(model)
#-------------------Regression----------------#
    print('Learn regression of the 5 clusters with different labels from k-means ' )
    label0=(np.asarray(np.where(dd == 0))).T
    label1=(np.asarray(np.where(dd == 1))).T
    label2=(np.asarray(np.where(dd == 2))).T
    label3=(np.asarray(np.where(dd == 3))).T
    label4=(np.asarray(np.where(dd == 4))).T

    print('set the output matrix')
    clementanswer=np.zeros((numrowstest,1))

    print('Start the regression')
    from sklearn.neural_network import MLPRegressor
##
    model0 = MLPRegressor(solver= 'lbfgs',max_iter=3000)
    a0=inputtrainclass[label0,:]
    a0=np.reshape(a0,(-1,9),'F')

    b0=ydami[label0,:]
    b0=np.reshape(b0,(-1,1),'F')
    model0.fit(a0, b0)

##
    model1 = MLPRegressor(solver= 'lbfgs',max_iter=3000)
    a1=inputtrainclass[label1,:]
    a1=np.reshape(a1,(-1,9),'F')

    b1=ydami[label1,:]
    b1=np.reshape(b1,(-1,1),'F')
    model1.fit(a1, b1)

##
    model2 = MLPRegressor(solver= 'lbfgs',max_iter=3000)
    a2=inputtrainclass[label2,:]
    a2=np.reshape(a2,(-1,9),'F')

    b2=ydami[label2,:]
    b2=np.reshape(b2,(-1,1),'F')
    model2.fit(a2, b2)

##
    model3 = MLPRegressor(solver= 'lbfgs',max_iter=3000)
    a3=inputtrainclass[label3,:]
    a3=np.reshape(a3,(-1,9),'F')

    b3=ydami[label3,:]
    b3=np.reshape(b3,(-1,1),'F')
    model3.fit(a3, b3)

##
    model4 = MLPRegressor(solver= 'lbfgs',max_iter=3000)
    a4=inputtrainclass[label4,:]
    a4=np.reshape(a4,(-1,9),'F')

    b4=ydami[label4,:]
    b4=np.reshape(b4,(-1,1),'F')
    model4.fit(a4, b4)


    print('Time for the prediction')
    labelDA0=(np.asarray(np.where(labelDA == 0))).T
    labelDA1=(np.asarray(np.where(labelDA == 1))).T
    labelDA2=(np.asarray(np.where(labelDA == 2))).T
    labelDA3=(np.asarray(np.where(labelDA == 3))).T
    labelDA4=(np.asarray(np.where(labelDA == 4))).T


    a00=inputtest[labelDA0,:]
    a00=np.reshape(a00,(-1,9),'F')
    if a00.shape[0] != 0:
        clementanswer[np.ravel(labelDA0),:]=np.exp(np.reshape(model0.predict(a00),(-1,1)))


    a11=inputtest[labelDA1,:]
    a11=np.reshape(a11,(-1,9),'F')
    if a11.shape[0] != 0:
        clementanswer[np.ravel(labelDA1),:]=np.exp(np.reshape(model1.predict(a11),(-1,1)))

    a22=inputtest[labelDA2,:]
    a22=np.reshape(a22,(-1,9),'F')
    if a22.shape[0] != 0:
        clementanswer[np.ravel(labelDA2),:]=np.exp(np.reshape(model2.predict(a22),(-1,1)))

    a33=inputtest[labelDA3,:]
    a33=np.reshape(a33,(-1,9),'F')
    if a33.shape[0] != 0:
        clementanswer[np.ravel(labelDA3),:]=np.exp(np.reshape(model3.predict(a33),(-1,1)))

    a44=inputtest[labelDA4,:]
    a44=np.reshape(a44,(-1,9),'F')
    if a44.shape[0] != 0:
        clementanswer[np.ravel(labelDA4),:]=np.exp(np.reshape(model4.predict(a44),(-1,1)))
    return clementanswer




