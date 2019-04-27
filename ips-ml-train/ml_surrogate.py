# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 23:21:53 2019

@author: Dr Clement Etienam
"""
from keras.layers import Dense, Activation
from keras.models import Sequential
import numpy as np
def surrogateforward(inputini,outputini,max_epocs):
    print ('deep neural netowrk using keras for the tauth data')
    np.random.seed(7)
    modelDNN = Sequential()

# Adding the input layer and the first hidden layer
    modelDNN.add(Dense(200, activation = 'relu', input_dim = 9))

# Adding the second hidden layer
    modelDNN.add(Dense(units = 420, activation = 'relu'))

# Adding the third hidden layer
    modelDNN.add(Dense(units = 21, activation = 'relu'))

# Adding the output layer

    modelDNN.add(Dense(units = 1))

    modelDNN.compile(optimizer = 'adam', loss = 'mean_squared_error')

# Fitting the ANN to the Training set
    modelDNN.fit(inputini, outputini, batch_size = 100, epochs = max_epocs)

    
    return modelDNN
