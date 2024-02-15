#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import h5py
import shutil

def transferGrid(networkFile, print_test=False):

	## Copy the old network file to the new network file
	shutil.copyfile(networkFile, networkFile.split('.')[0] + '_oldGrid.h5')

	## The number of grid points to be added
	nPoints = 25

	## Open the file
	f0 = h5py.File(networkFile, 'r+')
	## Read the grid to know which grid point is which depth
	gridDset0 = f0['headerGroup/grid']
	oldGridSize = len(gridDset0)

	## Look for the original grid position
	for i in range(1,len(gridDset0)):
		if (gridDset0[i] - gridDset0[i-1] > 0.11): break
	initSurf = i - 25
	## Compute the current void portion
	oldVoidPortion = 100.0 * float(initSurf) / float(len(gridDset0))
	if print_test:
		print('\t \t Old void portion:', oldVoidPortion, 'with', oldGridSize, 'grid points.')
        
	## Compute the new one
	newGridSize = oldGridSize + nPoints
	if print_test:
		print('\t \t oldGridSize=', oldGridSize, '; newGridSize=', newGridSize)
	newVoidPortion = 100.0 * float(initSurf + nPoints) / float(newGridSize) #len(gridDset0) + nPoints)
	if print_test:
		print('\t \t New void portion:', newVoidPortion, 'with', newGridSize, ' grid points.')

	## Add points to the new grid if nPoints>0
	newGrid = []
	if (nPoints > 0):
		for i in range(nPoints):
			newGrid.append(i * 0.1)
		for i in range(len(gridDset0)):
			newGrid.append(gridDset0[i] + (nPoints) * 0.1)
	## Remove them
	else:
		for i in range(-nPoints, len(gridDset0)):
			newGrid.append(gridDset0[i] + (nPoints) * 0.1)

	## Replace the grid dataset
	gridArray = np.array(newGrid, dtype=np.float)
	del f0['headerGroup/grid']
	dataset = f0.create_dataset('headerGroup/grid', (len(gridArray),), dtype=np.float)
	dataset[...] = gridArray

	## Open the header group
	headerGroup = f0['headerGroup']
	## Replace the nx attribute
	nx = headerGroup.attrs['nx']
	headerGroup.attrs['nx'] = nx + nPoints

	## Get the last time step saved in the file
	concGroup0 = f0['concentrationsGroup']
	timestep = concGroup0.attrs['lastTimeStep']
	## Open the concentration group at this time step
	if print_test:
		print('time step is', timestep)
	if (timestep == -1):
		print('WARNING: network file was empty (timestep -1)')
		print('\t will return old values')
		return [oldGridSize, oldVoidPortion]
	#if timestep exists, continue:		
	groupName ='concentration_' + str(timestep)
	subConcGroup0 = concGroup0[groupName]
	## Read the concentration and index datasets
	concDset0 = subConcGroup0['concs']
	indexDset0 = subConcGroup0['concs_startingIndices']
	## Replace the surface position
	surfacePos0 = subConcGroup0.attrs['iSurface']
	subConcGroup0.attrs['iSurface'] = surfacePos0 + nPoints

	## Read the last data in the conc dataset, it corresponds to the temperature and its index
	tempData = concDset0[len(concDset0)-1]
	## Create the new concentrations
	newConc = []
	if (nPoints > 0):
		for i in range(nPoints):
			tempConc = (tempData[0], tempData[1])
			newConc.append(tempConc)
		for i in range(len(concDset0)):
			tempConc = (concDset0[i][0], concDset0[i][1])
			newConc.append(tempConc)
	else:
		for i in range(-nPoints, len(concDset0)):
			tempConc = (concDset0[i][0], concDset0[i][1])
			newConc.append(tempConc)

	## Create a specific data type for concentrations
	dType = np.dtype([('ConcType.first', np.int32),
					  ('ConcType.second', np.float)])
	concArray = np.array(newConc, dtype=dType)
	## Replace the concentration dataset
	del subConcGroup0['concs']
	dataset = subConcGroup0.create_dataset('concs', (len(concArray),), dtype=dType)
	dataset[...] = concArray
	## Create the new conc indices
	newIndex = []
	if (nPoints > 0):
		for i in range(nPoints):
			newIndex.append(i)
		for i in range(len(indexDset0)):
			newIndex.append(indexDset0[i] + nPoints)
	else:
		for i in range(-nPoints, len(indexDset0)):
			newIndex.append(indexDset0[i] + nPoints)
	indexArray = np.array(newIndex, dtype=np.uint32)
	del subConcGroup0['concs_startingIndices']
	dataset = subConcGroup0.create_dataset('concs_startingIndices', (len(indexArray),), dtype=np.uint32)
	dataset[...] = indexArray

	if print_test:
		print('from transferGrid, returning newGridSize, voidPortion = ', newGridSize, newVoidPortion)
	return [newGridSize, newVoidPortion]
