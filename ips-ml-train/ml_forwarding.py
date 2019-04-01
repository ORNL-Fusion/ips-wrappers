##------------------------------------------------------------------------------
##'Forwarding code'  
##'Author: Dr Clement Etienam'
##'Supervisor: Professor Kody Law'
##'Collaborator: Dr Mark Cianciosa'

##------------------------------------------------------------------------------

import numpy as np

def ensembleforwarding(queried_clement,Ne):
    X_keepituncertain=queried_clement
   
    print('Determine how many uncertain samples came thorugh')
    rowsclement=len(X_keepituncertain)
    colsclement = len(X_keepituncertain[0]) # columns of input
    grand= np.ones((45,colsclement,colsclement))
    sgsimuse= np.ones((45,colsclement)) 
    for j in range(colsclement):
        Xuse=X_keepituncertain[:,j]
        for i in range(colsclement):
            outt=np.random.normal(0,0.5,(45))
            sumit=np.reshape(outt,(1,-1),'F')+Xuse
            sgsimuse[:,i] = sumit
        grand[:,:,j]=sgsimuse
    xensemble=grand[:,:,0]
    #Cal the forwarding model here
    print('  Programme executed  ')

    return( xensemble )
 

