import tensorflow

def aratio(r, a):
    return tensorflow.math.logical_and(r/a >= 2.5, r/a <= 3.5)

def r_limits(r):
    return tensorflow.math.logical_and(r >= 3.3, r <= 4.7)

def pinj_limits(pinj):
    return tensorflow.math.logical_and(pinj >= 20.0, pinj <= 50.0)

def fbs_limits(fbs):
    return tensorflow.math.logical_and(fbs >= 0.0, fbs <= 1.0)

def pfus_limits(pfus):
    return tensorflow.math.logical_and(pfus >= 0.0, pfus <= 1000.0)

manifest = [
    {'function' : aratio,      'args' : ['r','a'] },
    {'function' : r_limits,    'args' : ['r']     },
    {'function' : pinj_limits, 'args' : ['pinj']  },
    {'function' : fbs_limits,  'args' : ['fbs']   },
    {'function' : pfus_limits, 'args' : ['pfus']  }
]
