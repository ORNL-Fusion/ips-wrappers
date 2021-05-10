import tensorflow

def fbs(fbs_value):
    return fbs_value < 1.0

def aratio(r, a):
    return (2.33 <= r/a) & (r/a <= 3.5)

def q95_input(r, a, bt, ip, kappa, delta):
    eps = a/r
    q95 = 5.0*(a*eps*bt/ip)*((1.0+kappa*kappa*(1.0+2.0*delta*delta-1.2*delta*delta*delta))/2.0)*((1.17-0.65*eps)/(tensorflow.math.pow(1.0-tensorflow.math.pow(eps,2.5),2)))
    return (q95 >= 4.0) & (q95 <= 7.0)

def q95(q95_value):
    return (q95_value >= 4.0) & (q95_value <= 7.0)

manifest = [
    {'function' : fbs,       'args' : ['fbs']  },
    {'function' : aratio,    'args' : ['r','a']},
    {'function' : q95_input, 'args' : ['r','a','bt','ip','kappa','delta']},
    {'function' : q95,       'args' : ['q95_input']}
]
