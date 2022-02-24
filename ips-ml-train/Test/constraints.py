import tensorflow

def fbs(ip, ibs):
    return ibs/ip < 1.0

def aratio(r, a):
    return (2.33 <= r/a) & (r/a <= 3.5)

def q95_input(r, a, bt, ip, kappa, delta):
    eps = a/r
    q95 = 5.0*(a*eps*bt/ip)*((1.0+kappa*kappa*(1.0+2.0*delta*delta-1.2*delta*delta*delta))/2.0)*((1.17-0.65*eps)/(tensorflow.math.pow(1.0-tensorflow.math.pow(eps,2.5),2)))
    return (q95 >= 4.0) & (q95 <= 7.0)

def q95(q95_value):
    return (q95_value >= 4.0) & (q95_value <= 7.0)

def greenwald(ip,a,ne):
    limit = ip/(math.pi*a*a)
    return ne < limit

def r_limits(r):
    return (r >= 1.4) & (r <= 2.2)

def a_limits(a):
    return (a >= 0.4) & (a <= 0.94)

def kappa_limits(kappa):
    return (kappa >= 1.7) & (kappa <= 2.0)

def delta_limits(delta):
    return (delta >= 0.5) & (delta <= 0.65)

def bt_limits(bt):
    return (bt >= 3.0) & (bt <= 7.0)

def ip_limits(ip):
    return (ip >= 0.74) & (ip <= 12.0)

def fgw_ped_limits(fgw_ped):
    return (fgw_ped >= 0.2) & (fgw_ped <= 1.0)

def nepeak_limits(nepeak):
    return (nepeak >= 1.5) & (nepeak <= 2.0)

def betan_limits(betan):
    return (betan >= 2.5) & (betan <= 4.5)

def betan_ped_limits(betan_ped):
    return (betan_ped >= 0.5) & (betan_ped <= 1.2)

def ibs_limits(ibs):
    return (ibs >= 0.74) & (ibs <= 12.0)

manifest = [
    {'function' : fbs,              'args' : ['ip','ibs']                        },
    {'function' : aratio,           'args' : ['r','a']                           },
    {'function' : q95_input,        'args' : ['r','a','bt','ip','kappa','delta'] },
    {'function' : q95,              'args' : ['q95_input']                       },
    {'function' : greenwald,        'args' : ['ip','a','neped']                  },
    {'function' : r_limits,         'args' : ['r']                               },
    {'function' : a_limits,         'args' : ['a']                               },
    {'function' : kappa_limits,     'args' : ['kappa']                           },
    {'function' : delta_limits,     'args' : ['delta']                           },
    {'function' : bt_limits,        'args' : ['bt']                              },
    {'function' : ip_limits,        'args' : ['ip']                              },
    {'function' : fgw_ped_limits,   'args' : ['fgw_ped']                         },
    {'function' : nepeak_limits,    'args' : ['nepeak']                          },
    {'function' : betan_limits,     'args' : ['betan']                           },
    {'function' : betan_ped_limits, 'args' : ['betan_ped']                       },
    {'function' : ibs_limits,       'args' : ['ibs']                             }
]
