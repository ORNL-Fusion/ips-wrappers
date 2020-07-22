def fbs(fbs_value):
    return fbs_value < 1.0

def aratio(r, a):
    return 2.33 <= r/a & r/a <= 3.5

manifest = [
    {'function' : fbs,    'args' : ['fbs']  },
    {'function' : aratio, 'args' : ['r','a']}
]
