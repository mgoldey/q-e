import numpy as np

# read S
saaup = np.loadtxt('Saa1')
saadown = np.loadtxt('Saa2')

sabup = np.loadtxt('Sab1')
sabdown = np.loadtxt('Sab2')

sbaup = np.loadtxt('Sba1')
sbadown = np.loadtxt('Sba2')

sbbup = np.loadtxt('Sbb1')
sbbdown = np.loadtxt('Sbb2')

detsup = np.loadtxt('detS1')
detsdown = np.loadtxt('detS2')

# read cofac mat
caaup = np.loadtxt('Caa1')
caadown = np.loadtxt('Caa2')

cabup = np.loadtxt('Cab1')
cabdown = np.loadtxt('Cab2')

cbaup = np.loadtxt('Cba1')
cbadown = np.loadtxt('Cba2')

cbbup = np.loadtxt('Cbb1')
cbbdown = np.loadtxt('Cbb2')

# read wmat
waaup = np.loadtxt('Waa1')
waadown = np.loadtxt('Waa2')

wabup = np.loadtxt('Wab1')
wabdown = np.loadtxt('Wab2')

wbaup = np.loadtxt('Wba1')
wbadown = np.loadtxt('Wba2')

wbbup = np.loadtxt('Wbb1')
wbbdown = np.loadtxt('Wbb2')

wmatup = np.loadtxt('W1')
wmatdown = np.loadtxt('W2')

#
# compute S mat
#

# up
pysmatup = [[round(np.linalg.det(saaup), 8), round(np.linalg.det(sabup), 8)],
            [round(np.linalg.det(sbaup), 8), round(np.linalg.det(sbbup), 8)]]
print """== S up python vs QE ==
{}\t{}
{}\t{}
""".format(pysmatup[0], detsup[0], pysmatup[1], detsup[1])

# down
pysmatdown = [[round(np.linalg.det(saadown), 8),
               round(np.linalg.det(sabdown), 8)],
              [round(np.linalg.det(sbadown), 8),
               round(np.linalg.det(sbbdown), 8)]]
print """== S down python vs QE ==
{}\t{}
{}\t{}
""".format(pysmatdown[0], detsdown[0], pysmatdown[1], detsdown[1])
#
# compute cofac
#
pycup = [[np.transpose(np.linalg.inv(saaup)*pysmatup[0][0]),
          np.transpose(np.linalg.inv(sabup)*pysmatup[0][1])],
         [np.transpose(np.linalg.inv(sbaup)*pysmatup[1][0]),
          np.transpose(np.linalg.inv(sbbup)*pysmatup[1][1])]]

pycdown = [[np.transpose(np.linalg.inv(saadown)*pysmatdown[0][0]),
            np.transpose(np.linalg.inv(sabdown)*pysmatdown[0][1])],
           [np.transpose(np.linalg.inv(sbadown)*pysmatdown[1][0]),
            np.transpose(np.linalg.inv(sbbdown)*pysmatdown[1][1])]]
#
# compute Wmat
#
pywaaup = 0.0
pywabup = 0.0
pywbaup = 0.0
pywbbup = 0.0
for i in range(0, len(waaup)):
    for j in range(0, len(waaup)):
        pywaaup = pywaaup + waaup[i][j]*pycup[0][0][i][j]
        pywabup = pywabup + wabup[i][j]*pycup[0][1][i][j]
        pywbaup = pywbaup + wbaup[i][j]*pycup[1][0][i][j]
        pywbbup = pywbbup + wbbup[i][j]*pycup[1][1][i][j]
print """== W up python vs QE ==
[{}, {}]\t[{}, {}]
[{}, {}]\t[{}, {}]
""".format(round(pywaaup, 8),
           round(pywabup, 8),
           round(wmatup[0][0], 8),
           round(wmatup[0][1], 8),
           round(pywbaup, 8),
           round(pywbbup, 8),
           round(wmatup[1][0], 8),
           round(wmatup[1][1], 8))
#
#
pywaadown = 0.0
pywabdown = 0.0
pywbadown = 0.0
pywbbdown = 0.0
for i in range(0, len(waadown)):
    for j in range(0, len(waadown)):
        pywaadown = pywaadown + waadown[i][j]*pycdown[0][0][i][j]
        pywabdown = pywabdown + wabdown[i][j]*pycdown[0][1][i][j]
        pywbadown = pywbadown + wbadown[i][j]*pycdown[1][0][i][j]
        pywbbdown = pywbbdown + wbbdown[i][j]*pycdown[1][1][i][j]
print """== W down python vs QE ==
[{}, {}]\t[{}, {}]
[{}, {}]\t[{}, {}]
""".format(round(pywaadown, 8),
           round(pywabdown, 8),
           round(wmatdown[0][0], 8),
           round(wmatdown[0][1], 8),
           round(pywbadown, 8),
           round(pywbbdown, 8),
           round(wmatdown[1][0], 8),
           round(wmatdown[1][1], 8))
