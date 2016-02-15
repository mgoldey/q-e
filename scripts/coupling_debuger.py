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
pysmat = [[round(np.linalg.det(saaup), 8), round(np.linalg.det(sabup), 8)],
          [round(np.linalg.det(sbaup), 8), round(np.linalg.det(sbbup), 8)]]
print """== S up python vs QE ==
{}\t{}
{}\t{}
""".format(pysmat[0], detsup[0], pysmat[1], detsup[1])

# down
pysmat = [[round(np.linalg.det(saadown), 4), round(np.linalg.det(sabdown), 4)],
          [round(np.linalg.det(sbadown), 4), round(np.linalg.det(sbbdown), 4)]]
print """== S up python vs QE ==
{}\t{}
{}\t{}
""".format(pysmat[0], detsdown[0], pysmat[1], detsdown[1])

print sbaup[0, 1]
