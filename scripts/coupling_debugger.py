import numpy as np
import scipy as sp
from scipy import linalg
from os import sys

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

# read h, f and corrections
fandc = np.loadtxt('FandC')
hmat = np.loadtxt('H')

if (sys.argv[1] == 'w'):
    # read W, vecs of W and horth
    w = np.loadtxt('W')
    v = np.loadtxt('Weigvec')
    horth = np.loadtxt('Horth')

if (sys.argv[1] == 'l'):
    # read W, vecs of W and horth
    sinvoh = np.loadtxt('sinvoh')
    horth = np.loadtxt('HorthLow')

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
pycup = [[np.transpose(np.linalg.inv(saaup)*detsup[0][0]),
          np.transpose(np.linalg.inv(sabup)*detsup[0][1])],
         [np.transpose(np.linalg.inv(sbaup)*detsup[1][0]),
          np.transpose(np.linalg.inv(sbbup)*detsup[1][1])]]

pycdown = [[np.transpose(np.linalg.inv(saadown)*detsdown[0][0]),
            np.transpose(np.linalg.inv(sabdown)*detsdown[0][1])],
           [np.transpose(np.linalg.inv(sbadown)*detsdown[1][0]),
            np.transpose(np.linalg.inv(sbbdown)*detsdown[1][1])]]
#
# compute Wmat use w from QE and
# the cofactor computed from python
#

# up
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
# down
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


#
# compute h
#
# <a|H|b>
# <b|H|a>
# hab = ((fandc[0][0]*pysmatup[1][0]*pysmatdown[1][0] +
#        fandc[1][0]*pysmatup[0][1]*pysmatdown[0][1]) -
#       (pywabup*pysmatdown[0][1] + pywabdown*pysmatup[0][1] +
#        pywbaup*pysmatdown[1][0] + pywbadown*pysmatup[1][0]))
# hab = 0.5*hab
wtotab = pywabup*detsdown[0][1] + pywabdown*detsup[0][1]
wtotba = pywbaup*detsdown[1][0] + pywbadown*detsup[1][0]

fba = fandc[0][0]*detsup[1][0]*detsdown[1][0]
fab = fandc[1][0]*detsup[0][1]*detsdown[0][1]

pyhmat = [[fandc[0][0]+fandc[0][1],
           fab - wtotab],
          [fba - wtotba,
           fandc[1][0]+fandc[1][1]]]

pyhmat[1][0] = 0.5*(pyhmat[1][0] + pyhmat[0][1])
pyhmat[0][1] = pyhmat[1][0]

print """== H python vs QE ==
[{}, {}]\t[{}, {}]
[{}, {}]\t[{}, {}]
""".format(round(pyhmat[0][0], 8),
           round(pyhmat[0][1], 8),
           round(hmat[0][0], 8),
           round(hmat[0][1], 8),
           round(pyhmat[1][0], 8),
           round(pyhmat[1][1], 8),
           round(hmat[1][0], 8),
           round(hmat[1][1], 8))
#
# compute H orth using W
#
if (sys.argv[1] == 'w'):
    stot = np.array([[detsup[0][0]*detsdown[0][0],
                      detsup[0][1]*detsdown[0][1]],
                     [detsup[1][0]*detsdown[1][0],
                      detsup[1][1]*detsdown[1][1]]])
    stotinv = np.linalg.inv(stot)
    # W.V = S.V.L
    # S^-1.W.V = V.L
    sinvw = np.dot(stotinv, w)
    # w_eigval, w_eigenlvec, w_eigenvec = sp.linalg.eig(w, b=None,
    # left=False, right=True, overwrite_a=False, overwrite_b=False,
    # check_finite=True)
    # print sp.linalg.eig(w, b=stot, left=True, right=True,
    # overwrite_a=False, overwrite_b=False, check_finite=True)
    w_eigval, w_eigvec = np.linalg.eig(sinvw)
    #
    # V^dag.H.V
    #
    # w_eigvec = v
    w_eigvec_dag = np.transpose(w_eigvec)
    tmp = np.dot(w_eigvec_dag, hmat)
    pyhorth = np.dot(tmp, w_eigvec)
    print """
== H orth by W python vs QE ==
[{}, {}]\t[{}, {}]
[{}, {}]\t[{}, {}]
""".format(round(pyhorth[0][0], 8),
           round(pyhorth[0][1], 8),
           round(horth[0][0], 8),
           round(horth[0][1], 8),
           round(pyhorth[1][0], 8),
           round(pyhorth[1][1], 8),
           round(horth[1][0], 8),
           round(horth[1][1], 8))


if (sys.argv[1] == 'l'):
    # compute S^-1/2
    stot = np.array([[detsup[0][0]*detsdown[0][0],
                      detsup[0][1]*detsdown[0][1]],
                     [detsup[1][0]*detsdown[1][0],
                      detsup[1][1]*detsdown[1][1]]])
    seig, sv = np.linalg.eig(stot)
    sdiag = np.dot(np.linalg.inv(sv), stot)
    sdiag = np.dot(sdiag, sv)
    sdiagsinvoh = [[sdiag[0][0]**(-0.5),
                    sdiag[0][1]],
                   [sdiag[1][0],
                    sdiag[1][1]**(-0.5)]]
    sdiagsinvoh = np.dot(sdiagsinvoh, np.linalg.inv(sv))
    pysinvoh = np.dot(sv, sdiagsinvoh)
    #
    # compute S^-1/2 H S^-1/2
    #
    pyhorth = np.dot(pysinvoh, hmat)
    pyhorth = np.dot(pyhorth, pysinvoh)
    print """
== S^-1/2 python vs QE ==
[{}, {}]\t[{}, {}]
[{}, {}]\t[{}, {}]
""".format(round(pysinvoh[0][0], 8),
           round(pysinvoh[0][1], 8),
           round(sinvoh[0][0], 8),
           round(sinvoh[0][1], 8),
           round(pysinvoh[1][0], 8),
           round(pysinvoh[1][1], 8),
           round(sinvoh[1][0], 8),
           round(sinvoh[1][1], 8))

    print """
== H orth by L python vs QE ==
[{}, {}]\t[{}, {}]
[{}, {}]\t[{}, {}]
""".format(round(pyhorth[0][0], 8),
           round(pyhorth[0][1], 8),
           round(horth[0][0], 8),
           round(horth[0][1], 8),
           round(pyhorth[1][0], 8),
           round(pyhorth[1][1], 8),
           round(horth[1][0], 8),
           round(horth[1][1], 8))
