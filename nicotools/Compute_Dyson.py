import sys
import numpy as np
from determinants import *

nstate_dyson = int(sys.argv[3])

# READS THE DETERMINANTS OF CISTATES 1 AND 2 AND LIST THOSE DIFFERING BY ONLY ONE ORBITAL (THE ONE "MISSING")
fcistate1 = open(sys.argv[1],"r")
mo_num, ndets1, nstate1 = ( int(x) for x in fcistate1.readline().split())
if(nstate_dyson>nstate1):
 print "nstate_dyson should be smaller than nstate1"
 sys.exit()
print 'Number of MO', mo_num
print 'Number of Determinants', ndets1

dets_alp1 = []
dets_bet1 = []
dets_typ1 = []
for i in range(ndets1):
  fcistate1.readline()
  deta = [' ']
  l = fcistate1.readline().split()
  deta += [ x for x in l[0] if (x=='+' or x=='-')]
  nmo_occ = deta.count('+')
  dets_alp1.append(deta)
  detb = [' ']
  l = fcistate1.readline().split()
  detb += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_bet1.append(detb)
  nholes=0
  for i in range(nmo_occ):
   if(not(deta[i+1])=='+'):
     nholes+=1
   if(not(detb[i+1])=='+'):
     nholes+=1
  dets_typ1.append(nholes)

fcistate2 = open(sys.argv[2],"r")
mo_num, ndets2, nstate2 = ( int(x) for x in fcistate2.readline().split())
print 'Number of MO', mo_num
print 'Number of Determinants', ndets2
print 

dets_alp2 = []
dets_bet2 = []
dets_typ2 = []
for i in range(ndets2):
  fcistate2.readline()
  deta = [' ']
  l = fcistate2.readline().split()
  deta += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_alp2.append(deta)
  detb = [' ']
  l = fcistate2.readline().split()
  detb += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_bet2.append(detb)
  nholes=0
  for i in range(nmo_occ):
   if(not(deta[i+1])=='+'):
     nholes+=1
   if(not(detb[i+1])=='+'):
     nholes+=1
#   if(nholes==2 and deta[0:nmo_occ+2]==detb[0:nmo_occ+2]):
#     nholes+=1
  dets_typ2.append(nholes)

ldets = []
for i in range(ndets1):
  for j in range(ndets2):
    #print dets_alp1[i],dets_alp2[j]
    ndiff, ndiff_alp, ndiff_bet, orbdiff_alp, orbdiff_bet = compare2dets(dets_alp1[i],dets_bet1[i],dets_alp2[j],dets_bet2[j])
    if(ndiff==1):
       print j,i, ndiff, orbdiff_alp,orbdiff_bet
       ldets.append((j,i,orbdiff_alp,orbdiff_bet))
#print ldets

# READS THE CI COEFFS
esta1 = [] 
two_e_sta1 = [] 
dip_sta = [] 
cista1 = [] 
nucl_rep = 0.0
w_singles = []
w_doubles = []
for i in range(nstate1):
  d = fcistate1.readline().split()
  nucl_rep = float(fcistate1.readline().split()[0])
  esta1.append(float(d[0])-nucl_rep)
  two_e_sta1.append(float(d[0])-float(d[1])-nucl_rep)
  dip_sta.append(float(d[2]))
  ws = 0.0
  wd = 0.0
  cicoeff = [] 
  for j in range(ndets1):
    d = fcistate1.readline().split()
    cicoeff.append(float(d[0]))
    if(dets_typ1[j]==1):
      ws += float(d[0])**2
    elif(dets_typ1[j]==2):
      wd += float(d[0])**2
  w_singles.append(ws)
  w_doubles.append(wd)
  cista1.append(cicoeff)

#print esta1
#print cista1

esta2 = []
cista2 = []
for i in range(nstate2):
  d = fcistate2.readline().split()
  nucl_rep = float(fcistate2.readline().split()[0])
  esta2.append(float(d[0])-nucl_rep)
  cicoeff = []
  for j in range(ndets2):
    cicoeff.append(float(fcistate2.readline().split()[0]))
  cista2.append(cicoeff)

#print esta2
#print cista2

# COMPUTES NORM OF THE DYSON ORBITALS
fdyson = open('Dyson_norms.txt','w')
for i1 in range(nstate_dyson):
  normtot = 0.0
  print i1, esta1[i1], two_e_sta1[i1], w_singles[i1], w_doubles[i1]
  print >> fdyson, i1, esta1[i1],two_e_sta1[i1], dip_sta[i1], w_singles[i1], w_doubles[i1]
  for i2 in range(nstate2):
    print "(N-1)e & (N)e states = ",i2,i1
    mocoeffs = np.zeros(mo_num+1)
    for l in ldets:
      j2, j1, orbalp, orbbet = l
      print j2, j1, orbalp, orbbet
      print cista1[i1][j1], cista2[i2][j2]
      if(len(orbalp)==1):
        orb=orbalp[0]
      else:
        orb=orbbet[0]
      scal = 1.0
#      if(dets_typ1[j1]==0 and dets_typ2[j2]==1): # GS/1h
#        scal = np.sqrt(2.0)
#      elif(dets_typ1[j1]==1 and dets_typ2[j2]==1): # 1h1p/1h
#        scal = 2.0
#      print j1,j2,dets_typ1[j1],dets_typ2[j2],scal
#      print orb, cista1[i1][j1], cista2[i2][j2],  np.sum(np.square(mocoeffs))
#      print
      mocoeffs[orb]+=cista1[i1][j1]*cista2[i2][j2]*scal
#      c1 = cista1[i1]
    normtot+=np.sum(np.square(mocoeffs))
    print >> fdyson, np.sum(np.square(mocoeffs)),
#    print 'Dyson orb. norm',np.sum(np.square(mocoeffs))
#    print 
  print >> fdyson, normtot
#  print normtot
#  print 
print "Dyson norms in Dyson_norms.txt"




