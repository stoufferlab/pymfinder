#!/usr/bin/python

# system libraries
import copy
from math import sqrt
import os
import random
import subprocess
import sys

# global variables
motif_translation = {"12":"S1",
                     "38":"S2",
                     "98":"S3",
                     "36":"S4",
                     "6":"S5",
                     "46":"D1",
                     "108":"D2",
                     "14":"D3",
                     "74":"D4",
                     "102":"D5",
                     "238":"D6",
                     "110":"D7",
                     "78":"D8",}

from roles import unipartite_roles as motif_roles

"""
def normalizeZScores(mstats,normalize=True):
  if normalize:
    zsum = 0.0
    for m in mstats:
      zsum += mstats[m][-1]**2

    zsum = sqrt(zsum)

    for m in mstats:
      real, zscore = mstats[m]
      mstats[m] = (real, zscore, zscore/float(zsum))

    return mstats
  else:
    for m in mstats:
      real, zscore = mstats[m]
      mstats[m] = (real, zscore, 0)
    return mstats
"""

# print out the motif structure
# can print to file or to stdout and append or not
def printMotifAnalysis(originalFilename,ms,stdout=False,append=True):
  if not stdout:
    outputFilename = '/'.join(originalFilename.split("/")[:-1]) + '/motif_stats.dat'
    if append:
      outFile = open(outputFilename,'a')
    else:
      outFile = open(outputFilename,'w')
  else:
    outFile = sys.stdout

  for i in [0,1,2]:
    for m in ['S1','S2','S3','S4','S5','D1','D2','D3','D4','D5','D6','D7','D8',]:
      outFile.write(str(ms[m][i])+' ')

  outFile.write('\n')
  outFile.flush()

  if not stdout:
    outFile.close()

def realTime(motifstderr):
  motifstderr = motifstderr.strip().split()
  if 'real' in motifstderr:
    time = motifstderr[motifstderr.index('real') + 1]
    time = time.replace('m',' ').replace('s',' ').split()
    time = 60*float(time[0]) + float(time[1])
    return time
  else:
    return 'NA'

# read in a network
def readNetwork(innie):
  net = []
  if type(innie) == type(sys.stdin):
    closeme = False
  elif type(innie) == type(''):
    innie = open(innie,'r')
    closeme = True
  else:
    print >> sys.stderr, "What do you want me to do now?"
    
  for line in innie:
    sline = line.strip().split()
    net.append(tuple(sline))

  if closeme:
    innie.close()

  return net  

# write out a network so that mfinder works
def writeTempNetwork(net,outFilename):
  outFile = open(outFilename,'w')
  try:
    for pred, prey in net:
      if pred != prey:
        outFile.write(' '.join(map(str,[pred, prey, 1])))
        outFile.write('\n')
  except ValueError:
    for pred, prey, weight in net:
      if pred != prey:
        outFile.write(' '.join(map(str,[pred, prey, weight])))
        outFile.write('\n')
  outFile.close()

# calculate the motifs for a network calling C++ executable
def mfinder(inFilename,outFilename,motifsize,nrandom,directed=True,maxmem=1000,time=False,membership=True):
  command = "%s %s/../mfinder1.2/mfinder %s -s %s -r %s %s -f %s -q -nu %s"
  
  if membership:
    membershipstring = '-omem -maxmem %s' % maxmem
  else:
    membershipstring = ''

  if time:
    timestring = 'time'
  else:
    timestring = ''

  undirectedflag = ''
  if not directed:
    undirectedflag = '-nd'

  parameters = (timestring, os.path.dirname(__file__), inFilename, motifsize, nrandom, undirectedflag, outFilename, membershipstring)

  process = subprocess.Popen(command % parameters,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)

  return process.communicate()


# parse the output file with counts, z-scores, etc
def mfinderOutput(filename):
  inFile = open(filename,'r')

  stats = {}

  motifstats = False
  for line in inFile:
    if motifstats:
      if '+-' in line:
        motif, nreal, nrand, z, p, conc, uniq = line.strip().split()
        nrand, stdevrand = nrand.split('+-')
        stats[motif] = (nreal,nrand,stdevrand,z,p,conc,uniq)

    if 'Full list of subgraphs size ' in line:
      motifstats = True

  return stats

def printMotifStatistics(stats):
  motifs = [int(i) for i in stats.keys()]
  motifs.sort()

  for m in motifs:
    print str(m)+":", ' '.join([i for i in stats[str(m)]])
  
  return
  
# parse in which motifs, and in what combination, individual species appear
def mfinderMembership(filename):
  inFile = open(filename,'r')

  membership = {}
  memberstats = False
  for line in inFile:
    if 'subgraph id = ' in line:
      motif = int(line.strip().split()[-1])
      memberstats = False
    elif 'Nreal : ' in line:
      nreal = int(line.strip().split()[-1])

    if memberstats:
      sline = tuple([int(i) for i in line.strip().split()])
      sline = tuple(line.strip().split())
      if sline:
        try:
          membership[motif].append(sline)
        except:
          membership[motif] = [sline]

    if 'Full list of ' in line:
      memberstats = True

  inFile.close()
  return membership

# determine motif roles given network and membership list
def mfinderRoles(net,membership,motifsize=3):
  try:
    weights = dict([((i,j),float(w)) for i,j,w in net])
  except ValueError:
    weights = dict([((i,j),1) for i,j in net])
  unet = weights.keys()

  roles = {}
  for motif in membership:
    for motifnodes in membership[motif]:
      mindex = [i for i,j in motif_roles[motifsize]].index(motif)
      possible_roles = [tuple([motif]+list(r)) for r in motif_roles[motifsize][mindex][1]]
      for node in motifnodes:
        indegree  = sum([(othernode,node) in unet for othernode in motifnodes if othernode != node])
        outdegree = sum([(node,othernode) in unet for othernode in motifnodes if othernode != node])
        key = (motif,indegree,outdegree)

        if key not in possible_roles:
          print >> sys.stderr, key
          print >> sys.stderr, "Apparently there is a role you aren't accounting for in roles.py."
          sys.exit()

        if not key in roles:
          roles[key] = {}

        try:
          roles[key][node] += 1
        except:
          roles[key][node] = 1

  nodes = list(set([i[0] for i in net] + [i[1] for i in net]))
  nodes.sort()

  noderoleprofiles = {}
  for node in nodes:
    noderoleprofiles[node] = {}
    for motif, rs in motif_roles[motifsize]:
      for r in rs:
        noderoleprofiles[node][tuple([motif] + list(r))] = 0

  for role in roles:
    for node in roles[role]:
      noderoleprofiles[node][role] = roles[role][node]

  return noderoleprofiles

def printRoleProfiles(noderoleprofiles,motifsize=3):
  sortednodes = noderoleprofiles.keys()
  try:
    sortednodes = map(int,sortednodes)
  except:
    pass
  sortednodes.sort()
  sortednodes = map(str,sortednodes)
  
  roles = []
  for motif, rs in motif_roles[motifsize]:
    roles += [tuple([motif] + list(r)) for r in rs]

  for node in sortednodes:
    print str(node)+":",
    print ' '.join([str(noderoleprofiles[node][role]) for role in roles])


# 
def runMotifAnalysis(net,motifsize=3,nrandom=250,directed=True):
  pid = str(os.getpid())
  netFilename = '/tmp/kk.mfinder.in.' + pid
  outFilename = '/tmp/kk.mfinder.out.' + pid

  motifFilename = outFilename + '_OUT.txt'
  #memberFilename = outFilename + '_MEMBERS.txt'

  writeTempNetwork(net,netFilename)
  pipestdout, pipestderr = mfinder(netFilename,outFilename,motifsize,nrandom,directed,membership=False,)

  ms = mfinderOutput(motifFilename)
  printMotifStatistics(ms)

  #pipe = subprocess.Popen('rm %s %s %s' % (netFilename,motifFilename,memberFilename),
  pipe = subprocess.Popen('rm %s %s' % (netFilename,motifFilename),
                          shell=True)
  pipe.communicate()

  return


# calculate species motif role profiles and print them out species by species
def runMotifProfileAnalysis(net,motifsize=3,nrandom=250,directed=True):
  pid = str(os.getpid())
  netFilename = '/tmp/kk.mfinder.in.' + pid
  outFilename = '/tmp/kk.mfinder.out.' + pid

  motifFilename = outFilename + '_OUT.txt'
  memberFilename = outFilename + '_MEMBERS.txt'

  writeTempNetwork(net,netFilename)
  pipestdout, pipestderr = mfinder(netFilename,outFilename,motifsize,nrandom,directed)
  ms = mfinderOutput(motifFilename)

  maxnreal = max([int(ms[i][0]) for i in ms])

  if maxnreal > 1000:
    print >> sys.stderr, "upped maxmem to " + str(maxnreal+1)
    pipestdout, pipestderr = mfinder(netFilename,outFilename,motifsize,nrandom,directed,maxmem=maxnreal+1)
    ms = mfinderOutput(motifFilename)

  membership = mfinderMembership(memberFilename)

  pipe = subprocess.Popen('rm %s %s %s' % (netFilename,motifFilename,memberFilename),
                          shell=True)
  pipe.communicate()

  profiles = mfinderRoles(net,membership,motifsize)

  printRoleProfiles(profiles)

  return

