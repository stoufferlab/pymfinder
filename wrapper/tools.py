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


motif_roles = {3:[(12, 0, 1),
                  (12, 1, 1),
                  (12, 1, 0),
                  (38, 0, 2),
                  (38, 1, 1),
                  (38, 2, 0),
                  (98, 1, 1),
                  (36, 0, 1),
                  (36, 2, 0),
                  (6, 0, 2),
                  (6, 1, 0),
                  (46, 1, 2),
                  (46, 2, 0),
                  (108, 0, 2),
                  (108, 2, 1),
                  (14, 1, 2),
                  (14, 1, 1),
                  (14, 1, 0),
                  (74, 0, 1),
                  (74, 1, 1),
                  (74, 2, 1),
                  (102, 1, 1,),
                  (102, 2, 1,),
                  (102, 1, 2,),
                  (238, 2, 2),
                  (110, 1, 2),
                  (110, 2, 1),
                  (110, 2, 2),
                  (78, 1, 1,),
                  (78, 2, 2,),
                  ],

               4:[(6, 0, 2),
                  (6, 1, 0),
                  (12, 0, 1),
                  (12, 1, 0),
                  (12, 1, 1),
                  (14, 1, 0),
                  (14, 1, 1),
                  (14, 1, 2),
                  (36, 0, 1),
                  (36, 2, 0),
                  (38, 0, 2),
                  (38, 1, 1),
                  (38, 2, 0),
                  (46, 1, 2),
                  (46, 2, 0),
                  (74, 0, 1),
                  (74, 1, 1),
                  (74, 2, 1),
                  (108, 0, 2),
                  (108, 2, 1),
                  (110, 1, 2),
                  (110, 2, 1),
                  (110, 2, 2),
                  (238, 2, 2),
                  ],
               }


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
  for pred, prey in net:
    if pred != prey:
      outFile.write(' '.join([str(i) for i in [pred, prey, 1]]))
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
  membership = {}

  memberstats = False
  inFile = open(filename,'r')
  for line in inFile:
    if 'subgraph id = ' in line:
      #try:
      #  print motif, nreal, len(membership[motif])
      #except:
      #  pass
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
        #print motif,sline
      #if '+-' in sline:
      #  print sline


    if 'Full list of ' in line:
      memberstats = True
      
  return membership

# determine motif roles given network and membership list
def mfinderRoles(net,membership,motifsize=3):
  roles = {}
  for motif in membership:
    for motifnodes in membership[motif]:
      for node in motifnodes:
        indegree  = sum([(othernode,node) in net for othernode in motifnodes if othernode != node])
        outdegree = sum([(node,othernode) in net for othernode in motifnodes if othernode != node])

        if not (motif,indegree,outdegree) in roles.keys():
          roles[(motif,indegree,outdegree)] = {}

        try:
          roles[(motif,indegree,outdegree)][node] += 1
        except:
          roles[(motif,indegree,outdegree)][node] = 1

  #print roles

  nodes = list(set([i[0] for i in net] + [i[1] for i in net]))
  nodes.sort()

  noderoleprofiles = {}
  for node in nodes:
    noderoleprofiles[node] = {}
    for role in motif_roles[motifsize]:
      noderoleprofiles[node][role] = 0

  for role in roles:
    for node in roles[role]:
      noderoleprofiles[node][role] = roles[role][node]

  return noderoleprofiles

def printRoleProfiles(noderoleprofiles,motifsize=3):
  sortednodes = noderoleprofiles.keys()
  sortednodes.sort()
  
  for node in sortednodes:
    print str(node)+":",
    print ' '.join([str(noderoleprofiles[node][role]) for role in motif_roles[motifsize]])


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

