#!/usr/bin/python

# system libraries
import os
import subprocess
import sys

# calculate the time to run the code (I don't really ever use this)
def realTime(motifstderr):
  motifstderr = motifstderr.strip().split()
  if 'real' in motifstderr:
    time = motifstderr[motifstderr.index('real') + 1]
    time = time.replace('m',' ').replace('s',' ').split()
    time = 60*float(time[0]) + float(time[1])
    return time
  else:
    return 'NA'

# read in a network from stdin or from a data file
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

# make sure mfinder is compiled
def mfinderMake():
  try:
    process = subprocess.check_call(['make','-s','-C',os.path.dirname(__file__) + '/../mfinder1.2/'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,)
    return process
  except subprocess.CalledProcessError:
    print("error in compiling the mfinder1.2 package")
    sys.exit()

# calculate the motifs for a network by calling mfinder C++ executable
def mfinder(inFilename,outFilename,motifsize,nrandom,directed=True,maxmem=1000,time=False,membership=True):
  mfinderMake()

  command = list(map(str,['%s/../mfinder1.2/mfinder' % os.path.dirname(__file__),
                          inFilename,
                          '-s', motifsize,
                          '-r', nrandom,
                          '-f', outFilename,
                          '-q',
                          '-nu',]))

  if not directed:
    command.append('-nd')
    
  if membership:
    for i in map(str,['-omem','-maxmem',maxmem]):
      command.append(i)

  process = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,)

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
  motifs = map(int,stats.keys())
  motifs.sort()
  for m in map(str,motifs):
    print m+":", ' '.join([i for i in stats[m]])
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
def mfinderRoles(net,membership,motifsize=3,bipartite=False,weighted=False):
  try:
    weights = dict([((i,j),float(w)) for i,j,w in net])
  except ValueError:
    weights = dict([((i,j),1) for i,j in net])
  unet = weights.keys()

  if not bipartite:
    from unipartite_roles import unipartite_roles as motif_roles
  else:
    sys.stderr.write("Sorry, calculation of node roles in bipartite networks has not yet been implemented.\n")
    sys.exit(1)

  roles = {}
  for motif in membership:
    for motifnodes in membership[motif]:
      mindex = [i for i,j in motif_roles[motifsize]].index(motif)
      possible_roles = [tuple([motif]+list(r)) for r in motif_roles[motifsize][mindex][1]]

      if weighted:
        w = 0
        for node in motifnodes:
          for n in motifnodes:
            if node != n:
              try:
                w += weights[(node,n)]
              except KeyError:
                pass
        if w == int(w):
          w = int(w)
      else:
        w = 1  

      for node in motifnodes:
        indegree  = sum([(othernode,node) in unet for othernode in motifnodes if othernode != node])
        outdegree = sum([(node,othernode) in unet for othernode in motifnodes if othernode != node])
        key = (motif,indegree,outdegree)

        # if the node's in and out degrees are insufficient to discern its role
        # we will add the degrees of the nodes it interacts with (its neighbors)
        if key not in possible_roles:
          if indegree > 0:
            connected_to = [othernode for othernode in motifnodes if othernode != node and (othernode,node) in unet]
            outdegrees = [sum([(i,j) in unet for j in motifnodes if j != i]) for i in connected_to]
            outdegrees.sort()
            key = tuple(list(key) + [tuple(outdegrees)])
          else:
            connected_to = [othernode for othernode in motifnodes if othernode != node and (node,othernode) in unet]
            indegrees = [sum([(j,i) in unet for j in motifnodes if j != i]) for i in connected_to]
            indegrees.sort()
            indegrees = tuple(indegrees)
            key = tuple(list(key) + [tuple(indegrees)])

        # apparently the degrees of the node the focal node interacts with is still insufficient
        if key not in possible_roles:
          print >> sys.stderr, key
          print >> sys.stderr, "Apparently there is a role you aren't accounting for in roles.py."
          sys.exit()

        if not key in roles:
          roles[key] = {}

        try:
          roles[key][node] += w
        except:
          roles[key][node] = w

  nodes = list(set([i[0] for i in net] + [i[1] for i in net]))
  nodes.sort()

  noderoleprofiles = {}
  for node in nodes:
    noderoleprofiles[node] = {}
    for motif,rs in motif_roles[motifsize]:
      for r in rs:
        noderoleprofiles[node][tuple([motif] + list(r))] = 0

  for role in roles:
    for node in roles[role]:
      noderoleprofiles[node][role] = roles[role][node]

  return noderoleprofiles

# print out the role profiles species by species
def printRoleProfiles(noderoleprofiles,motifsize=3,bipartite=False):
  sortednodes = noderoleprofiles.keys()
  try:
    sortednodes = map(int,sortednodes)
  except:
    pass
  sortednodes.sort()
  sortednodes = map(str,sortednodes)
  
  if not bipartite:
    from unipartite_roles import unipartite_roles as motif_roles
  else:
    sys.stderr.write("Sorry, calculation of node roles in bipartite networks has not yet been implemented.\n")
    sys.exit(1)

  roles = []
  for motif,rs in motif_roles[motifsize]:
    roles += [tuple([motif] + list(r)) for r in rs]

  for node in sortednodes:
    print str(node)+":",
    print ' '.join([str(noderoleprofiles[node][role]) for role in roles])

  return

# run motif structure analysis (and, if desired, compare to randomized networks)
def runMotifAnalysis(net,motifsize=3,nrandom=250,directed=True):
  pid = str(os.getpid())
  netFilename = '/tmp/kk.mfinder.in.' + pid
  outFilename = '/tmp/kk.mfinder.out.' + pid

  motifFilename = outFilename + '_OUT.txt'

  writeTempNetwork(net,netFilename)
  pipestdout, pipestderr = mfinder(netFilename,outFilename,motifsize,nrandom,directed,membership=False,)

  ms = mfinderOutput(motifFilename)
  printMotifStatistics(ms)

  pipe = subprocess.Popen('rm %s %s' % (netFilename,motifFilename),
                          shell=True)
  pipe.communicate()

  return

# calculate species motif role profiles and print them out species by species
def runMotifProfileAnalysis(net,motifsize=3,bipartite=False,weighted=False):
  pid = str(os.getpid())
  netFilename = '/tmp/kk.mfinder.in.' + pid
  outFilename = '/tmp/kk.mfinder.out.' + pid

  motifFilename = outFilename + '_OUT.txt'
  memberFilename = outFilename + '_MEMBERS.txt'

  writeTempNetwork(net,netFilename)
  pipestdout, pipestderr = mfinder(netFilename,outFilename,motifsize,nrandom=0,directed=True)
  ms = mfinderOutput(motifFilename)

  maxnreal = max([int(ms[i][0]) for i in ms])
  if maxnreal > 1000:
    print >> sys.stderr, "upped maxmem to " + str(maxnreal+1)
    pipestdout, pipestderr = mfinder(netFilename,outFilename,motifsize,nrandom=0,directed=True,maxmem=maxnreal+1)
    ms = mfinderOutput(motifFilename)

  membership = mfinderMembership(memberFilename)
  profiles = mfinderRoles(net,membership,motifsize,bipartite,weighted)
  printRoleProfiles(profiles,motifsize,bipartite)

  pipe = subprocess.Popen('rm %s %s %s' % (netFilename,motifFilename,memberFilename),
                          shell=True)
  pipe.communicate()

  return

