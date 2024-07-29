# The different classes used to construct a population in the genetic
# algorithm.
from random import uniform, randint, shuffle, choice, sample
from math import sqrt, ceil
from types import IntType
from copy import deepcopy
from bisect import bisect
from structure import int2bp, bp2int, int2str, Structure
from individual import BaseDistribution, Individual
import mfe
from sys import stderr
import sys
import time
from shutil import rmtree
from subprocess import Popen,PIPE
from tempfile import mkdtemp,mkstemp
import multiprocessing
import os

wRNAfold=True
try:
    import RNAfold
except ImportError:
    wRNAfold=False

# Mutator class with multiprocessing option
class Mutator(multiprocessing.Process):
    
    def __init__(self,task,results):
        multiprocessing.Process.__init__(self)
        self.name = multiprocessing.current_process().name
        self.id=os.getpid()
        self.mytask=task
        self.daemon=True
        self.myresults=results
    
    def run(self):
        proc_name=self.name
        while True:
           
            task=self.mytask.get()
            result=task()
            self.myresults.put(result)
        print 'Finished: '+self.name
        return
    
global inputq
inputq=multiprocessing.Queue()
global resultq
resultq=multiprocessing.Queue()
global mutators
mutators=[]  
global server
  
### POPULATION; general class containing the core of a GA population
class Population:
  limit=2
  # Initialise population with n individuals. The target structure(s)
  # for the population is given by target, the class structure for
  # individuals by indclass, and basedist is used to draw nucleotides
  # when initialising and mutating sequences.

  def __init__(self, target, indclass, softtype = "RNAfold" , basedist = BaseDistribution(), n = 0, plimit=limit, motif = None):  


    self.members = set([]) # Current set of individuals
    self.new = set([]) # Set of individuals created waiting to be added
    self.target = target # Target structure(s) for population
    self.__individualclass__ = indclass # Class for creating individuals
    self.softtype=softtype #make sure all members of the population sure the same software type
    self.__basedistribution__ = basedist # Distribution to draw bases from
    self.plimit=plimit
    self.motif = motif
    
    global mutators
    global inputq
    global resultq
    
    # Start mutations
    for mut in mutators:
        if mut.is_alive()==False:
            mut.start()
    m=len(mutators)
    #Remove mutators if there are too many, add in more and start if there are not enough
    if m>plimit:
        for i in xrange(m-plimit):
            mutators.pop()

    # Add mutators if there aren't enough
    if m<plimit:
        for i in xrange(plimit-m):
            mut=Mutator(inputq,resultq)
            mutators.append(mut)
            
            
    # Add specified number of individuals
    for i in xrange(n):
        self.members.add(indclass(target=target,dist=self.__basedistribution__, softtype=self.softtype))

  # Return target structure(s) for this population
  def gettarget(self):
    return self.target
  # Return distribution for drawing new nucleotides
  def getbasedistribution(self):
    return self.__basedistribution__
    # Return software type
  def getsofttype(self):
    return self.softtype
  # Add n more individuals to population. If stoppingcriteria is not
  # None, it should be a function taking an iterator over the added
  # individuals as single argument, returning True if criteria for
  # stopping GA has been met and False otherwise, in which case
  # addrandom will return the same Boolean value.
  def addrandom(self, n, idfunc = None, stoppingcriteria = None):
    if idfunc == None:
      idfunc = (lambda x: None)
    
    def _addrandom():
        for i in xrange(n):
            new = self.__individualclass__(target=self.target, dist=self.__basedistribution__, id = idfunc(i), softtype=self.softtype, motif = self.motif)
            new.initialise()
            self.new.add(new)
            yield new
    if stoppingcriteria != None:
      return stoppingcriteria(_addrandom())
    else:
      for i in _addrandom():
        pass

  # Add waiting new individuals to current population
  def addnew(self):
    self.members.update(self.new)
    self.new = set([])
  # Size of population
  def __len__(self):
    return len(self.members)

  # Iterator for population
  def __iter__(self):
    for i in self.members:
      yield i

    
  #This function is to seperate the individuals in to groups and 
  #calculate their boltzman probabilities all at once with RNAfold 
  #rather than individually.      
  def massprecompute(self,individuals):
      if "__getprobabilityprediction__" in dir(self.__individualclass__):
          global mutators
          global resultq
          global inputq
          for mutator in mutators:
                  if mutator.is_alive()==False:
                      mutator.start()
          nmuts=len(mutators)
          m=list(individuals)
          indpermut=len(m)/nmuts
          remainder=len(m)%nmuts
        
          # Split the individuals evenly across the mutators.  The 
          for i in xrange(nmuts):
              
              # Find the indices in m of the individuals for each mutant
              # since slice takes up but not including the left index we want it one greater on the right side than necessary
              if i<remainder:
                  lower=i*(indpermut)+i*1
                  upper=lower+indpermut+1
              else:
                  if indpermut==0:
                      break
                  lower=i*indpermut+remainder
                  upper=(i+1)*indpermut+remainder
            
              #print lower,' ',upper,' ',m[lower:upper]
              inputq.put(MassPrecompute(set(m[lower:upper])))
          newinds=set([])
          totptimes=0
          while len(newinds)<len(m):
              x=resultq.get()
              newinds.update(x)
          
          return newinds
      else:
          return None
      
   # Set up function for making sure that we do not fold before
   # computing Boltzmann probabilities (as this would be a waste of
   # time). 
  def pprecompute(self,individuals):
      if "__getprobabilityprediction__" in dir(self.__individualclass__):
          
          global mutators
          global resultq
          global inputq
          for mutator in mutators:
              if mutator.is_alive()==False:
                  mutator.start()
                  #print 'starting processor'
          temperatures = None
          for i in individuals:
            if temperatures == None:
              temperatures = set([t.temperature for t in i.target])
              
            for t in temperatures:
              inputq.put(PreCompute(i,t))
          
          calculated=set([])
          avggetprobtime=0
          while len(calculated)<len(individuals):
            x=resultq.get()
            calculated.add(x)
            #avggetprobtime+=ptime
          return calculated
          
      else:
          return None

  # Calculate Boltzmann probabilities using specfied software
  def precompute(self,individuals):
      if "__getprobabilityprediction__" in dir(self.__individualclass__):
          temperatures = None
          for i in individuals: 
            if temperatures == None:
              temperatures = set([t.temperature for t in i.target])
            for t in temperatures:
                i.getprobability(t)
                
      else:
            return None

# Initiate class for Fitness, Mutation, Recombination, Crossover scoring,
# and Precompute
class Fitness():
    def __init__(self,monster,idnum):
        self.monster=monster
        self.num=idnum
    def __call__(self):
        return (self.monster.getfitness(),self.num)
    
class Mutation():
    def __init__(self,monster):
        self.monster=monster

    def __call__(self):
        return self.monster.mutate(self.monster.getposition())
        
class Recombination():
    
    def __init__(self,monster1,monster2,crosspoint):
        self.monster1=monster1
        self.monster2=monster2
        self.cp=crosspoint
    def __call__(self):
        return self.monster1.crossover(self.monster2,self.cp)
    

class CrossOverScore():
    
    def __init__(self,c1,id1,c2,id2,cuts):
        self.c1=c1
        self.c2=c2
        self.cuts=cuts
        self.id1=id1
        self.id2=id2
        
    
    def __call__(self):
        score=0
        
        
        for c in self.cuts:
            for i in c:
                for j in c:
                    if i != j:
                        score += pow(float(self.c1[(i, j)] + self.c2[(j, i)]), 2) / len(c)
            
        return self.id1,self.id2,score
class PreCompute():
    
    def __init__(self,monster,t):
        self.temp=t
        self.monster=monster
    def __call__(self):
        
        t1=time.time()
        self.monster.getprobability(self.temp)
        t2=time.time()
        ptime= (t2-t1)
        return self.monster
#Class used at collable in massprecompute to do all of the boltzman
#probabilities at once
class MassPrecompute():
    def __init__(self,monsters):
        self.monsters=monsters

    #We want to call RNAfold once and then pass in the sequences one at a time
    #Then collect the output and match it back to the monsters
    def __call__(self):
        
        global wRNAfold
        m=list(self.monsters)
        t1=time.time()
        sequences=[]
        m=list(self.monsters)
        #All of the monsters passed should be in the same population and thus
        #should have the same target(s),   Thus get the temps from the target of
        #a single member of the population 
        temperatures=set([tar.temperature for tar in m[0].target])
        
        for monster in m:
            sequences.append(monster.getseq())
        
        
        # Run sequences in parallel RNAfold  
        for T in temperatures:   
            if wRNAfold: 
                tempid,tempfile=mkstemp()
                args=" "+tempfile
                if T==None:
                    # Input into RNA fold must be a single string that is 
                    # delimited by spaces. (there should not be a space at the beginning of the string.
                    # The file name to write to should come first, followed by running flags.  
                    # Then the string 'sequences'  followed by all of the sequences.  Each sequences should be preceded by '>'
                    # and a unique index for referencing the output back to the group
                    args+=' -p -noPS sequences'
                else:
                    args+=' -p -noPS -T '+str(T)+' sequences'
                
                for i,seq in enumerate(sequences):
                    args+=' >'+str(i)
                    args+=' '+seq

                RNAfold.runRNAfold(args)
                outfile=open(tempfile,'r')
                output=outfile.read()
                os.remove(tempfile)
                seqoutput=output.split('>')
                for i in xrange(seqoutput.count('')):
                    seqoutput.remove('')
            else:
                if T==None:
                    cmd=['RNAfold','-p','-noPS']
                else:
                    cmd=['RNAfold','-p','-noPS','-T',str(T)]
                    
                p=Popen(cmd,stdin=PIPE,stdout=PIPE)
                for i,seq in enumerate(sequences):
                    print>>p.stdin,'>'+str(i)
                    print>>p.stdin,seq
                p.stdin.close()
                output=p.stdout.read()
                seqoutput=output.split('>')
                for i in xrange(seqoutput.count('')):
                    seqoutput.remove('')
            
            for out in seqoutput:
                lines=out.split('\n')
                ind=int(lines[0])
                try:
                    energies = lines[2].strip().split(None, 1)
                    ensemble = float(lines[3].strip().split(None, 1)[1][1:-1])
                    if energies == [] or len(energies[0]) != len(m[ind].seq):
                        # Did not receive expected output from RNAfold
                        raise IndexError
                except IndexError:
                    raise IndexError,'\n'+str(energies)+'\n'+str(ensemble)+'\n'+"Precompute did not receive appropriate output from RNAfold"
                
                p = [mfe._ProbVec() for j in xrange(len(m[ind].seq))]
                for j in p:
                    j[None] = 1.0

                # Read output
                datahere=lines[7:]
                for t in datahere:
                    if "ubox" in t:
                        t = t.split()
                        j = int(t[0]) - 1
                        k = int(t[1]) - 1
                        q = pow(float(t[2]), 2)
                        p[j][k] = q
                        p[k][j] = q
                        p[j][None] -= q
                        p[k][None] -= q             
                
                m[ind].prob[T]=(p,) 
                m[ind].struc[T]=Structure(energies[0],T)
                m[ind].optscore[T]=(float(energies[1][1:-1]),)
                m[ind].ensemble[T]=ensemble
        return set(m)
                    
            
            
        
### MUTATE; classes for selecting a set of mutation events creating
### new individuals.


# Simple class that uses each sequence in current population equally many times
class Mutate:
  # Add n mutated sequences to population (one for each individual in
  # current population if n is undefined), using each current
  # individual as starting point the same number of times.
  def mutate(self, n = None):
    if n == None:
      n = len(self)
    for i in self.members:
      for j in xrange(n / len(self)):
        self.new.add(i.mutate(i.getposition()))
        
    if n % len(self) != 0:
      # Some current individuals have to be used as staring point one
      # more time - choose these at random.
      c = sample(range(len(self)), n % len(self))
      c.sort()
      j = 0
      k = 0
      for i in self.members:
        if c[j] == k:
          self.new.add(i.mutate(i.getposition()))
          j += 1
          if j == len(c):
            break
        k += 1
    
  def pmutate(self, n = None, pp=False):
   
   
    global inputq
    global resultq
    global mutators
    m=list(self.members)
    if len(multiprocessing.active_children())==0:
        for mutator in mutators: 
            print 'starting processors'
            mutator.start()   
        
    if n == None:
      n = len(self)
      
    for j in xrange(n / len(self)):
        for i in xrange(len(self)):  
            inputq.put(Mutation(m[j]))
    if n % len(self) != 0:
      # Some current individuals have to be used as staring point one
      # more time - choose these at random.
      c = sample(range(len(self)), n % len(self))
      c.sort()
      j = 0
      k = 0
      
      for i in xrange(len(self.members)):
        if c[j] == k:
          inputq.put(Mutation(m[i]))
          j += 1
          if j == len(c):
            break
        k += 1
    mutants=set([])
    while len(mutants)<n:
        x=resultq.get()
        mutants.add(x)
    self.new.update(mutants)

# Class choosing starting points independently uniformly at random
class RandomMutate:
  ### Add n mutated sequences (the same number as in the current
  ### population if n is undefined)
  def mutate(self, n = None):
    if n == None:
      n = len(self)
    m = list(self.members)
    for j in xrange(n):
      i = choice(m)
      self.new.add(i.mutate(i.getposition()))
      
  def pmutate(self, n=None ):
    
    global inputq
    global resultq
    global mutators
    #We want to use parallelization to compute all of the fitnesses as well
    fitness={}
    if len(multiprocessing.active_children())==0:
        for mutator in mutators: 
            print 'starting processors'
            mutator.start()    
    if n == None:
        n = len(self)
    m = list(self.members)
   
    randominds=sample(range(len(self)),n) 
    
      
   
    for j in randominds: 
       x=Mutation(m[j])
       inputq.put(x)
    
  
    new=[]
    while len(new)<n:
        x=resultq.get()
        new.append(x)
        
            
    self.new.update(set(new))
    return
        
# Class choosing starting points according to fitness
class FitnessWeightedMutate:
  ### Add n mutated sequences (the same number as in the current
  ### population if n is undefined)
  def mutate(self, n = None):
   
    y=filter(lambda x: x > 0, map(lambda y: y.getfitness(), self.members))
    if len (y)>0:
        z=min(y)
    else:
         z=1
    def _reciprocal(x):
      if x <= 0:
        # A sequence with fitness of 0 is set to be twice as likely to
        # be picked as the sequence with the best fitness larger than
        # 0.
        return 2.0 / z
      else:
        return 1.0 / x
    if n == None:
      n = len(self)
      
    m = list(self.members)
    # Generate prefix sum array of fitnesses
    w = [_reciprocal(m[0].getfitness())] + (len(self) - 1) * [0]
    
    for i in xrange(1, len(self)):
      w[i] = w[i - 1] + _reciprocal(m[i].getfitness())
    
    
    new=set([])
    mutategenes=[]
    for j in xrange(n):  
        i=m[bisect(w, uniform(0, w[-1]), hi = len(self) - 1)]
        new.add(i.mutate(i.getposition()))
    self.new.update(new)    
    
  def pmutate(self, n = None):
    
    global inputq
    global resultq
    global mutators
    # We want to use parallelization to compute all of the fitnesses as well
    m=list(self.members)
    fit=[0]*len(m)
    if len(multiprocessing.active_children())==0:
        for i,mutator in enumerate(mutators): 
            print 'starting processors'
            mutator.start()
            
    
    for i,member in enumerate(m):
        inputq.put(Fitness(member,i))
    
    
    counter=0
    while counter<len(self.members):
                x=resultq.get()
                fit[x[1]]=x[0]
                counter=counter+1

         
    y=filter(lambda(x):x>0,fit)
    if len(y)>0:
        z=min(y)
    else:
        z=1
    def _reciprocal(x):
      if x <= 0:
        # A sequence with fit of 0 is set to be twice as likely to
        # be picked as the sequence with the best fit larger than
        # 0.
        return 2.0 / z
      else:
        return 1.0 / x
    
    if n == None:
      n = len(self)
    
   
    # Generate prefix sum array of fitnesses
    w = [_reciprocal(fit[0])] + (len(self) - 1) * [0]
    for i in xrange(1, len(self)):
      w[i] = w[i - 1] + _reciprocal(fit[i])

    for j in xrange(n): 
      x=Mutation(m[bisect(w, uniform(0, w[-1]), hi = len(self) - 1)])
      inputq.put(x)
    
    
    new=set([])
  
    counter=0
    
    while len(new)<n:
        x=resultq.get()
        new.add(x)
        
            
    self.new.update(new)
   
### RECOMBINE; classes for selecting a set of crossover events creating
### new individuals.

# Simple class just choosing pairs at random
class SelectPair:
  def __selectpairs__(self, n):
    m = list(self.members)
    return [sample(m, 2) for i in xrange(n)]
  #same as __selectpair__ for random
  def __pselectpairs__(self,n):
      m=list(self.member)
      return [sample(m, 2) for i in xrange(n)]
 

def score(c1,c2,cuts):
        score=0
        for c in cuts:#self.gettarget().getcuts():
          for i in c:
            for j in c:
              if i != j:  
               score += pow(float(c1[(i, j)] + c2[(j, i)]), 2) / len(c)
        return score     
# Class choosing n pairs based on fitness - it is assumed that fitness
# is a non-negative number, with 0 being a perfectly fitness individual.
class WeightedSelectPair:
  def __selectpairs__(self, n):
    
    z = min(filter(lambda x: x > 0, map(lambda y: y.getfitness(), self.members)))
    
    def _reciprocal(x):
      if x <= 0:
        # A sequence with fitness of 0 is set to be twice as likely to
        # be picked as the sequence with the best fitness larger than
        # 0.
        return 2.0 / z
      else:
        return 1.0 / x

    pairs = []
    m = list(self.members)
    # Generate prefix sum array of fitnesses
    w = [_reciprocal(m[0].getfitness())] + (len(self) - 1) * [0]
    for i in xrange(1, len(self)):
      w[i] = w[i - 1] + _reciprocal(m[i].getfitness())
    for i in xrange(n):
      # Choose first individual
      a = bisect(w, uniform(0, w[-1]), hi = len(self) - 1)
      # Choose second individual
      x = uniform(0, w[-1] - _reciprocal(m[a].getfitness()))
      if x >= w[a] - _reciprocal(m[a].getfitness()):
        # Choose beyond a
        b = a + 1 + bisect(w[a + 1:], x + _reciprocal(m[a].getfitness()), hi = len(self) - a - 2)
      else:
        b = bisect(w[:a], x, hi = a - 1)
      pairs.append((m[a], m[b]))
    return pairs
    
  def __pselectpairs__(self, n):
      
    
    
    
    
    global inputq
    global resultq
    global mutators
    # We want to use parallelization to compute all of the fitnesses as well
    fit={}
    if len(multiprocessing.active_children())==0:
        for mutator in mutators: 
            print 'starting processors'
            mutator.start()
    
        
    for i,member in enumerate(self.members):
        inputq.put(Fitness(member,i))
    
    
    
    while len(fit.keys())<len(self.members):
                x=resultq.get()
                fit[x[1]]=x[0]
            
    
    y=filter(lambda(x):x>0,fit.values())
    if len(y)>0:
        z=min(y)
    else:
        z=1
            
    
    def _reciprocal(x):
      if x <= 0:
        # A sequence with fit of 0 is set to be twice as likely to
        # be picked as the sequence with the best fit larger than
        # 0.
        return 2.0 / z
      else:
        return 1.0 / x

    pairs = []
    m = list(self.members)
    # Generate prefix sum array of fitnesses
    w = [_reciprocal(fit[0])] + (len(self) - 1) * [0]
    for i in xrange(1, len(self)):
      w[i] = w[i - 1] + _reciprocal(fit[i])
    for i in xrange(n):
      # Choose first individual
      a = bisect(w, uniform(0, w[-1]), hi = len(self) - 1)
      # Choose second individual
      x = uniform(0, w[-1] - _reciprocal(fit[a]))
      if x >= w[a] - _reciprocal(fit[a]):
        # Choose beyond a
        b = a + 1 + bisect(w[a + 1:], x + _reciprocal(fit[a]), hi = len(self) - a - 2)
      else:
        b = bisect(w[:a], x, hi = a - 1)
      pairs.append((m[a], m[b]))
    return pairs


# Class choosing pairs based on sum of squares of scores of all
# possible cross-overs between the pair.def cross_score(c1,c2):
class CombinedSelectPair:
  def __selectpairs__(self, n):
    
    # Function computing one dimensional index from pair of indeces
    def pair2idx(a, b):
      return a * (a - 1) / 2 + b
    # Function computing pair corresponding to one dimensional index
    def idx2pair(index):
      a = int(sqrt(2 * index)) # Either a or a + 1
      a = int(sqrt(2 * index + a)) # 2x is between a^2 - a and a^2 + a - 2
      b = index - a * (a - 1) / 2  # Once we know a, b is easy
      return a, b
    # Compute weight for combining all pairs
    m = list(self.members)
    w = len(m) * (len(m) - 1) / 2 * [0.0]

    indxa = sample(xrange(1, len(m)), int(ceil(len(m)/2)) )
    for a in indxa:
      indxb = sample(xrange(a), int(ceil(a/2)))
      for b in indxb:
        # For each pair, the weight is the combined fitness summed
        # over all crossover points.
        idx = pair2idx(a, b)

        for c in self.gettarget().getcuts():
          for i in c:
            for j in c:
              if i != j:  
                w[idx] += pow(float(m[a].cweight((i, j)) + m[b].cweight((j, i))), 2) / len(c)
    
    
    # Choose n pairs according to the computed weights
    pairs = []
    # Change w to prefix sum array
    for i in xrange(1, len(w)):
      w[i] += w[i - 1]
    for i in xrange(n):
      a, b = idx2pair(bisect(w, uniform(0, w[-1]), hi = len(w) - 1))
      pairs.append((m[b], m[a]))
    
    return pairs
  
  def __pselectpairs__(self, n):
    # Function computing one dimensional index from pair of indeces
    t=time.time()
    def pair2idx(a, b):
      return a * (a - 1) / 2 + b
    # Function computing pair corresponding to one dimensional index
    def idx2pair(index):
      a = int(sqrt(2 * index)) # Either a or a + 1
      a = int(sqrt(2 * index + a)) # 2x is between a^2 - a and a^2 + a - 2
      b = index - a * (a - 1) / 2  # Once we know a, b is easy
      return a, b
    
      
    
    
    
    global inputq
    global resultq
    global mutators
    
    if len(multiprocessing.active_children())==0:
        for mutator in mutators: 
            print 'starting children processes'
            mutator.start()
            
    m = list(self.members)
    w = len(m) * (len(m) - 1) / 2 * [0.0]
    cuts=self.gettarget().getcuts()
    for i in xrange(len(w)):
        a,b=idx2pair(i)
        inputq.put(CrossOverScore(m[a].cweight(all=True),a,m[b].cweight(all=True),b,cuts))
        
    counter=0
    while counter<len(w):
        x=resultq.get()
        counter+=1
        idx=pair2idx(x[0],x[1])
        # Go from index to the number of members minus 1 for the one being crossed
        w[idx]=x[2]
    
    pairs = []

    # Change w to prefix sum array
    for i in xrange(1, len(w)):
      w[i] += w[i - 1]
    for i in xrange(n):
      a, b = idx2pair(bisect(w, uniform(0, w[-1]), hi = len(w) - 1))
      pairs.append((m[b], m[a]))

    return pairs
    
# Class for choosing cross over points uniformly at random
class SelectCut:
  # Select a cut for recombining a and b such that positions that can
  # be part of a cut are chosen uniformly at random.
  def __selectcut__(self, a, b):
    cuts = filter(lambda x: len(x) > 1, self.gettarget().getcuts())
    w = map(lambda x: len(x), cuts)
    for i in xrange(1, len(w)):
      w[i] += w[i - 1]
    # Select dependency class at random weighted by size
    i = bisect(w, randint(0, w[-1] - 1), hi = len(w) - 1)
    # Choose random pair of positions in dependency class
    return sample(cuts[i], 2)

# Class for choosing cross over points weighted by current fitness
class WeightedSelectCut:
  def __selectcut__(self, a, b):
    # Compute prefix sums of fitnesses
    w = []
    cuts = []
   
    for c in self.gettarget().getcuts():
      for i in c:
        for j in c:
          if i<j:
            w.append(a.cweight((j, i)) + b.cweight((i, j)))
            cuts.append((i, j))
    
    for i in xrange(1, len(w)):
      w[i] += w[i - 1]
    # Select a cut proportional to its weight
    i = bisect(w, uniform(0, w[-1]), hi = len(w) - 1)
    return cuts[i]

### RECOMBINE; classes for creating recombinants

class Recombine:
  # Add n sequences obtained by crossover (the same number as in the
  # current population if n is undefined).
  
  def recombine(self, n = None):
    # Choose pairs to create recombinants from
    if n == None:
      n = len(self)
    pairs = self.__selectpairs__(n)
    
    # Choose recombination point for each pair, and create recombinant
    for p in pairs:
      c = self.__selectcut__(p[0], p[1])
      self.new.add(p[0].crossover(p[1], c))
      
   
      
         
  def precombine(self, n=None ):
      # Choose pairs to create recombinants from
    if n == None:
      n = len(self)  
    
    pairs = self.__pselectpairs__(n)
    
    #Have to make a list of cut point before calling parallelized version
    # Commenting here - William
    paircuts=[]
    for i in range(len(pairs)):
        paircuts.append(self.__selectcut__(pairs[i][0],pairs[i][1]))
        
    
    

    global inputq
    global resultq
    global mutators
    
    for mutator in mutators:
              if mutator.is_alive()==False:
                  mutator.start()
    
    
      
    for i,pair in enumerate(pairs):
         inputq.put(Recombination(pair[0],pair[1],paircuts[i]))  
    
    
    
    new=[]
    
    
    while len(new)<(n):
        x=resultq.get()
        new.append(x)
         
    self.new.update(new)
    
    
### REDUCE; Classes for reducing the population size

# Method for eliminating duplicates in m, as long as we still have
# at least n remaining elements.
def _removeduplicates(m, n):
  

  s = set([])
  k = len(m) - 1
  i = 0
  while i <= k:
    if tuple(m[i].seq) in s:
      # Seen this sequence before, eliminate from m
      m[i], m[k] = m[k], m[i]
      k -= 1
      if k < n:
        # Can't afford to eliminate any more
        break
    else:
      # New sequece, remember and proceed with next one
      s.add(tuple(m[i].seq))
      i += 1
  del m[k + 1:]
  return m

# Simple class just retaining the fittest individuals
class Reduce:
  def reduce(self, n):
    if len(self) > n:
      m = _removeduplicates(list(self.members), n)
      shuffle(m)
      m.sort(lambda x, y: int(x.getfitness() >= y.getfitness()) - int(x.getfitness() < y.getfitness()))
      self.members = set(m[:n])
      
# Class reducing the population using a Pareto method
class ParetoReduce:
  def reduce(self, n):   
    if len(self) > n:
        m = _removeduplicates(list(self.members), n)
        shuffle(m)
        F = []
        # Hold members in S based on their rank
        S = dict()
        junk =[]
        for p in m:
            S[hash(p)] = []
            p.n = 0
            # Count number of dominating members
            for q in m:
                if p.dominated(q):
                    S[hash(p)].append(q)
                elif q.dominated(p):
                    p.n += 1
            if p.n == 0:
                p.rank = 1
                F.append([])
                F[0].append(p)
                junk.append(p)
        i = 0
        # Append best members to junk and delete until n remain
        while len(F[i]) != 0:
            Q = []
            for p in F[i]:
                for q in S[hash(p)]:
                    q.n -= 1
                    if q.n == 0:
                        q.rank = i + 1
                        Q.append(q)
                        junk.append(q)
            i += 1
            F.append([])
            F[i] = Q # Contains ranks as the index and holds all members of that rank inside       
    self.members = set(junk[:n])
      
  
# Class selecting based on a weighted combination of fitness and
# diversity. This class is returned by a function that takes the
# weight as parameter, where the score of an individual in a
# particular stage of the culling is computed as the fitness plus the
# weight times the average fraction of positions where a sequence
# is identical to the already selected sequences.
def DiversityReduce(weight):
  class _DiversityReduce:
    def reduce(self, n):
      if len(self) > n:
        def sim(a, b):
          c = 0
          for i in xrange(len(a.seq)):
            if a.seq[i] == b.seq[i]:
              c += 1
          return float(c) / len(a.seq)
        m = _removeduplicates(list(self.members), n)
        shuffle(m)
        w = [i.getfitness() for i in m]
        # Find most fitness individual
        keep = 0
        for i in xrange(1, len(m)):
          if w[i] < w[keep]:
            keep = i
        m[keep], w[keep], keep = m[-1], w[-1], [m[keep]]

        # Select remaining n - 1 individuals
        d = [sim(i, keep[0]) for i in m[:-1]]
        for i in xrange(1, n):
          # Find best score
          best = w[0] + d[0] * weight / float(i)
          besti = 0
          for j in xrange(1, len(m) - i):
            s = w[j] + d[j] * weight / float(i)
            if s < best:
              best = s
              besti = j
          # Move best scoring individual to keep and remove it from remaining
          keep.append(m[besti])
          if i != n - 1:
            m[besti] = m[-1 - i]
            w[besti] = w[-1 - i]
            d[besti] = d[-1 - i]
            # Update diversity scores
            for j in xrange(len(m) - 1 - i):
              d[j] += sim(m[j], keep[-1])
        # All done selecting the keepers
        self.members = set(keep)
  return _DiversityReduce
