#Adaptation of mfe to use pknotsRG
from subprocess import Popen, PIPE, check_call
from tempfile import mkdtemp, mkstemp
import re
from shutil import rmtree
from os.path import join
from structure import int2str, Structure, pairing, pkpairing
from time import sleep
from types import StringType
import os ,time,sys
from random import choice

# Check shortened path call works, otherwise extend
wRNAfold=True
from time import time
try:
    import RNAfold
except ImportError:
    wRNAfold=False
ViennaRNA_location = None  # Only needs to be specified if the ViennaRNA collection of executables are not on the system path
_ViennaRNA_add_location = lambda executable: join(ViennaRNA_location, executable) if ViennaRNA_location else executable
os.environ["PATH"]=os.environ["PATH"]+":/usr/local/bin"

# Function defining behaviour if a folding program fails
def onfail(f, tries, *args):
    print>>sys.stderr, "Unexpected output from folding program - attempt", tries+1 , "of 10 to recover"
    if tries < 10:
    # Don't retry more than 10 times
        sleep(10) # Take ten
        return f(*args, tries = tries + 1)
    else:
        raise RuntimeError, "Function " + f.__name__ + " did not receive result from external program"

# Compute MFE structure of RNA sequence s, return pair of structure in
# bracket notation and MFE as float (pk... uses pknotsRG, ... uses default)
def fold(s, T = None, tries = 0):
  cmd = [_ViennaRNA_add_location("RNAfold")]
  tmpdir = mkdtemp()
  if T != None:
    cmd.extend(["-T", str(T)])
  p = Popen(cmd, stdin = PIPE, stdout = PIPE, cwd = tmpdir)
  print >> p.stdin, int2str(s)
  p.stdin.close()
  t = p.stdout.readlines()[-1].strip().split(None, 1)
  p.stdout.close()
  rmtree(tmpdir, ignore_errors = True)
  if t == [] or len(t[0]) != len(s):
    # Did not receive expected output from RNAfold
    return onfail(fold, tries, s, T)
  #print "from normal fold: struct", t[0]
  #print "energy:", t
  return t[0], float(t[1][1:-1])

# Compute MFE structure of RNA sequence s, return pair of structure in
# bracket notation and MFE as float

def pfold(s, T = None, tries = 0):
  #make file for ppfold to read 
  cmd=['java','-jar','PPfold3.1_special.jar', '--seq', int2str(s)]

  p=Popen(cmd, stdin = PIPE, stdout = PIPE, cwd=os.getcwd())
  t = p.stdout.readlines()
  while p.stdout.readline() != "":
    # Make sure RNAfold has properly finished before reading probability plot
    pass
  p.stdout.close()

  struct = t[0].strip("][''\n")
  
  cmd = [_ViennaRNA_add_location("RNAeval")]
  tmpdir = mkdtemp()
  pop = Popen(cmd, stdin = PIPE, stdout = PIPE, cwd = tmpdir)
  print >> pop.stdin, int2str(s)
  print >> pop.stdin, struct
  pop.stdin.close()
  u = pop.stdout.readlines()[-1].strip().split(None, 1)
  pop.stdout.close()
  rmtree(tmpdir, ignore_errors = True)
  if u == [] or len(u[0]) != len(s):
    # Did not receive expected output from RNAfold
    return onfail(score, tries, s, t)
  
  return struct, float(u[1][1:-1])



def pkfold(s, T = None, tries = 0):
    

  #os.environ['PATH']=os.environ['PATH']+":/Users/whweir/Documents/Python/CompBioProject/src/frnakenstein"
  cmd = ["pknotsRG"]
  tstrt=time()
  try:
      p = Popen(cmd, stdin = PIPE, stdout = PIPE, bufsize=0)
  except OSError:
      print OSError.args,'Trying again'
      return onfail(pkfold,tries,s,T)
      
  
  print >> p.stdin, int2str(s)
  p.stdin.flush()
  p.stdin.close()
  p.wait()
  out=p.stdout.read()
  p.stdout.close()
  lines=out.split('\n')
  q=lines[1]
  temp=q.split()
  u = re.search("-?\d+\.+\d?", temp[1]).group()
  r=temp[0]
  t = [r, u]
  t5=time()
  if t == [] or len(t[0]) != len(s):

    # Did not receive expected output from pknotsRG
    return onfail(pkfold, tries, s, T)
  return t[0], float(t[1])


# NEED TO CHANGE THIS - RNAfold CAN'T TAKE PSEUDOKNOTS!!!!
# Compute base pair probabilities and return as a list of dictionaries
# with the None entry giving the unpaired probability

class _ProbVec(dict):
  def __getitem__(self, key):
    if self.has_key(key):
      return super(_ProbVec, self).__getitem__(key)
    else:
      return 0.0
    
def boltzmann(s, T = None, struct = None, tries = 0):
  # Class for maintaining probabilities, if none is present return 0
  #tmpdir = mkdtemp(
  
  global wRNAfold
  if wRNAfold:
      tempid,tempfile=mkstemp()
      inputstring=str(tempfile)
      """Input into RNA fold must be a single string that is 
        delimited by spaces. (there should not be a space at the beginning of the string.
        The file name to write to should come first, followed by running flags.  
        T hen the string 'sequences'  followed by all of the sequences.  Each sequences should be preceded by '>'
        and a unique index for referencing the output back to the group """
      if T==None:
          inputstring+=" -p --noPS sequences"
      else:
          inputstring+=" -p --noPS -T "+str(T)+" sequences"
            
      inputstring+=' >'+(int2str(s))+' '+int2str(s)
      RNAfold.runRNAfold(inputstring)
      f=open(tempfile,'r')
      t=f.read()
      f.close()
      os.remove(tempfile)
  else:
        if T==None and struct == None:
            cmd=[_ViennaRNA_add_location('RNAfold'),'-p','--noPS']          
        elif T==None and struct != None:
            cmd=[_ViennaRNA_add_location('RNAfold'),'-p','--noPS', '-C']
        elif T!=None and struct != None:
            cmd=[_ViennaRNA_add_location('RNAfold'),'-p','--noPS','-C','-T',str(T)]
        else:
            cmd=[_ViennaRNA_add_location('RNAfold'),'-p','--noPS','-T',str(T)]
            
        p=Popen(cmd,stdin=PIPE,stdout=PIPE)
        if struct == None:
            print >>p.stdin,">"+int2str(s)
            print >>p.stdin,int2str(s)
        else:
            print >>p.stdin,">"+int2str(s)
            print >>p.stdin,int2str(s)
            print >>p.stdin, struct
            
        p.stdin.close()
        t=p.stdout.read()
  
  lines=t.split('\n')

  for i in xrange(lines.count('')):
      lines.remove('')

  try:
      energies = lines[2].strip().split(None, 1)
      ensemble = float(lines[3].strip().split(None, 1)[1][1:-1]) 
      if energies == [] or len(energies[0]) != len(s):
          # Did not receive expected output from RNAfold
          raise IndexError
  except IndexError:
      raise IndexError,"Precompute did not receive appropriate output from RNAfold"
            
  p = [_ProbVec() for j in xrange(len(s))]
  for j in p:
      j[None] = 1.0
  #f = open(os.path.join(tmpdir, str(i)+'_dp.ps'))
  #t = f.readline()
  
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

  #rmtree(tempfile,ignore_errors=True)
  return (p, energies[0], float(energies[1][1:-1]), ensemble)


def pboltzmann(s, T = None, struct = None, tries = 0):

  #run PPfold and get the dot  matrix out of the file created
  cmd=['java','-jar','PPfold3.1_special.jar', '--seq', int2str(s)]

  p=Popen(cmd, stdin = PIPE, stdout = PIPE, cwd=os.getcwd())
  t = p.stdout.readlines()
  while p.stdout.readline() != "":
    # Make sure RNAfold has properly finished before reading probability plot
    pass
  p.stdout.close()

  p = [_ProbVec() for i in xrange(len(s))]
  for i in p:
    i[None]=1.0
    
  try:
    struct = t[0].strip("][''\n")
    pv = [] 
    for line in t :
      if re.search("([0-9])", line) != None and re.search("file",line) == None:
        pv.append(map(float,line.split()))

    for i in range(len(pv)):
      for j in range(len(pv[i])):
        if pv[i][j] > float(.0000001):
          q=pv[i][j]
          p[i][j] = q
          p[j][i] = q
          p[i][None] -= q
          p[j][None] -= q
          
    if struct == [] or len(p) == 0:
      # Did not receive expected output from RNAfold
      raise IndexError
  except IndexError:
    return onfail(pboltzmann, tries, s, T)

  cmd2 = [_ViennaRNA_add_location("RNAeval")]
  tmpdir = mkdtemp()
  pop = Popen(cmd2, stdin = PIPE, stdout = PIPE, cwd = tmpdir)
  print >> pop.stdin, int2str(s)
  print >> pop.stdin, struct
  pop.stdin.close()
  u = pop.stdout.readlines()[-1].strip().split(None, 1)
  pop.stdout.close()
  rmtree(tmpdir, ignore_errors = True)
  if u == [] or len(u[0]) != len(s):
    # Did not receive expected output from RNAfold
    return onfail(score, tries, s, t)
  return p, struct, float(u[1][1:-1])

def pkboltzmann(s, T = None, tries = 0):
  (fd,tname) = mkstemp(suffix='.seq')
  tfile = os.fdopen(fd,'w')
  tfile.write(';.seq file')
  tfile.write('Sequence')
  tfile.write(int2str(s) + "1")
  tfile.close()

  cmd = ["pfunction", tname, "output.pfs"]
  p = Popen(cmd, stdin = PIPE, stdout = PIPE, cwd = os.getcwd())
  


# RNAeval works with two-bracket pseudoknots
# Evaluate the energy of structure t on sequence s
def energy(s, t, tries = 0):
  cmd = [_ViennaRNA_add_location("RNAeval")]
  tmpdir = mkdtemp()
  p = Popen(cmd, stdin = PIPE, stdout = PIPE, cwd = tmpdir)
  print >> p.stdin, int2str(s)
  print >> p.stdin, t.bracket()
  p.stdin.close()
  u = p.stdout.readlines()[-1].strip().split(None, 1)
  p.stdout.close()
  rmtree(tmpdir, ignore_errors = True)
  if u == [] or len(u[0]) != len(s):
    # Did not receive expected output from RNAfold
    return onfail(score, tries, s, t)
  return float(u[1][1:-1])


# Function for generating initial sequence for a given structure - randomly
# choose a nucleotide for each left bracket and point, and randomly choose
# suitable pairing nucleotide for each right bracket (pseudoknot structure
# takes the pk... version)

def design(t):
    s = pkpairing(t)
    u = ""
    for i in xrange(len(t)):
        if s[i]:
            if t[i] in "{([<.":
                u += choice("ACGU")
            else:
                if u[s[i]] == "A":
                    u += "U"
                elif u[s[i]] == "C":
                    u += "G"
                elif u[s[i]] == "G":
                    u += choice("CU")
                elif u[s[i]] == "U":
                    u += choice("AG")
        else:
            u += choice("ACGU")
    return u
