from subprocess import Popen, PIPE
from tempfile import mkdtemp
from shutil import rmtree
from os.path import join
from structure import int2str, Structure
from time import sleep

# Function defining behaviour if a folding program fails
def onfail(f, tries, *args):
  if tries < 10:
    # Don't retry more than 10 times
    sleep(10) # Take ten
    return f(*args, tries = tries + 1)
  else:
    raise RuntimeError, "Function " + f.__name__ + " did not receive result from external program"

# Compute list of base pairs from bracket representation of structure in s
def bracket2bp(s):
  stack = []
  struc = set([])
  for i in xrange(len(s)):
    if s[i] == '(':
      stack.append(i)
    elif s[i] == ')':
      struc.append((stack.pop(), i))
  return struc

# Compute list of loops from bracket representation of structure in s,
# where each loop is reported as the list of backbone edges in the
# loop represented by the 5' index of the edge
def bracket2loops(s):
  stack = [[-1]]
  loops = []
  for i in xrange(len(s)):
    if s[i] == '(':
      stack.append([i])
    else:
      if s[i] == ')':
        loops.append(stack.pop())
      stack[-1].append(i)
  if len(stack[-1]) > 1:
    loops.append(stack.pop())
  return loops

# Compute MFE structure of RNA sequence s, return pair of structure in
# bracket notation and MFE as float
def fold(s, T = None, tries = 0):
  cmd = ["RNAfold"]
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
  return t[0], float(t[1][1:-1])

# Compute base pair probabilities and return as a list of dictionaries
# with the None entry giving the unpaired probability
def boltzmann(s, T = None, tries = 0):
  # Class for maintaining probabilities, if none is present return 0
  class _ProbVec(dict):
    def __getitem__(self, key):
      if self.has_key(key):
        return super(_ProbVec, self).__getitem__(key)
      else:
        return 0.0
  cmd = ["RNAfold", "-p"]
  tmpdir = mkdtemp()
  if T != None:
    cmd.extend(["-T", str(T)])
  p = Popen(cmd, stdin = PIPE, stdout = PIPE, cwd = tmpdir)
  print >> p.stdin, int2str(s)
  p.stdin.close()
  mfe = p.stdout.readlines()[1].strip().split(None, 1)
  while p.stdout.readline() != "":
    pass
  p.stdout.close()
  if mfe == [] or len(mfe[0]) != len(s):
    # Did not receive expected output from RNAfold
    return onfail(boltzmann, tries, s, T)
  p = [_ProbVec() for i in xrange(len(s))]
  for i in p:
    i[None] = 1.0
  f = open(join(tmpdir, "dot.ps"))
  t = f.readline()
  while "data starts here" not in t:
    t = f.readline()
  t = f.readline()
  while t != "" and "showpage" not in t:
    if "ubox" in t:
      t = t.split()
      i = int(t[0]) - 1
      j = int(t[1]) - 1
      q = pow(float(t[2]), 2)
      p[i][j] = q
      p[j][i] = q
      p[i][None] -= q
      p[j][None] -= q
    t = f.readline()
  f.close()
  rmtree(tmpdir, ignore_errors = True)
  return (p, mfe[0], float(mfe[1][1:-1]))

# Evaluate the energy of structure t on sequence s
def energy(s, t, tries = 0):
  cmd = ["RNAeval"]
  tmpdir = mkdtemp
  if t.temperature != None:
    cmd.extend(["-T", str(t.temperature)])
  p = Popen(cmd, stdin = PIPE, stdout = PIPE, cwd = tmpdir)
  print >> p.stdin, int2str(s)
  print >> p.stdin, t.bracket()
  p.stdin.close()
  u = p.stdout.readlines()[-1].strip().split(None, 1)
  p.stdout.close()
  rmtree(tmpdir, ignore_errors = True)
  if u == [] or len(u[0]) != len(s):
    # Did not receive expected output from RNAfold
    return onfail(energy, tries, s, t)
  return float(u[1][1:-1])
