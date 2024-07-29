from copy import deepcopy
from types import IntType, ListType, TupleType, StringTypes
### Mapping from integers to bases and base pairs, and from integer
### representations of base pairs to pairs of integer representations
### of bases.
BASES = "ACGUT"
def BASESindex(x):
  return min(BASES.index(x), 3)
  
BASEPAIRS = ["AU", "CG", "GC", "UA", "GU", "UG"]

def int2str(s):
  return "".join(map(lambda x: BASES[int(x)], s))

def str2int(s):
  return map(lambda x: BASESindex(x.upper()), s)

def int2bp(a):
  if a < 4:
    # Canonical base pair
    return (a, 3 - a)
  else:
    # Wobble base pair
    return (a - 2, 7 - a)

def bp2int(a):
  return 2 * a[0] + a[1] - 3

# Base pairs differing between two structures
def bpdif(t, u):
    
  return t.basepairs().symmetric_difference(u.basepairs())

# Base pair distance between two structures
def bpdist(t, u):
  return len(bpdif(t, u))

# Positions where two stuctures differ
def posdif(t, u):
  symdiff=t.basepairs().symmetric_difference(u.basepairs())
  poss=reduce(lambda x,y:x.union(y),symdiff,set([]))                                           
  return poss

# Distance measured as number of positions where structures differ
def posdist(t, u):
  return len(posdif(t, u))

# Return array of pairing position (pk... for pseudoknots), None indicating unpaired position
def pkpairing(s):
  p = len(s) * [None]
  par_stack = []
  cur_stack = []
  squ_stack = []
  ang_stack = []
  for i in xrange(len(s)):
      if s[i] == "(":
          par_stack.append(i)
      elif s[i] == ")":
          j = par_stack.pop()
          p[i] = j
          p[j] = i
      elif s[i] == "{":
          cur_stack.append(i)
      elif s[i] == "}":
          j = cur_stack.pop()
          p[i] = j
          p[j] = i
      elif s[i] == "[":
          squ_stack.append(i)
      elif s[i] == "]":
          j = squ_stack.pop()
          p[i] = j
          p[j] = i
      elif s[i] == "<":
          ang_stack.append(i)
      elif s[i] == ">":
          j = ang_stack.pop()
          p[i] = j
          p[j] = i
  return p

def pairing(s):
    p = len(s) * [None]
    stack = []
    for i in xrange(len(s)):
      if s[i] == '(':
        stack.append(i)
      elif s[i] == ')':
        j = stack.pop()
        p[j] = i
        p[i] = j
    return p

# Functions generating a merge of two bracket notations (pk... for pseudoknots):
#  .(){} : structures agree
#  <>  : . in s, (){} in t
#  {}  : (){} in s, . in t
#  []  : (){} in both s and t, but paired with different positions
#  \/  : (){} in s, )(}{ in t
def pkcombinestructures(s, t):
  sp = pkpairing(s)
  tp = pkpairing(t)
  c = []
  for i in xrange(len(s)):
    if sp[i] == tp[i]:
      # Match between structures
      c.append(s[i])
    elif s[i] == '.':
      # Is unpaired, should have been paired
      c.append({'(': '<', "{": "<", "[": "<", "<": "<", ')': '>', "}": ">", "]": ">", ">": ">"}[t[i]])
    elif s[i] == t[i]:
      # Base paired the right way, but to wrong position
      c.append({'(': '[', "{": "[", "[": "[", "<": "[", ')': ']', "}": "]", "]": "]", ">": "]"}[s[i]])
    elif t[i] == '.':
      # Is base paired, should have been unpaired
      c.append({'(': '{', "{": "{", "[": "{", "<": "{", ')': '}', "}": "}", "]": "}", ">": "}"}[s[i]])
    else:
      # Base paired, but in wrong direction
      c.append({'(': '\\', "{": "\\", "[": "\\", "<": "\\", ')': '/', "}": "/", "]": "/", ">": "/"}[s[i]])
  return "".join(c)

def combinestructures(s, t):
  sp = pairing(s)
  tp = pairing(t)
  c = []
  for i in xrange(len(s)):
    if sp[i] == tp[i]:
      # Match between structures
      c.append(s[i])
    elif s[i] == '.':
      # Is unpaired, should have been paired
      c.append({'(': '<', ')': '>'}[t[i]])
    elif s[i] == t[i]:
      # Base paired the right way, but to wrong position
      c.append({'(': '[', ')': ']'}[s[i]])
    elif t[i] == '.':
      # Is base paired, should have been unpaired
      c.append({'(': '{', ')': '}'}[s[i]])
    else:
      # Base paired, but in wrong direction
      c.append({'(': '\\', ')': '/'}[s[i]])
  return "".join(c)

# Compute list of base pairs from bracket representation of structure in s
def bracket2bp(s):
  par_stack = []
  cur_stack = []
  squ_stack = []
  ang_stack = []
  struc = set([])
  for i in xrange(len(s)):
    if s[i] == '(':
      par_stack.append(i)
    elif s[i] == ')':
      struc.add((par_stack.pop(), i))
    elif s[i] == '{':
      cur_stack.append(i)
    elif s[i] == '}':
      struc.add((cur_stack.pop(), i))
    elif s[i] == '[':
      squ_stack.append(i)
    elif s[i] == ']':
      struc.add((ang_stack.pop(), i))
    elif s[i] == '<':
      ang_stack.append(i)
    elif s[i] == '>':
      struc.add((ang_stack.pop(), i))
  return struc

# Compute list of loops from bracket representation of structure in s,
# where each loop is reported as the list of backbone edges in the
# loop represented by the 5' index of the edge(not used - currently only
# compatible with non-pseudoknotted structures)
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

# Class for maintaining a structure in various formats
class Structure:
  # Function for initialising a structure
  def __init__(self, struct = None, T = None):
    if struct != None:
      if type(struct) in StringTypes:
        # Structure passed as bracket notation
        self.__bracket = struct
      elif type(struct) in [ListType, TupleType] and len(struct) == 2 and type(struct[0]) == IntType and "__iter__" in dir(struct[1]):
        # Structure passed as pair of sequence length and set/list of base pairs
        self.__length = struct[0]
        self.__basepairs = struct[1]
    self.temperature = T
  # Length of structure is number of positions
  def __len__(self):
    if "_Structure__length" not in dir(self):
      self.__length = len(self.bracket())
    return self.__length
  # Generate representation of structure
  def __repr__(self):
    return self.bracket() + int(self.temperature != None) * " {0}".format(self.temperature)
  # Validate that structure information is well formed. If seq is not
  # None, it is assumed to be a sequence and it is also checked that
  # all base pairs are canonical.
  def __getitem__(self,i):
      return self.bracket()[i]
  def validate(self, seq = None):
    # Set up function for checking base pairs
    if seq != None:
      # Sequence to check structure agains has been specified
      if type(seq) in StringTypes:
          
        # Sequence is passed as string, 
        def checkbp(b):
          # Convert characters to upper case and replace T with U
          for i in xrange(len(b)):
            if seq[b[i]] in "tT":
              b[i] = 'U'
            else:
              b[i] = seq[b[i]].upper()
          return "".join(b) in BASEPAIRS
      else:
        # Sequence is passed as list of integers
        def checkbp(b):
          return "".join(map(lambda x: BASES[seq[x]], b)) in BASEPAIRS
    else:
      # No sequence to check structure against
      def checkbp(b):
        return True
    # Check validity of bracket representation
    try:
      if self.__bracket != None:
        count = 0
        for i in xrange(len(self.__bracket)):
          if self.__bracket[i] in '({[<':
            count += 1
          elif self.__bracket[i] in ')}]>':
            count -= 1   
        if count != 0:
          raise ValueError, "Structure " + self.__bracket + " is not well formed"
    except AttributeError:
      # Structure does not have bracket representation
      pass
    # Check validity of basepair representation
    try:
      for p in self.__basepairs:
        for i in p:
          if i < 0 or i >= self.__length:
            # Position i is outside range from 0 to length - 1
            raise ValueError, "Structure " + str(self.__basepairs) + " does not have all base paired positions within range of its length of " + str(self.__length)
        if not checkbp(p):
          if type(seq) not in StringTypes:
            seq = int2str(seq)
          raise ValueError, "Structure " + str(self.__basepairs) + " and sequence " + seq + " are not compatible"
    except AttributeError:
      # Structure does not have basepair representation
      pass

  # Function for returning structure in bracket format - this should
  # be the go-to function for all other structure formats (can handle both
  # knotted and non-knotted structures)
  def bracket(self):
    try:
      if self.__bracket == None:
        raise AttributeError
    except AttributeError:
      if self.__basepairs != None:
        b = self.__length * ['.']
        s = sorted(self.basepairs, key = lambda x: int(x[0]))
        brackets = "({[<)}]>"
        nested_count=0
        b[s[0][0]] = "("
        b[s[0][1]] = ")"
        for i in xrange(1, len(s)):
          if s[i-1][0] < s[i][0] < s[i-1][1] < s[i][1]:
            nested_count += 1
            if nested_count >= 4:
              raise IndexError, "Structure has more than 4-nested pseudoknot"
          elif s[i][0] > s[i-1][1]:
            nested_count = 0
          b[s[i][0]] = brackets[nested_count]
          b[s[i][1]] = brackets[nested_count + 4]
        self.__bracket = "".join(b)
        print self.__bracket
      else:
        # No source of structure information
        raise LookupError, "Attempt to look up structure that is undefined"
    return self.__bracket

  # Function for returning structure as set of base pairs (can handle both
  # knotted and non-knotted strucutures)
  def basepairs(self):
    try:
      if self.__basepairs == None:
        raise AttributeError
    except AttributeError:
      try:
        self.__basepairs = set([])
        par_stack = []
        cur_stack = []
        squ_stack = []
        ang_stack = []
        s = self.bracket()
        self.__length = len(s)
        for i in xrange(len(s)):
          if s[i] == '(':
            par_stack.append(i)
          elif s[i] == ')':
            self.__basepairs.add((par_stack.pop(), i))
          elif s[i] == "{":
            cur_stack.append(i)
          elif s[i] == "}":
            self.__basepairs.add((cur_stack.pop(), i))
          elif s[i] == '[':
            squ_stack.append(i)
          elif s[i] == ']':
            self.__basepairs.add((squ_stack.pop(), i))
          elif s[i] == "<":
            ang_stack.append(i)
          elif s[i] == ">":
            self.__basepairs.add((ang_stack.pop(), i))
        if par_stack != [] or cur_stack != [] or squ_stack != [] or ang_stack != []:
          raise IndexError
      except IndexError:
        raise ValueError, "Invalid bracket representation:" + s
    return self.__basepairs
  # Function for returning unpaired positions in structure - this does
  # not uniquely define the structure and is not considered a valid
  # structure representation.
  def unpaired(self):
    u = set(range(len(self)))
    for b in self.basepairs():
      u.difference_update(b)
    return u

# Class for maintaining target consisting of one or more structures
class Target:
  # Initialisation of target - structures should be list of target structures
  def __init__(self, structures = None):
    if structures == None:
      self.structures = []
    else:
      self.structures = [s for s in structures]
  # Length of a target is number of structures it contains
  def __len__(self):
    return len(self.structures)
  # Generate representation of target
  def __repr__(self):
    return "\n".join([repr(x) for x in self.structures])
  # Size of a target is number of positions it contains
  def size(self):
    if self.structures == []:
      return -1
    else:
      return len(self.structures[0])
  # Validate that target consists of well formed structures of the same length
  def validate(self, seq = None):
    for s in self.structures:
      s.validate(seq)
      if len(s) != len(self.structures[0]):
        raise ValueError, "Target structures not of the same length"
 
  def __getitem__(self,i):
      return self.structures[i]
  # Iterator for target structures
  def __iter__(self):
    if self.structures != None:
      for s in self.structures:
        yield s
  # Add one more structure to target
  def addstructure(self, structure):
    self.structures.append(structure)
    # Previously computed dependencies and cuts are now void
    self.dependencies = None
  # Compute the dependency classes imposed by the structures in this target
  def computedependencies(self):
    if len(self.structures) == 0:
      raise RuntimeError, "Cannot build dependency structure with no targets"
    n = len(self.structures[0].bracket())
    self.dependencies = [set([]) for i in xrange(n)]
    for t in self.structures:
      for bp in t.basepairs():
        self.dependencies[bp[0]].add(bp[1])
        self.dependencies[bp[1]].add(bp[0])
    
    # Check consistency and separate into components
    self.components = []
    colour = n * [False]
    for i in xrange(n):
      if colour[i] is False:
        # Node i is not part of previously visited component; assign
        # it first colour and start building a new component
        stack = [i]
        colour[i] = 1
        self.components.append(set([]))
        # Iterate through elements in component
        while stack != []:
          j = stack.pop()
          self.components[-1].add(j)
          # Check neighbours to this element
          for k in self.dependencies[j]:
            if colour[k] is False:
              # First time we see this element, colour it with
              # opposite colour of current element and insert it as
              # one that has to be visited
              colour[k] = 3 - colour[j]
              stack.append(k)
            else:
              # We have seen this element before, make sure that its
              # colour is consistent with current element
              if colour[k] == colour[j]:
                raise ValueError, "Targets impose incompatible restrictions on sequence"
    
  # Retrieve dependency component containing position i
  def getdependency(self, i):
    if "dependencies" not in dir(self) or self.dependencies == None:
      self.computedependencies()
      print "getdependency", self.dependencies
    return self.dependencies[i]
  # Iterator for all dependency components
  def getdependencies(self):
    if "components" not in dir(self) or self.dependencies == None:
      self.computedependencies()
    return self.components
  # Retrieve all crossover cuts compatible with the dependencies. A
  # cut is given by the positions before which the sequence changes,
  # and is represented as sets from which any two elements constitute
  # a cut compatible with the dependencies.
  def getcuts(self):
    if "cuts" not in dir(self) or self.dependencies == None:
      # Make sure we have some dependencies to compute the cuts from
      if "dependencies" not in dir(self) or self.dependencies == None:
        self.computedependencies()
      # Build cut classes - initially every position is in the same class
      self.cuts = set([tuple(range(len(self.dependencies) + 1))])
      # Each dependency break classes that have positions both inside
      # and outside the base pair
      for i in xrange(len(self.dependencies)):
        for j in self.dependencies[i]:
          if i < j: # Each base pair is registered in two places
            # We can't delete and add to classes while we iterate
            # through its members, so store information about classes
            # to be deleted and to be added
            new = []
            old = []
            for A in self.cuts:
              # For each class, determine positions inside and outside
              # current dependency
              inside = []
              outside = []
              for k in A:
                if k > i and k <= j:
                  inside.append(k)
                else:
                  outside.append(k)
              if len(inside) > 0 and len(outside) > 0:
                # Class split into two classes
                old.append(A)
                if len(inside) > 1:
                  # Classes reduced to a single position do not define
                  # a non-empty set of valid cuts.
                  new.append(tuple(inside))
                if len(outside) > 1:
                  # Classes reduced to a single position do not define
                  # a non-empty set of valid cuts.
                  new.append(tuple(outside))
            # Update cut classes with the new splits
            for A in old:
              self.cuts.remove(A)
            for A in new:
              self.cuts.add(A)
    # Now we know we have cut classes
    return self.cuts
