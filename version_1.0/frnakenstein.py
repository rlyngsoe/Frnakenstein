#!/home/stat-compbio-summer-school/cbss0030/Programs/bin/python
#This is the svn version on the server
from sys import stdout, platform
from individual import BaseDistribution, Individual
from types import IntType, TupleType, StringType
import multiprocessing
from argparse import ArgumentParser, ArgumentError, FileType, Action, SUPPRESS
from sys import stdin, stderr
from individual import Individual, RandomInitialise, DesignInitialise, FileInitialise, RandomPosition, WrongPosition, ProbabilityPosition, MinProbabilityPosition, LogProbabilityPosition, ThresholdProbabilityPosition, RelativeProbabilityPosition, CombinedPosition, RandomCut, CorrectCut, ProbabilityCut, MinProbabilityCut, ProductProbabilityCut, ThresholdProbabilityCut, RelativeProbabilityCut, CombinedCut, MutateDependent, WeightedMutateDependent, RetainMutateDependent, WrongObjective, ProbabilityObjective, MinProbabilityObjective, ProductProbabilityObjective, ThresholdProbabilityObjective, RelativeProbabilityObjective, EnsembleProbabilityObjective, ConjunctiveObjective, DisjunctiveObjective, WrongFitness, ProbabilityFitness, MinProbabilityFitness, ProductProbabilityFitness, ThresholdProbabilityFitness, RelativeProbabilityFitness, EnsembleProbabilityFitness, CombinedFitness
from population import Population, Mutate, RandomMutate, FitnessWeightedMutate, SelectPair, WeightedSelectPair, CombinedSelectPair, SelectCut, WeightedSelectCut, Reduce, DiversityReduce, ParetoReduce, Recombine
from structure import Structure, Target, int2str, combinestructures
from types import FunctionType, LambdaType, ListType, ClassType
from time import time


# Function outputting the members of population pop sorted by
# objective. Output is written to file, and the function format is
# used to format each individual. The format function takes five
# arguments: the sequence of the individual, the objective function
# value for the individual, the fit function value for the
# individual, the predicted structure for the individual, and the ID
# for the individual.
def output_population(pop, file = stdout, format = None):
  if format == None:
    # Set up default formatting
    def format(seq, obj, fit, struc, iid):
      f = "Sequence: {0}\n Objective: {1}, Fitness: {2}\n Structure{3}".format(seq, obj, fit, int(len(struc) > 1) * "s")
      for s in struc:
        if s[0] != None:
          f += "\n  {0}, T = {1}".format(s[0], s[1])
        else:
          f += "\n  {0}".format(s[0])
      if iid != None:
        # Parse ID tree and find summary statistics: first element is
        # set of base sequences involved in the recombination history
        # of sequence (note: it is not checked whether their
        # contribution has survived until the end); second element is
        # number of mutations; third element is number of
        # recombinations.
        def parse_id(i):
          if type(i) == IntType or type(i) == StringType:
            return (set([i]), 0, 0)
          elif type(i) == TupleType:
            if len(i) == 1:
              return (lambda x: (x[0], x[1] + 1, x[2]))(parse_id(i[0]))
            elif len(i) == 2:
              return (lambda x, y: (x[0].union(y[0]), x[1] + y[1], x[2] + y[2] + 1))(parse_id(i[0]), parse_id(i[1]))
          raise ValueError, "Could not parse ID"
        try:
          s = parse_id(iid)
          f += "\n # base seqs: {0}, # mutations: {1}, # recombinations {2}".format(len(s[0]), s[1], s[2])
        except ValueError:
          f += "\n ID: {0}".format(i)
      return f

  # Create list of predicted structures at all target temperatures
  target = pop.gettarget()
  temps = set([])
  for t in target:
    temps.add(t.temperature)
  def getstructures(ind):
    tseq = []
    for t in temps:
      tseq.append((ind.getstruct(t).bracket(), t))
    return tseq
  
  
  m = [(i.getobjective(), i.getfitness(), i.getseq(), getstructures(i), i.getid(), i.getgccontent()) for i in pop.members]
  m.sort()
  for i in m:
    print >> file, format(i[2], i[0], i[1], i[3], i[4])

# Function running a GA of ngen generation to construct sequences
# folding into target, which should be a Target object. Population
# size is fixed at npop, in each generation nmut mutants are created,
# nxov crossovers are performed, and nrandom new random sequences are
# added. Population and Individual are classes for representing the
# population and an individual. For every generation the current
# population is written to logfile if this is not None. When
# nucleotides need to be chosen randomly, specifically for mutations
# and possibly for initialisation of sequences, these are chosen
# according to basedist, which can either be a BaseDistribution object
# or an object or structure with attributes or keys with names in the
# set {unpaired, basepair, paired} which hold vectors of 4, 6, and 4
# probabilities, respectively, used for drawing unpaired nucleotides,
# base pairs, and nucleotides in paired positions. If not None,
# stoppingcriteria is applied to an iterator over the population
# members at the end of each generation, and if the return value is
# interpreted as True the GA is terminated. If report is not None,
# this is called with population and generation number at the end of
# each generation, uses possibly including output of progress
# information or update of information in GUI. If track is interpreted
# as True, information is stored in individuals that allow tracking of
# their ancestry. If precompute is not None, it is applied to the set
# of new individuals created in each generation, e.g. allowing
# computation on each individual to be carried out in parallel.

def ga (npop, ngen, nmut, nxov, pop, logfile=None, stoppingcriteria=None, nrandom=None, report=None, track=False, usepprocess=False):
  # Set up distribution for drawing nucleotides
  if track:
    stop = pop.addrandom(npop, idfunc = (lambda x: str(x)), stoppingcriteria = stoppingcriteria)
  else:
    stop = pop.addrandom(npop, stoppingcriteria = stoppingcriteria)
  
  avgtimepre=0
  avgptime=0

  # Precompute in parallel or on single core
  if usepprocess:
      new=pop.massprecompute(pop.new)
      if new!=None:
          pop.new=new
  else:
      pop.precompute(pop.new)

  pop.addnew()
  
  
  
  if stop:
    # Stopping criteria met during initialisation
    if logfile != None:
      output_population(pop, logfile)
      print >> logfile, "Stopping criteria met after adding {0} initial sequences to population".format(len(pop.members))
    return pop
  # Evolve population for ngen generations
  avgtimec=0
  avgtimem=0
  avgtimered=0
  avgtimepre=0
  for g in xrange(ngen):
    print '*****',g,'*****'
    # Create nmut mutated sequences
    if usepprocess:
        t1=time()
        pop.pmutate(nmut)
        t2=time()
        avgtimem+=(t2-t1)
        
        # Create nxov sequences generated by crossover
        t1=time()
        pop.precombine(nxov)
        t2=time()
        avgtimec+=(t2-t1)
    else:
        t1=time()
        pop.mutate(nmut)
        t2=time()
        avgtimem+=(t2-t1)
        t1=time()
        pop.recombine(nxov)
        t2=time()
        avgtimec+=(t2-t1)
    # Add newly created sequences to population
    t1=time()
    if usepprocess:
      new=pop.massprecompute(pop.new)
      #avgptime+=ptime
      if new!=None:
          pop.new=new
    else:
      pop.precompute(pop.new)
    t2=time()
    avgtimepre+=(t2-t1)
    pop.addnew()
    
    
    # Random sequences should not be considered for elimination until
    # next round, but need to be part of the precomputation
    if nrandom and g != ngen - 1:
      # Add random sequences to maintain diversity
      if track:
        pop.addrandom(nrandom, idfunc = (lambda x: str(g) + "." + str(x)))
      else:
        pop.addrandom(nrandom)
        
    # Precompute structural information (can be used to parallellise
    # computations or prevent initial computation of MFE structure
    # followed by Boltzmann computation that also computes MFE
    # structure).
    t1=time()
    if usepprocess:
      new=pop.massprecompute(pop.new)
      #avgptime+=ptime
      if new!=None:
          pop.new=new
    else:
       pop.precompute(pop.new)
    
    t2=time()
    avgtimepre+=(t2-t1)
    
    # Reduce population to fixed size
    t1=time()
    pop.reduce(npop)
    t2=time()
    avgtimered+=(t2-t1)
    
    
    if logfile != None:
      # Output current status to log file
      print >> logfile, "=== Population after generation {0} ===".format(g)
      print >> logfile, "Stopping criteria met after generation", g
      print >> logfile, 'total time mutating: ', avgtimem
      print >> logfile, 'average time mutating: ',avgtimem/ngen
      print >> logfile, 'total time recombining: ', avgtimec
      print >> logfile, 'average time recombining: ',avgtimec/ngen
      print >> logfile, 'total time precomputing:',avgtimepre
      print >> logfile, 'total time reducing:',avgtimered  
      output_population(pop, logfile)

    if report != None:
      report(pop, g)

    if stoppingcriteria != None and stoppingcriteria(pop.members.__iter__()):
      # Stop genetic algorithm
      if logfile != None:
        print >> logfile, "Stopping criteria met after generation", g
        print >> logfile, 'total time mutating: ', avgtimem
        print >> logfile, 'average time mutating: ',avgtimem/ngen
        print >> logfile, 'total time recombining: ', avgtimec
        print >> logfile, 'average time recombining: ',avgtimec/ngen
        print >> logfile, 'total time precomputing:',avgtimepre
        print >> logfile, 'total time reducing:',avgtimered  
      break

    # Promote random sequences to full members
    pop.addnew()
    
  print 'total time mutating: ', avgtimem
  print 'average time mutating: ',avgtimem/ngen
  print 'total time recombining: ', avgtimec
  print 'average time recombining: ',avgtimec/ngen
  print 'total time precomputing:',avgtimepre
  print 'total time reducing:',avgtimered 
   
  return pop


  # Set up command line interface for running GA

def add_structure(namespace, s, T=None):
    # Function for transforming to a canonical dot-bracket notation, converting
    # all brackets to parentheses if the pseudoknotted option is not specified
    def transform(s):
      global ChangeSoftware
      ChangeSoftware = False
      if namespace.softtype.lower() != "pknotsrg":
        t = ""
        bracket_types = [0,0,0,0]
        for i in s:
          if i == "(":
            t += '('
            bracket_types[0] = 1
          elif i == "{":
            t += '{'
            bracket_types[1] = 1
          elif i == "[":
            t += '['
            bracket_types[2] = 1
          elif i == "<":
            t += '<'
            bracket_types[3] = 1
          elif i == ".":
            t += "."
          elif i in ")>}]":
            t += i
          elif not i.isspace():
            raise ValueError, "Unexpected character " + i + "in target specification"
        if sum(bracket_types) > 1:
          ChangeSoftware = True
          print "Pseudoknot detected, folding with pknotsRG. If structure does not contain a pseudoknot, please use dot-parenthesis format for target structure input"
          options.softtype == "pknotsRG"
      if namespace.softtype.lower() == "pknotsrg":
        t = ""
        for i in s:
          if i in ".-_:":
            t += '.'
          elif i in "({[<)}]>":
            t += i
          elif not i.isspace():
            raise ValueError, "Unexpected character " + i + "in target specification"
      return t

    struc = Structure(transform(s), T)
    try:
      struc.validate()
    except ValueError:
      raise ValueError, "Structure " + str(s) + " not well formed"
    try:
      if namespace.target.size() != len(struc):
        raise ValueError, "Structure " + str(s) + " not of the same length as previous target structure(s)"
    except AttributeError:
      namespace.target = Target()
    namespace.target.addstructure(struc)
    
def set_defaults(options):
    """Default when target is a single structure corresponds to -g N --designinit --retaindependent --mutatefit --combinedposition --wrongposition 2.0 --boltzmannposition 1.0 --combinedxover --correctxover 1.0 --boltzmannxover 1.0 --combinedpairs --wrongobjective --combinedfitness --wrongfitness 1.0 --boltzmannfitness 1.0 --diversity 1.0. Default when target consists of multiple structures corresponds to -g N --randominit --retaindependent --mutatefit --combinedposition --wrongposition 1.0 --boltzmannposition 1.0 --logboltzmannposition 1.0 --combinedxover --correctxover 1.0 --boltzmannxover 1.0 --prodboltzmannxover 2.0 --combinedpairs --wrongobjective --combinedfitness --wrongfitness 1.0 --boltzmannfitness 1.0 --prodboltzmannfitness 2.0 --diversity 1.0, except if the target is meta-stable, i.e. two or more structures are specified for the same temperature, where --minprobabilityobjective T --ensemblefitness 1.0 replaces the objective and fit options. N is the product of the length of the target and the number of target structures, and T is two thirds of the reciprocal of the number of target structures. When only some aspects are specified, these default values take effect for the remaining aspects, e.g. if the only command line option is -g 50 the default choices are still used for mutations, recombinations, etc.; the --diversity option is only added when no fit scheme has been specified, so when specifying a fit scheme you do not have to add --diversity 0 to cancel the diversity augmentation. If no target structure is specified, a simple input interface queries the user for one or more target structures."""

    # Output s with line length l (default 78) to outfile (default stdout)
    def format_output(s, l = 78, outfile = stdout):
      while s != "":
        if len(s) <= l:
          # Not enough left for multiple lines
          print >> outfile, s
          break
        if s[l].isspace():
          # Natural break point after first l characters
          print >> outfile, s[:l]
          s = s[l + 1:].strip()
        else:
          # Find latest natural breakpoint among first l characters
          i = l + 1 - len(s[:l + 1].split()[-1])
          if i > 0:
            print >> outfile, s[:i - 1]
            s = s[i:].strip()
          else:
            print >> outfile, s[:l]
            s = s[l:].strip()
    # Check whether target structure has been specified
    if "target" not in options:
      format_output("Input target structure(s) in dot-bracket notation. Each structure can be entered over several lines. The end of a structure is indicated by either an empty line or the presence of a numerical value, separated by a white space character from dot-brackets, at the end of a line. In the latter case, the numerical value is interpreted as the temperature specification for the target structure. Terminate input with a single @ at the end of the line, or on an input line of its own.")
      cont = True
      while cont:
        # Input structures
        # Output ruler
        for i in xrange(1, 9):
          stdout.write(4 * '.' + ',' + 4 * '.' + str(i))
        stdout.write("\n")
        t = []
        temp = None
        while cont:
          # Input lines specifying current structure
          u = stdin.readline()
          if len(u) <= 1:
            if len(u) == 0:
              # EOF
              cont = False
            break
          v = u.strip()
          if v != "":
            # Got a live one
            if v[-1] == '@':
              # End of input
              cont = False
              v = v[:-1].split()
            else:
              v = v.split()
            if v != []:
              if v[-1][-1].isdigit():
                # Temperature specification at end of line
                try:
                  temp = float(v[-1])
                  if temp == int(temp):
                    temp = int(temp)
                except ValueError:
                  raise ValueError, "Invalid temperature specification: " + v[-1]
                
                t.extend(v[:-1])
                break
              else:
                t.extend(v)

    
        # Reached end of structure specification
        t = "".join(t)
        if t != "" or temp != None:
          add_structure(options, t, temp)
      if ChangeSoftware:
        options.softtype = "pknotsRG"
      format_output("Structure{0} specified:".format(int(len(options.target) == 1) * "s"))
      print options.target

    # Check what kind of target has been specified
    temps = set([])
    metastable = False
    for t in options.target:
      if t.temperature in temps:
        metastable = True
        break
      temps.add(t.temperature)
    # Set default values for options not specified
    if "g" not in options:
      # Set number of generations to number of target structures times
      # length of target
      options.__setattr__("g", len(options.target) * options.target.size())
      if options.v >= 2:
        format_output("Number of generations set to {0}".format(options.g))
    
    if options.motif:
        options.__setattr__("initialise", RandomInitialise)
        
    if "initialise" not in options:
      if len(options.target) > 1 and not options.motif:
        options.__setattr__("initialise", RandomInitialise)
        if options.v >= 2:
          format_output("Using random compatible sequences for initialisation")
      else:
        options.__setattr__("initialise", DesignInitialise)
        if options.v >= 2:
          format_output("Using RNAinverse designs for initialisation")
    if "wdependent" not in options:
      # Set dependency induced mutations to depend on positional weight
      options.__setattr__("wdependent", RetainMutateDependent)
      if options.v >= 2:
        format_output("Mutation propagation scheme set to retaindependent")
    if "mutate" not in options:
      # Set selection of sequence for mutation to be based on fit
      options.__setattr__("mutate", FitnessWeightedMutate)
      if options.v >= 2:
        format_output("Scheme for selecting sequences for mutations set to mutatefit")
    if "position" not in options:
      # Set selection of position for mutation to be based on weighted
      # combinations of whether predicted structure matches target
      # structure, Boltzmann probability of incorrect structure, and
      # if multiple structures are specified the sum of negative
      # logarithms of probability of correct structure.
      if options.softtype == "pknotsRG":
        options.__setattr__("position", [(2.0, WrongPosition), (1.0, RandomPosition)])
        if options.v >= 2:
          format_output("Scheme for selecting position to mutate set to 2-1 weighted combination of wrongposition and randomposition")
      else:
        if len(options.target) > 1:
          # More than one target structure
          options.__setattr__("position", [(1.0, WrongPosition), (1.0, ProbabilityPosition), (1.0, LogProbabilityPosition)])
          if options.v >= 2:
            format_output("Scheme for selecting position to mutate set to 1-1-1 weighted combination of wrongposition, boltzmannposition, and logboltzmannposition")
        else:
          # Single target structure
          options.__setattr__("position", [(2.0, WrongPosition), (1.0, ProbabilityPosition)])
          if options.v >= 2:
            format_output("Scheme for selecting position to mutate set to 2-1 weighted combination of wrongposition and boltzmannposition")
      options.__setattr__("combinedposition", CombinedPosition)
    if "xover" not in options:
      # Set selection of cross over points for recombination to be based
      # on weighted combinations of whether predicted structure matches
      # target structure, Boltzmann probability of incorrect structure,
      # and if multiple target structures are specified the product of
      # probabilities of correct structure.
      if options.softtype == "pknotsRG":
        options.__setattr__("xover", [(2.0, (CorrectCut, WeightedSelectCut)), (1.0, (RandomCut, WeightedSelectCut))])
        if options.v >= 2:
          format_output("Scheme for selecting cross over point set to 2-1 weighted combination of correctcut and randomcut")
      else:
        if len(options.target) > 1:
          # More than one target structure
          options.__setattr__("xover", [(1.0, (CorrectCut, WeightedSelectCut)), (1.0, (ProbabilityCut, WeightedSelectCut)), (2.0, (ProductProbabilityCut, WeightedSelectCut))])
          if options.v >= 2:
            format_output("Scheme for selecting cross over point set to 1-1-2 weighted combination of correctxover, boltzmannxover, and prodboltzmannxover")
        else:
          # Single target structure
          options.__setattr__("xover", [(1.0, (CorrectCut, WeightedSelectCut)), (1.0, (ProbabilityCut, WeightedSelectCut))])
          if options.v >= 2:
            format_output("Scheme for selecting cross over point set to 1-1 weighted combination of correctxover and boltzmannxover")
      options.__setattr__("combinedxover", (CombinedCut, WeightedSelectCut))
    if "pairs" not in options:
      # Set selection of pair for recombination to be based on how well
      # their position fitnesses recombine
      options.__setattr__("pairs", CombinedSelectPair)
      if options.v >= 2:
        format_output("Scheme for selecting pairs of sequences for recombination set to combinedpairs")
    if "objective" not in options:
      # Set objective to be number of wrong positions for predicted
      # structure, except in the case of a metastable target where each
      # position just needs to be correct with probability two thirds of
      # the reciprocal of number of target structures.
      if options.softtype == "pknotsRG":
        options.__setattr__("objective", WrongObjective)
      else:
        if metastable:
          # More than one target at same temperature
          options.__setattr__("objective", MinProbabilityObjective(1 - .67 / len(options.target)))
          if options.v >= 2:
            format_output("Objective function set to minboltzmannobjective {0}".format(1 - .67 / len(options.target)))
        else:
          # No repeated temperature, perfect solution theoretically possible
          options.__setattr__("objective", WrongObjective)
          if options.v >= 2:
            format_output("Objective function set to wrongobjective")
    if "fit" not in options:
      # Set fit to be a weighted combination of whether predicted
      # structure matches target structure, Boltzmann probability of
      # incorrect structure, and if a metastable sequence is specified
      # the product of probabilities of correct structure.
      if options.softtype == "pknotsRG":
        options.__setattr__("fit", WrongFitness)
      else:
        if metastable and options.softtype != "PPfold":
          # Metastable sequence specified
          options.__setattr__("fit", EnsembleProbabilityFitness(1.0))
          if options.v >= 2:
            format_output("Fitness function set to ensemblefitness 1")
        else:
          if len(options.target) > 1:
            # More than one target structure
            options.__setattr__("fit", [(1.0, WrongFitness), (1.0, ProbabilityFitness), (2.0, ProductProbabilityFitness)])
            if options.v >= 2:
              format_output("Fitness function set to 1-1-2 weighted combination of wrongfitness, boltzmannfitness, and prodboltzmannfitness")
          else:
            # Single target structure
            options.__setattr__("fit", [(1.0, WrongFitness), (1.0, ProbabilityFitness)])
            if options.v >= 2:
              format_output("Fitness function set to 1-1 combination of wrongfitness and boltzmannfitness")
          options.__setattr__("combinedfitness", CombinedFitness)
        if "diversity" not in options:
          options.__setattr__("diversity", 1.0)
          if options.v >= 2:
            format_output("Augmenting fit with diversity 1")

class TargetAction(Action):
    def __call__(self, parser, namespace, values, option_string = None):
      if len(values) < 1:
        raise ValueError, "Missing argument for target specification"
      if len(values) > 2:
        raise ValueError, "Too many arguments in target specification"

      if len(values) == 2:
        # Determine temperature
        try:
          T = float(values[1])
          if int(T) == T:
            T = int(T)
        except ValueError:
          raise ValueError, "Temperature in target specification should be a real or integer value"
      else:
        T = None
      add_structure(namespace, values[0], T)
      
# Function generating Action subclasses for parsing possibly
  # weighted options. The value is appended to list in the namespace
  # attribute dest of the Action subclass with the specified weight,
  # if weighted, otherwise the attribute dest of Action is set to
  # value.
def WeightedAction(dest, cls):
    class _WeightedAction(Action):
      def __call__(self, parser, namespace, values, option_string = None):
        # Function for retrieving class when it is really a class
        # generating function.
        def _getcls(cls, values):
          if type(cls) in (FunctionType, LambdaType):
            return cls(values)
          elif type(cls) in (TupleType, ListType):
            newcls = []
            for c in cls:
              if type(c) in (FunctionType, LambdaType):
                newcls.append(c(values))
              else:
                newcls.append(c)
            return newcls
          else:
            return cls
        # cls is combined and should be weighted iff the destination is a list
        if type(values) != ListType:
          if values == None:
            values = []
          else:
            values = [values]
        if dest in namespace and type(namespace.__getattribute__(dest)) == ListType:
          if len(values) == 0:
            raise ValueError, "Missing weight for {0}class {1}".format(int(option_string != None) * ("option " + option_string + "/"), cls.__doc__)
          else:
            namespace.__getattribute__(dest).append((float(values[0]), _getcls(cls, values[1:])))
        else:
          namespace.__setattr__(dest, _getcls(cls, values))
    return _WeightedAction

# Function generating Action subclasses for parsing combination
  # options. ccls should be a function taking a list of schemes as
  # single argument and returning the combination class, dest should
  # be the namespace attribute name for storing the combination class,
  # and subdest the namespace attribute name for storing the classes
  # that are combined which is initialised to an empty list
def CombineAction(dest, ccls, subdest):
    class _CombineAction(Action):
      def __call__(self, parser, namespace, values, option_string = None):
        namespace.__setattr__(dest, ccls)
        namespace.__setattr__(subdest, [])
    return _CombineAction

# Generates function ensuring absence of arguments before returning cls
def noargs(cls):
    def _noargs(values):
      if len(values) > 0:
        if type(cls) == ClassType:
          raise ValueError, "Expected no arguments for option specifying class {0} but got {1}".format(cls.__name__, values)
        else:
          raise ValueError, "Expected no arguments but got {0}".format(values)
      return cls
    _noargs.__doc__ = cls.__name__
    return _noargs

# Generate function ensuring the presence of two numerical
  # arguments in the ranges rang1 and rang2 (default from 0 to 1). The
  # cls argument should be a binary function creating a class from the
  # numerical arguments. opt should be a pair, where a value different
  # from False is taken as indication of an optional argument and used
  # as default. If both arguments are optional and a single value
  # available, this will be assigned to the first argument.
def numnumarg(cls, rang1 = (0, 1), rang2 = (0, 1), opt = (False, False)):
    def _clsname(c):
      if type(c) in (FunctionType, LambdaType):
        return "specifying class {0} ".format(c(max(rang1[0], 0), max(rang2[0], 0)).__name__)
      else:
        return ""
    if type(opt) not in [ListType, TupleType] or len(opt) != 2:
      # opt not in recognised format
      opt = (False, False)
    def _numnumarg(values):
      # Determine number of optional arguments that have been provided
      l = reduce(lambda x, y: x + int(y is False), opt, 0)
      o = len(values) - l
      if len(values) > 2:
        raise ValueError, "Expected (at most) two arguments for option {0} but got {1}".format(_clsname(cls), len(values))
      elif o < 0:
        raise ValueError, "Expected (at least) {0} arguments for option {1} but got {2}".format(l, _clsname(cls), len(values))
      # Build array of arguments
      a = []
      rang = [rang1, rang2]
      i = 0
      for j in xrange(2):
        if opt[j] is False or o > 0:
          # Use value provided for argument
          try:
            p = float(values[i])
            i += 1
          except ValueError:
            raise ValueError, "Expected numerical argument for option {0} but got {1}".format(_clsname(cls), values[i])
          if rang[j][0] != None and p < rang[j][0] or rang[j][1] != None and p > rang[j][0]:
            raise ValueError, "Expected value in range {0} to {1} for option {2} but got {3}".format(rang[j][0], rang[j][1], _clsname(cls), p)
          a.append(p)
          o -= 1
        else:
          # Use default value for optional argument
          a.append(opt[j])
      # Generate class from arguments
      return cls(*a)
    _numnumarg.__doc__ = _clsname(cls)
    return _numnumarg

# Generate function ensuring the presence of a single numerical
  # argument in the range rang (default from 0 to 1). The cls argument
  # should be a unary function creating a class from the numerical
  # argument. If opt is not False, the argument is assumed optional
  # with the value of opt as default.
def numarg(cls, rang = (0, 1), opt = False):
    def _clsname(c):
      if type(c) in (FunctionType, LambdaType):
        return "specifying class {0} ".format(c(max(rang[0], 0)).__name__)
      else:
        return ""
    def _numarg(values):
      try:
        if len(values) != 1:
          if len(values) == 0 and opt is not False:
            p = float(opt)
          else:
            raise ValueError
        else:
          p = float(values[0])
          if (rang[0] != None and p < rang[0]) or (rang[1] != None and p > rang[1]):
            raise ValueError, "Argument ({0}) to option {1}should be between {2} and {3}".format(p, _clsname(cls), rang[0], rang[1])
        return cls(p)
      except ValueError:
        if len(values) == 0 and opt is False:
          raise ValueError, "Expected exactly one numerical argument in the range from 0 to 1 for option {0} but got none".format(_clsname(cls))
        elif len(values) > 1:
          raise ValueError, "Expected exactly one numerical argument in the range from 0 to 1 for option {0} but got {1}".format(_clsname(cls), values)
        else:
          raise ValueError, "Expected numerical argument to option {0} but got {1}".format(_clsname(cls), values[0])
    _numarg.__doc__ = _clsname(cls)
    return _numarg
  

def createParser():
      parser = ArgumentParser(description = "Runs a genetic algorithm attempting to create RNA sequences folding into the specified secondary structure targets. Read usage descriptions for arguments taking a variable number of parameters carefully, as there may be restrictions not apparent from the basic specification; this is due to limitations in the module used for parsing the command line arguments.", epilog = set_defaults.__doc__, argument_default = SUPPRESS)
      parser.add_argument("--version", action = 'version', version = "%(prog)s 2.25/07/12")  
      parser.add_argument("-t", action=TargetAction, nargs="+", default=None, help="Add structure S as target structure. If a temperature T is specified, the target is assumed to be for this temperature, otherwise folding at default temperature is performed. Takes one mandatory and one optional argument.", metavar=("S", "T"))
      parser.add_argument("-g", type=int, help="Run Genetic Algorithm for G generations")
      parser.add_argument("-s", default=50, type=int, help="Set population size to N", metavar="N")
      parser.add_argument("-m", default=50, type=int, help="Create N mutants in each generation", metavar="N")
      parser.add_argument("-x", default=50, type=int, help="Create N recombinants in each generation", metavar="N")
      cpu = multiprocessing.cpu_count()
      parser.add_argument("--limit", default=cpu, type=int, dest='limit', help="Limit the number of cores to occupy when utilizing the pprocess option -p", metavar="limit")
      def positiveint(x):
        x = int(x)
        if x <= 0:
          raise ValueError, x + " is a non-positive integer"
        return x
      parser.add_argument("--addrandom", default = 0, type = positiveint, help = "Add N random sequences in each generation. Newly added sequences can immediately be used for mutations and recombinations. Sequences are created according to the same scheme as the initial population", metavar = "N")
      parser.add_argument("--randominit", action = "store_const", const = RandomInitialise, dest = "initialise", help = "Create initial population by drawing random sequences compatible with all target structures.")
      parser.add_argument("--designinit", action = "store_const", const = DesignInitialise, dest = "initialise", help = "Create initial population by repeated applications of RNAinverse. With more than one target structure, a random one is chosen for each initialisation. This is the default behaviour.")
      parser.add_argument("--fileinit", action = "store", type = FileInitialise, dest = "initialise", help = "Read initial sequences from file F - any failed read is replaced with a randomly initialised sequence, cf. --randominit.", metavar = "F")
      parser.add_argument("-i", dest = "track", default = False, action = "store_true", help = "Assign unique identifier to all founder sequences, keep track of which founders have contributed to each individual, and output this information with the individuals.")
      parser.add_argument("-v", default = 0, const = 2, type = int, nargs = '?', help = "Set verbosity to level N", metavar = "N")
      parser.add_argument("-l", default = None, type = FileType('w'), help = "Write log information to file F", metavar = "F")
      parser.add_argument("-b", default = True, action = "store_true", help = "Break off search when a perfect solution (as defined by the objective) has been found. This is the default behaviour.")
      parser.add_argument("-c", action = "store_false", dest = "b", help = "Continue search, even when a perfect solution (as defined by the objective) has been found.")
      parser.add_argument("-p", default=False, dest="usepprocess", action="store_true", help="Use usepprocess module to allow parallelization of Genetic algorithm across multiple cores")
      parser.add_argument("--predsoft", default = 'RNAfold', type = str, dest='softtype', action="store", nargs='?', help=" Choose folding software from RNAfold, PPfold, Both (combination of RNAfold and PPfold), and pknotsRG. pknotsRG must be selected for pseudoknotted sequences.", metavar = "")
      parser.add_argument("--motif", dest = "motif", default = False, help = "Enter motif as string of dots for unspecified indices and nucleotides, e.g. ..AC..G.")
      group = parser.add_argument_group(title = "Nucleotide sampling distributions", description = "Set distributions from which bases are drawn when mutating (and possibly initialising, see option --randominit) sequences. For nucleotides the order is A, C, G, U, for base pairs it is AU, CG, GC, UA, GU, UG.")
      group.add_argument("--unpaired", nargs = 4, type = float, default = None, help = "Set distribution from which to draw bases in unpaired positions to p1,...,p4", metavar = "p")
      group.add_argument("--basepair", nargs = 6, type = float, default = None, help = "Set distribution from which to draw base pairs to p1,...,p6", metavar = "p")
      group.add_argument("--paired", nargs = 4, type = float, default = None, help = "Set distribution from which to draw bases in paired positions to p1,...,p4", metavar = "p")
      group = parser.add_argument_group(title = "Dependency induced mutations", description = "Set strategy for propagating effects of a mutation across the dependency structure it is part of")
      group.add_argument("--uninformeddependent", dest = "wdependent", const = MutateDependent, action = "store_const", help = "When mutating positions because of dependency on position of original mutation, simply draw from the paired sampling distribution")
      group.add_argument("--weightdependent", dest = "wdependent", const = WeightedMutateDependent, action = "store_const", help = "When mutating positions because of dependency on position of original mutation, bias draw by the scores used to choose mutation position (see mutation site selection schemes)")
      group.add_argument("--retaindependent", dest = "wdependent", const = RetainMutateDependent, action = "store_const", help = "When mutating positions because of dependency on position of original mutation, retain the current nucleotide if possible")
      group = parser.add_argument_group(title = "Mutation sequence selection schemes")
      group.add_argument("--mutateall", dest = "mutate", default = Mutate, const = Mutate, action = "store_const", help = "Create equally many mutants from each sequence in current population")
      group.add_argument("--mutaterandom", dest = "mutate", const = RandomMutate, action = "store_const", help = "Choose sequences for mutation uniformly at random with replacement")
      group.add_argument("--mutatefit", dest = "mutate", const = FitnessWeightedMutate, action = "store_const", help = "Choose sequences for mutation at random based on fit")
      group = parser.add_argument_group(title = "Mutation site selection schemes", description = "In addition to choosing a single scheme, schemes can be combined using the --combinedposition option.")
      group.add_argument("--randomposition", action = WeightedAction("position", noargs(RandomPosition)), nargs = "?", help = "Use uniformly random selection of position to mutate", metavar = "")
      group.add_argument("--wrongposition", action = WeightedAction("position", noargs(WrongPosition)), nargs = "?", help = "Choose a random position with wrong predicted structure to mutate", metavar = "")
      group.add_argument("--boltzmannposition", action = WeightedAction("position", noargs(ProbabilityPosition)), nargs = "?", help = "Choose position to mutate based on the Boltzmann probability that it has an incorrect structure", metavar = "")
      group.add_argument("--minboltzmannposition", action = WeightedAction("position", noargs(MinProbabilityPosition)), nargs = "?", help = "Choose position to mutate based on maximum Boltzmann probability over all targets that it has an incorrect structure", metavar = "")
      group.add_argument("--logboltzmannposition", action = WeightedAction("position", noargs(LogProbabilityPosition)), nargs = "?", help = "Choose position to mutate based on negative logarithm of it having correct structure", metavar = "")
      group.add_argument("--thresholdboltzmannposition", action = WeightedAction("position", numarg(ThresholdProbabilityPosition)), nargs = "+", help = "Choose position to mutate based on the fraction of target structures for which the Boltzmann probability is not at least T. Takes a single argument.", metavar = "T")
      group.add_argument("--relativeboltzmannposition", action = WeightedAction("position", noargs(RelativeProbabilityPosition)), nargs = "?", help = "Choose position to mutate based on difference between it having correct structure and the most probable incorrect structure.", metavar = "")
      group.add_argument("--combinedposition", nargs = 0, action = CombineAction("combinedposition", CombinedPosition, "position"), help = "Set position selection scheme to be a weighted combination of all ensuing positional selection schemes. Any positional selection scheme previously specified is ignored, and all ensuing positional selection schemes need to have an added weight argument as first argument.")
      group = parser.add_argument_group(title = "Cross over points selection schemes", description = "In addition to choosing a single scheme, schemes can be combined using the --combinedxover option.")
      group.add_argument("--randomxover", action = WeightedAction("xover", (RandomCut, SelectCut)), nargs = "?", help = "Choose cross over points uniformly at random", metavar = "")
      group.add_argument("--correctxover", action = WeightedAction("xover", (CorrectCut, WeightedSelectCut)), nargs = "?", help = "Choose cross over points according to the fraction of positions with correct structure in the relevant regions of the recombining sequences", metavar = "")
      group.add_argument("--boltzmannxover", action = WeightedAction("xover", (ProbabilityCut, WeightedSelectCut)), nargs = "?", help = "Choose cross over points according to the fraction of expected number of positions with correct structure under the Boltzmann distribution, taken over the relevant regions of the recombining sequences", metavar = "")
      group.add_argument("--minboltzmannxover", action = WeightedAction("xover", (MinProbabilityCut, WeightedSelectCut)), nargs = "?", help = "Choose cross over points according to sum of minimum probability of having correct structure under the Boltzmann distribution, taken over all targets", metavar = "")
      group.add_argument("--prodboltzmannxover", action = WeightedAction("xover", (ProductProbabilityCut, WeightedSelectCut)), nargs = "?", help = "Choose cross over points according to sum of products of probability of having correct structure under the Boltzmann distribution, taken over all targets", metavar = "")
      group.add_argument("--thresholdboltzmannxover", action = WeightedAction("xover", (numarg(ThresholdProbabilityCut), WeightedSelectCut)), nargs = "+", help = "Choose cross over points according to the fraction of positions with probability of correct structure exceeding threshold T in the relevant regions of the recombining sequences. Takes a single argument.", metavar = "T")
      group.add_argument("--relativeboltzmannxover", action = WeightedAction("xover", (RelativeProbabilityCut, WeightedSelectCut)), nargs = "?", help = "Choose cross over points according to sum of scores based on difference between Boltzmann probability of correct structure and maximum probability of incorrect structure, taken over the relevant regions of the recombining sequences", metavar = "")
      group.add_argument("--combinedxover", nargs = 0, action = CombineAction("combinedxover", (CombinedCut, WeightedSelectCut), "xover"), help = "Set cross over point selection scheme to be a weighted combination of all ensuing cross over point selection schemes. Any cross over point selection scheme previously specified is ignored, and all ensuing cross over point selection schemes need to have an added weight argument as first argument.")
      
      group = parser.add_argument_group("Schemes for selecting recombining pairs")
      group.add_argument("--randompairs", dest="pairs", const=SelectPair, action="store_const", help="Select random pairs of individuals for cross overs")
      group.add_argument("--weightedpairs", dest="pairs", const=WeightedSelectPair, action="store_const", help="Select pairs by independently choosing individuals relative to their fit")
      group.add_argument("--combinedpairs", dest="pairs", const=CombinedSelectPair, action="store_const", help="Select pairs based on the average weight of a cross over point squared for the pair") 
      group = parser.add_argument_group("Objectives", description="Objectives define the sorting order of solutions and stopping criteria for the genetic algorithm. In addition to choosing a single objective, objectives can be combined using either --disjunctiveobjective or --conjunctiveobjective options.")
      group.add_argument("--wrongobjective", action=WeightedAction("objective", WrongObjective), nargs="?", help="Use number of positions with incorrect predicted structure as objective, with a sequence with predicted structure(s) perfectly matching the target structure(s) considered a perfect solution", metavar="")
      group.add_argument("--boltzmannobjective", action=WeightedAction("objective", numarg(ProbabilityObjective, opt=0)), nargs="*", help="Use expected number of incorrect positions under the Boltzmann distribution as objective. If a precision p is specified, the requirement for a perfect solution is that all structural features (base pairs and unpaired positions) have probability at most p of being wrong, with default behaviour corresponding to p = 0. Takes one optional argument.", metavar="p")
      group.add_argument("--minboltzmannobjective", action=WeightedAction("objective", numarg(MinProbabilityObjective, opt=0)), nargs="*", help="Use sum over all positions of minimum Boltzmann probability of correct structure taken over all target structures as objective. If a precision p is specified, the requirement for a perfect solution is for each position the target structure that has the least Boltzmann probability of being correct is incorrect with probability at most p, with default behaviour corresponding to p = 0. Takes one optional argument.", metavar="p")
      group.add_argument("--prodboltzmannobjective", action=WeightedAction("objective", numarg(ProductProbabilityObjective, opt=0)), nargs="*", help="Use sum over positions of 1 minus product over target structures of having correct structure, i.e. the probability of not having all target structures correct under assumption of independence, as objective. If a precision p is specified, the requirement for a perfect solution is for probability of not having all target structures correct at any position does not exceed p, with default behaviour corresponding to p = 0. Takes one optional argument.", metavar="p")
      group.add_argument("--thresholdboltzmannobjective", action=WeightedAction("objective", numarg(ThresholdProbabilityObjective)), nargs="+", help="Use number of positions where Boltzmann probability of correct structure does not exceed threshold T as objective. A perfect solution is one where every position has a Boltzmann probability of at least T of the correct structure for all target structures. Takes a single argument.", metavar="T")
      group.add_argument("--relativeboltzmannobjective", action=WeightedAction("objective", numarg(RelativeProbabilityObjective, rang=(-1, 1), opt=0)), nargs="*", help="Use score based on difference between Boltzmann probability of correct structure and maximum probability of incorrect structure as objective. If a precission p is specified, the requirement for a perfect solution is that this difference is at least p for all positions, with default behaviour corresponding to p = 0. Takes one optional argument.", metavar="p")
      group.add_argument("--ensembleobjective", action=WeightedAction("objective", numnumarg(EnsembleProbabilityObjective, rang1=(0, None), rang2=(0, None))), nargs="+", help="Use sum of log probabilities of target structures in Boltzmann ensemble and variation of log probabilities weighted by w as objective. A perfect solution is a solution where this objective does not exceed threshold T. Takes mandatory arguments w and T.", metavar="w T")
      group.add_argument("--disjunctiveobjective", action=CombineAction("combinedobjective", DisjunctiveObjective, "objective"), nargs=0, help="Set objective to be a weighted combination of all ensuing objective schemes, with a perfect solution being a solution that is perfect under at least one of the schemes. Any objective previously specified is ignored, and all ensuing objective schemes need to have an added weight argument as first argument.")
      group.add_argument("--conjunctiveobjective", action=CombineAction("combinedobjective", ConjunctiveObjective, "objective"), nargs=0, help="Set objective to be a weighted combination of all ensuing objective schemes, with a perfect solution being a solution that is perfect under all of the schemes. Any objective previously specified is ignored, and all ensuing objective schemes need to have an added weight argument as first argument.")
      group = parser.add_argument_group("Fitness functions", description="The fit of a solution determines whether it is carried over to the next generation in the genetic algorithm. Often the objective and the fit will be very similar. However, there may be cases where are more detailed or alternative fit is specified to better direct the search or explore the search space. For example, the objective may be to find a sequence with a predicted optimal structure perfectly matching a specified target, while a fit function using the Boltzmann distribution's information about probabilities of target structure elements may better identify the most ideal positions for mutations and recombinations.")
      group.add_argument("--wrongfitness", action=WeightedAction("fit", WrongFitness), nargs="?", help="Use fraction of positions with incorrect predicted structure as fit", metavar="")
      group.add_argument("--boltzmannfitness", action=WeightedAction("fit", ProbabilityFitness), nargs="?", help="Use fraction of expected number of positions with incorrect structure under the Boltzmann distribution as fit", metavar="")
      group.add_argument("--minboltzmannfitness", action=WeightedAction("fit", MinProbabilityFitness), nargs="?", help="Use the length normalised sum over all positions of the maximum Boltzmann probability of incorrect structure taken over all target structures as fit", metavar="")
      group.add_argument("--prodboltzmannfitness", action=WeightedAction("fit", ProductProbabilityFitness), nargs="?", help="Use length normalised sum over all positions of 1 minus the product of the Boltzmann probabilities of having correct structure, taken over all target structures", metavar="")
      group.add_argument("--thresholdboltzmannfitness", action=WeightedAction("fit", numarg(ThresholdProbabilityFitness)), nargs="+", help="Use fraction of positions where Boltzmann probability of correct structure does not exceed threshold T as fit. Takes a single argument.", metavar="T")
      group.add_argument("--relativeboltzmannfitness", action=WeightedAction("fit", RelativeProbabilityFitness), nargs="?", help="Use differences at each position between probability of correct structure and maximum probability of incorrect structure under the Boltzmann distribution as fit", metavar="")
      group.add_argument("--ensemblefitness", action=WeightedAction("fit", numarg(EnsembleProbabilityFitness, rang=(0, None))), nargs="+", help="Use sum of log probabilities of target structures in Boltzmann ensemble and variation of log probabilities weighted by w as fit. This number is divided by target length to obtain a `per site contribution', to make it more comparable to the rest of the fit possibilities that are all defined in a per site manner. Takes a single argument.", metavar="w")
      group.add_argument("--combinedfitness", nargs=0, action=CombineAction("combinedfitness", CombinedFitness, "fit"), help="Set fit function to be a weighted combination of all ensuing fit schemes. Any fit scheme previously specified is ignored, and all ensuing fit schemes need to have an added weight argument as first argument.")
      def nonnegnumber(x):
        try:
          x = float(x)
          if x == int(x):
            x = int(x)
        except ValueError:
          raise ValueError, "Expected a non-negative number, got {0}".format(x)
        return x
      group.add_argument("--diversity", type = nonnegnumber, nargs = "?", const = 1.0, help = "Augment fit with average fraction of positions identical to sequences already selected when choosing individuals progressing to the next generation, weighted by W (default is W = 1). This augmentation is done independent of whether the --combinedfitness option has been specified.", metavar = "W")
      return parser
    
def createMonster(options):
  class Monster(options.initialise, options.wdependent):
      pass
  # Add classes from aspects that can have a variable number of classes
  if "combinedposition" in options:
    # Make Monster a subclass of all specified mutation position
    # selection schemes.
    for c in options.position:
      class Monster(c[1], Monster):
        pass
    class Monster(options.combinedposition(options.position), Monster):
      pass
  else:
    # Only a single mutation position selection scheme specified
    class Monster(options.position, Monster):
      pass
  if "combinedxover" in options:
    # Make Monster a subclass of all specified mutation position
    # selection schemes.
    schemes = []
    for c in options.xover:
      class Monster(c[1][0], Monster):
        pass
      schemes.append((c[0], c[1][0]))
    class Monster(options.combinedxover[0](schemes), Monster):
      pass
    selectcut = options.combinedxover[1]
  else:
    # Only a single cross over selection scheme specified
    class Monster(options.xover[0], Monster):
      pass
    selectcut = options.xover[1]
  if "combinedobjective" in options:
    # Make Monster a subclass of all specified objectives
    for c in options.objective:
      class Monster(c[1], Monster):
        pass
    class Monster(options.combinedobjective(options.objective), Monster):
      pass
  else:
    # Only a single objective specified
    class Monster(options.objective, Monster):
      pass
  if "combinedfitness" in options:
    # Make Monster a subclass of all specified fit functions
    for c in options.fit:
      class Monster(c[1], Monster):
        pass
    class Monster(options.combinedfitness(options.fit), Monster):
      pass
  else:
    class Monster(options.fit, Monster):
      pass
  
  class Monster(Individual,Monster): 
      pass
  return Monster    

# Platform independent laugh player
def laugh():
    if platform == "win32":
      import winsound
      winsound.PlaySound("laugh_sound.wav", winsound.SND_FILENAME)
    elif platform == "darwin":
      import subprocess
      audio_file = "laugh_sound.wav"
      return_code = subprocess.call(["afplay", audio_file])


    
if __name__=='__main__':

  tstart = time()

  # Parse command line arguments
  parser=createParser()
  try:
        options = parser.parse_args()
  except ValueError, m:
    parser.print_help(stderr)
    print >> stderr, "\n", m
    exit(2)
  set_defaults(options)
  
  # Parse cutting method
  if "combinedxover" in options:
    selectcut = options.combinedxover[1]
  else:
    selectcut = options.xover[1]
  # Select standard diversity reduction scheme
  if "diversity" in options and options.diversity != 0:
    reduceclass = DiversityReduce(options.diversity)
  # Pareto method for reduction selected
  elif "diversity" in options and options.diversity == 0:
    reduceclass = ParetoReduce
  else:
    reduceclass = Reduce
 
  
  class Monster(createMonster(options)):
      pass

  # Only a single cross over selection scheme specified
  # Well, two classes actually... Should inherit classes controling how populations behave
  class Frankenstein(Population, Recombine, options.mutate, selectcut, options.pairs, reduceclass):
    pass

  basedist = BaseDistribution(options.unpaired, options.basepair, options.paired)
  if options.b:
    # Set up function for determining when a perfect design has been reached
    def perfection(members):
      for i in members:
        if i.perfectsolution():
          return True
      return False
  else:
    perfection = None

  

  # Set up function for printing progress
  if options.v <= 0:
    report = None
  elif options.v == 1:
    def report(pop, g, final = False):
      stdout.write('#')
  else:
    def report(pop, g):
      print "Finished generation", g + 1
      m = sorted([(i.getobjective(), i.getfitness(), i) for i in pop.members])[0]
      print "Best design is", int2str(m[2].seq)
      print "Objective: {0}; Fitness: {1}".format(m[0], m[1])
      for t in pop.gettarget():
        print combinestructures(m[1].getstruc(t.temperature).bracket(), t.bracket()) + ", T =", t.temperature

    
  # Standardise base distribution selection
  if basedist != None and not isinstance(basedist, BaseDistribution):
    w = []
    for name in ["unpaired", "basepair", "paired"]:
      if name in basedist:
        w.append(basedist[name])
      elif hasattr(name, basedist):
        w.append(getattr(name, basedist))
      else:
        w.append(None)
    basedist = BaseDistribution(*w)
    
  # Set up population
  pop = Frankenstein(options.target, Monster, options.softtype, basedist, plimit=options.limit, motif = options.motif)
  
  
  t1=time()
  # Run genetic algorithm over population
  finpop = ga(options.s, options.g, options.m, options.x, pop, logfile=options.l, stoppingcriteria=perfection, nrandom=options.addrandom, report=report, track=options.track, usepprocess=options.usepprocess)
  t2=time()
  print 'total time in ga: ',t2-t1

  # Finished running so output results
  laugh()
  output_population(finpop)
  tfin = time()
