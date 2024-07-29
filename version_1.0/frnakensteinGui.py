from Tkinter import *
import frnakenstein
import os, multiprocessing
from time import time
from cStringIO import StringIO
from individual import BaseDistribution, Individual
from individual import Individual, RandomInitialise, DesignInitialise, FileInitialise, RandomPosition, WrongPosition, ProbabilityPosition, MinProbabilityPosition, LogProbabilityPosition, ThresholdProbabilityPosition, RelativeProbabilityPosition, CombinedPosition, RandomCut, CorrectCut, ProbabilityCut, MinProbabilityCut, ProductProbabilityCut, ThresholdProbabilityCut, RelativeProbabilityCut, CombinedCut, MutateDependent, WeightedMutateDependent, RetainMutateDependent, WrongObjective, ProbabilityObjective, MinProbabilityObjective, ProductProbabilityObjective, ThresholdProbabilityObjective, RelativeProbabilityObjective, EnsembleProbabilityObjective, ConjunctiveObjective, DisjunctiveObjective, WrongFitness, ProbabilityFitness, MinProbabilityFitness, ProductProbabilityFitness, ThresholdProbabilityFitness, RelativeProbabilityFitness, EnsembleProbabilityFitness, CombinedFitness
from population import Population, Mutate, RandomMutate, FitnessWeightedMutate, SelectPair, WeightedSelectPair, CombinedSelectPair, SelectCut, WeightedSelectCut, Reduce, DiversityReduce, ParetoReduce, Recombine
from structure import Structure, Target, int2str, combinestructures
import ProgressBar
from idlelib.WidgetRedirector import WidgetRedirector
import TabBarClass
import tkFont
import tkFileDialog 
import re
import gc


#Color Constants used throughout the setup.  Can be changed to create desired effect
color1='#C71713'
color2='black'
green2='#1D5717'
green3='#6EC664'
red1='#981C1F'
rootcolor='black'
maincolor=color1
maintabcolor='black'
inputandoutputcolor=green3

#used this to stop the program from insided the method FrankGui.ga
global stopfrank
stopfrank=False
global allfranks
allfranks=set()


    
#These next classes are customized visual elements of the GUI
class GuiFrame(Frame):
    def __init__(self,master=None,color='black'):
        Frame.__init__(self,master=master)
        self['bd']=1
        self['relief']=SUNKEN
        self['padx']=3
        self['pady']=3
        self['bg']=color
class FrankButton(Button):
    def __init__(self,frame,text,command):
        Button.__init__(self,frame,text=text,command=command)
        self.config(highlightbackground=color1)
        self.config(bd=1)
    
class FrankLabel(Label):
    def __init__(self,frame,text='',bd=1,width=None,relief=NORMAL,anchor=W,color=None,textcolor='black'):
        self.text=text
        if width==None:
            Label.__init__(self,frame,text=text,bd=bd,anchor=anchor)
        else:
            Label.__init__(self,frame,text=text,bd=bd,anchor=anchor)
            
        bgcolor=frame['bg']
        if bgcolor=='black':
            self.config(fg='white')
        if color==None:
            self.config(bg=frame['bg'])
        else:
            self.config(bg=color)
        
        self.config(bd=0)
        
class FrankEntry(Entry):
    def __init__(self,master=None,width=1,color=green3,state=NORMAL):
        Entry.__init__(self,master=master,width=width,state=state)
        if color!=None:
            self.config(bg=color)
        else:
            self.config(bg=master['bg'])
        self.config(highlightbackground=color2)
       
        
class FrankMenubutton(Menubutton):
    def __init__(self,master=None,width=None,height=None,color=green3):
        Menubutton.__init__(self,master=master)
        if height!=None:
            self.config(height=height) 
        if width!=None:
            self.config(width=width)
        self.config(bg=color)
        


class FrankOptionMenu(OptionMenu):
    def __init__(self,master, choice, command=None, color=color2,*choices):
        
        OptionMenu.__init__(self,master,choice,*choices,command=command)
        self.config(bg=color)
        self.config(highlightbackground=color1)
class FrankRadiobutton(Radiobutton):
    def __init__(self, master, variable, value, text, width=None, color=green3):
        
        
        Radiobutton.__init__(self,master,var=variable,value=value,text=text)
        if width!=None:
            self.config(width=width)
        self.config(bg=color)

class FrankCheckbutton(Checkbutton):
    def __init__(self, master, variable, onvalue=True, offvalue=False, color=None):
        
        
        Checkbutton.__init__(self,master,var=variable,onvalue=onvalue,offvalue=offvalue)
        if color!=None:
            self.config(bg=color)
        
        


class ReadOnlyText(Text):
     def __init__(self, *args, **kwargs):
        Text.__init__(self, *args, **kwargs)
        self.redirector = WidgetRedirector(self)
        self.insert = \
            self.redirector.register("insert", lambda *args, **kw: "break")
        self.delete = \
            self.redirector.register("delete", lambda *args, **kw: "break")
#imports the progress bar class.    
class FrankProgress(ProgressBar.Meter):
    def __init__(self,*argvs):
        ProgressBar.Meter.__init__(self,*argvs)
        self['bd']=1
        self['relief']=SUNKEN
        self['padx']=3
        self['pady']=3


#This is the main gui class that laysout all of the other visual components within frames and 
#handles user input with the appropriate responses.  I've tried to select names for the various components
#based on what their roles are. 
class FrankGui():
     
     

    #Different options for fitness, objective, mutation position selection, and cross over selection
    fitargkeys=["--wrongfitness","--boltzmannfitness","--minboltzmannfitness","--prodboltzmannfitness","--thresholdboltzmannfitness","--relativeboltzmannfitness","--ensemblefitness",'--combinedfitness']
    objargkeys=["--wrongobjective","--boltzmannobjective","--minboltzmannobjective","--prodboltzmannobjective","--thresholdboltzmannobjective","--relativeboltzmannobjective","--ensembleobjective","--disjunctiveobjective","--conjunctiveobjective"]
    posargkeys=["--randomposition","--wrongposition","--boltzmannposition","--minboltzmannposition","--logboltzmannposition","--thresholdboltzmannposition","--relativeboltzmannposition","--combinedposition"]
    xoverargkeys=["--randomxover","--correctxover","--boltzmannxover","--minboltzmannxover","--logboltzmannxover","--thresholdboltzmannxover","--relativeboltzmannoxver","--combinedxover"]
    #This Method goes through and sets all of the default values for each of the
    #input widgets present in the interface
   
    def setDefaults(self):
        
        
        self.RanNumInt.set(0)
        self.addranEntry.delete(0, END)
        self.addranEntry.insert(0, 0)
        self.PopSizeInt.set(50)
        self.addPopSizeEntry.delete(0,END)
        self.addPopSizeEntry.insert(0,50)
        self.addGenEntry.delete(0,END)
        self.addGenEntry.insert(0,'-')
        self.initialStr.set('design')
        self.initial2.select()
        self.mutateVal.set('fitness')
        self.mutate3.select()
        self.combineVal.set('combined')
        self.combine3.select()
        self.mutchoice.set('random')
        self.recchoice.set('combined')
        self.softchoice.set('RNAfold')
        self.objchoice.set('wrong')
        self.fitchoice.set('combined')
        self.softchoice.set('RNAfold')
        
        #We want the wrong objective 
        self.CombinedEntries[0].config(state=NORMAL)
        self.CombinedEntries[0].delete(0,END)
        self.CombinedEntries[0].insert(0,1)
        self.CombinedValues[0].set(1)
        
        for i in xrange(1,len(self.CombinedEntries)):
            self.CombinedEntries[i].config(state=NORMAL)
            self.CombinedEntries[i].delete(0,END)
            self.CombinedEntries[i].insert(0,0)
            self.CombinedValues[i].set(0)
            
        #These end values are for the combined options and should always be zero!
        self.CombinedValues[-1].set(0)    
        self.CombinedValues[-2].set(0)
         
        
        self.CombinedFitValues[0].set(1)
        self.CombinedFitEntries[0].config(state=NORMAL) 
        self.CombinedFitEntries[0].delete(0,END)
        self.CombinedFitEntries[0].insert(0,1)
        self.CombinedFitValues[1].set(1)
        self.CombinedFitEntries[1].config(state=NORMAL) 
        self.CombinedFitEntries[1].delete(0,END)
        self.CombinedFitEntries[1].insert(0,1)
        self.CombinedFitValues[2].set(0)
        self.CombinedFitEntries[2].config(state=NORMAL) 
        self.CombinedFitEntries[2].delete(0,END)
        self.CombinedFitEntries[2].insert(0,0)
        self.CombinedFitValues[3].set(2)
        self.CombinedFitEntries[3].config(state=NORMAL) 
        self.CombinedFitEntries[3].delete(0,END)
        self.CombinedFitEntries[3].insert(0,2)
        for i in xrange(4,len(self.CombinedFitEntries)):
            self.CombinedFitEntries[i].config(state=NORMAL) 
            self.CombinedFitEntries[i].delete(0,END)
            self.CombinedFitEntries[i].insert(0,0)
            self.CombinedFitValues[i].set(0)
        self.CombinedFitValues[-1].set(0)
        
        #parallel processing won't work on windows platform.
        system=sys.platform
        if system=='win32':
            self.useparallellimitentry.config(state=DISABLED)
            self.useparallelbox.config(state=DISABLED)
        else:
            self.useparallelbox.select()
            self.useparallellimitentry.delete(0,END)
            self.useparallelLimit.set(multiprocessing.cpu_count())
            self.useparallellimitentry.insert(0,self.useparallelLimit.get())
        
        self.keepsearching.set(False)
        self.containpseudo.set(False)
        
        self.CombinedRecValues[0].set(0)
        self.CombinedRecEntries[0].config(state=NORMAL) 
        self.CombinedRecEntries[0].delete(0,END)
        self.CombinedRecEntries[0].insert(0,0)
        self.CombinedRecValues[1].set(1)
        self.CombinedRecEntries[1].config(state=NORMAL) 
        self.CombinedRecEntries[1].delete(0,END)
        self.CombinedRecEntries[1].insert(0,1)
        self.CombinedRecValues[2].set(1)
        self.CombinedRecEntries[2].config(state=NORMAL) 
        self.CombinedRecEntries[2].delete(0,END)
        self.CombinedRecEntries[2].insert(0,1)
        for i in xrange(2,len(self.CombinedRecEntries)):
            self.CombinedRecEntries[i].config(state=NORMAL) 
            self.CombinedRecEntries[i].delete(0,END)
            self.CombinedRecEntries[i].insert(0,0)
            self.CombinedRecValues[i].set(0)
        self.CombinedRecValues[-1].set(0)
        
        self.CombinedMutValues[0].set(0)
        self.CombinedMutEntries[0].config(state=NORMAL) 
        self.CombinedMutEntries[0].delete(0,END)
        self.CombinedMutEntries[0].insert(0,0)
        self.CombinedMutValues[1].set(2)
        self.CombinedMutEntries[1].config(state=NORMAL)
        self.CombinedMutEntries[1].delete(0,END)
        self.CombinedMutEntries[1].insert(0,2)
        self.CombinedMutValues[2].set(1)
        self.CombinedMutEntries[2].config(state=NORMAL)
        self.CombinedMutEntries[2].delete(0,END)
        self.CombinedMutEntries[2].insert(0,1)
        for i in xrange(2,len(self.CombinedRecEntries)):
            self.CombinedMutEntries[i].config(state=NORMAL)
            self.CombinedMutEntries[i].delete(0,END)
            self.CombinedMutEntries[i].insert(0,0)
            self.CombinedMutValues[i].set(0)
        self.CombinedMutValues[-1].set(0)
        
        
        
       
    #If some one wants to specify a particular nucleotide distribution
    #for creating a sequence or mutating it
    def getDistChoice(self,choice):
            if choice=='unpaired':
              distchoices=['A','C','G','U']
              for i in xrange(6):
                    if i>3:
                        self.bpentrys[i].grid_remove()
                        self.bpdistlabels[i].grid_remove()
                    else:
                        self.bpdistlabels[i]['text']=distchoices[i]+':'
                        self.bpdistlabels[i].grid(row=1,column=2*i)
                        self.bpentrys[i].grid(row=1,column=2*i+1)
              self.distributionsframe.grid(row=2,padx=3,pady=3)
                
            elif choice=='paired':
                
                distchoices=['A','C','G','U']
                for i in xrange(6):
                    if i>3:
                        self.bpentrys[i].grid_remove()
                        self.bpdistlabels[i].grid_remove()
                    else:
                        self.bpdistlabels[i]['text']=distchoices[i]+':'
              
                        self.bpdistlabels[i].grid(row=0,column=2*i)
                        self.bpentrys[i].grid(row=0,column=2*i+1)
                self.distributionsframe.grid(row=2,padx=3,pady=3)
                
            elif choice=='basepair':
                distchoices=['AU','CG','GC','UA','GU','UG']
                for i in xrange(3):
                    self.bpentrys[i].grid(row=0,column=2*i+1)
                    self.bpdistlabels[i]['text']=distchoices[i]+':'
                    self.bpdistlabels[i].grid(row=0,column=2*i)
                for i in xrange(3):
                    self.bpdistlabels[i+3]['text']=distchoices[i+3]+':'
                    self.bpdistlabels[i+3].grid(row=1,column=2*i)
                    self.bpentrys[i+3].grid(row=1,column=2*i+1)
                self.distributionsframe.grid(row=2,padx=3,pady=3)
    
    #These "get" functions just take the input from the boxes.  most are activating by
    #hitting the "RunGA" button            
    def getSoftware(self,choice):
        self.softchoice.set(choice)
    
        
    def getObjective(self,choice):
        ind=self.objchoices.index(choice)
        if choice!='conjunctive' and choice !='disjunctive':
            for i in xrange(len(self.CombinedEntries)):
                self.CombinedEntries[i].config(state=NORMAL)
                if i==ind:
                    self.CombinedValues[i].set(1)
                    self.CombinedEntries[i].delete(0,END)
                    self.CombinedEntries[i].insert(0,1)
                    self.CombinedEntries[i].config(state='readonly')
                    
                else:
                    self.CombinedValues[i].set(0)
                    self.CombinedEntries[i].delete(0,END)
                    self.CombinedEntries[i].insert(0,0)
                    self.CombinedEntries[i].config(state=DISABLED)
        else:
             for i in xrange(len(self.CombinedEntries)):
                self.CombinedEntries[i].config(state=NORMAL)
                
    def getFitness(self,choice):
        ind=self.fitchoices.index(choice)
        if choice!='combined':
            for i in xrange(len(self.CombinedFitEntries)):
                self.CombinedFitEntries[i].config(state=NORMAL)
                if i==ind:
                    self.CombinedFitValues[i].set(1)
                    self.CombinedFitEntries[i].delete(0,END)
                    self.CombinedFitEntries[i].insert(0,1)
                    self.CombinedFitEntries[i].config(state='readonly')
                    
                else:
                    self.CombinedFitValues[i].set(0)
                    self.CombinedFitEntries[i].delete(0,END)
                    self.CombinedFitEntries[i].insert(0,0)
                    self.CombinedFitEntries[i].config(state=DISABLED)
        else:
             for i in xrange(len(self.CombinedFitEntries)):
                self.CombinedFitEntries[i].config(state=NORMAL)
             
    def getMutationSite(self,choice):
        ind=self.mutchoices.index(choice)
        if choice!='combined':
            for i in xrange(len(self.CombinedFitEntries)):
                self.CombinedFitEntries[i].config(state=NORMAL)
                if i==ind:
                    self.CombinedMutValues[i].set(1)
                    self.CombinedMutEntries[i].delete(0,END)
                    self.CombinedMutEntries[i].insert(0,1)
                    self.CombinedMutEntries[i].config(state='readonly')
                    
                else:
                    self.CombinedMutValues[i].set(0)
                    self.CombinedMutEntries[i].delete(0,END)
                    self.CombinedMutEntries[i].insert(0,0)
                    self.CombinedMutEntries[i].config(state=DISABLED)
        else:
             for i in xrange(len(self.CombinedFitEntries)):
                self.CombinedMutEntries[i].config(state=NORMAL)
                
    def getRecombinationSite(self,choice):
        ind=self.recchoices.index(choice)
        if choice!='combined':
            for i in xrange(len(self.CombinedRecEntries)):
                self.CombinedRecEntries[i].config(state=NORMAL)
                if i==ind:
                    self.CombinedRecValues[i].set(1)
                    self.CombinedRecEntries[i].delete(0,END)
                    self.CombinedRecEntries[i].insert(0,1)
                    self.CombinedRecEntries[i].config(state='readonly')
                    
                else:
                    self.CombinedRecValues[i].set(0)
                    self.CombinedRecEntries[i].delete(0,END)
                    self.CombinedRecEntries[i].insert(0,0)
                    self.CombinedRecEntries[i].config(state=DISABLED)
        else:
             for i in xrange(len(self.CombinedFitEntries)):
                self.CombinedFitEntries[i].config(state=NORMAL)
    
    #Time  statements are threaded through the GA and it can calculate the average
    #time remaining based on the number of generations left to go through. 
    def timeupdate(self,time):
        
        time=round(time,2)
        self.timeentry.config(state=NORMAL)
        self.timeentry.delete(0,END)
        self.timeentry.insert(0, str(time))
        self.timeentry.config(state=DISABLED)
    
    #writes full output to selected file
    def writeout(self):
        f=open(self.outputfile.get(),'a')
        try:
            f.write(self.fileout.getvalue())
        except IOError:
            print IOError.strerror,'Cannot write out to file'
    
    #resultprint puts text content into the main text box, FrankGui.resultsbox
    #selecting clear=True erases everything first.    
    def resultprint(self,str,clear=False):
        
        """if self.writetofile.get():
            f=open(self.outputfile.get(),'a')
            try:
                f.write(str)
            except IOError:
                print IOError.strerror,'Cannot write out to file' """
        self.resultsbox.config(state=NORMAL)
        if clear:
            self.resultsbox.delete(1.0, END)
        self.resultsbox.insert(END, str)
        self.resultsbox.insert(END,'\n')
        self.resultsbox.config(state=NORMAL)
    
    #Had to take the ga method out of frnakenstein and define it as an bound method in the GUI class
    #so that we could update the progress components while the GA was running.  
    def ga (self,npop, ngen, nmut, nxov, pop, logfile=None, stoppingcriteria=None, nrandom=None, report=None, track=False, usepprocess=False):
  # Set up distribution for drawing nucleotide
      
      
      if track:
        stop = pop.addrandom(npop, idfunc = (lambda x: str(x)), stoppingcriteria = stoppingcriteria)
      else:
        stop = pop.addrandom(npop, stoppingcriteria = stoppingcriteria)
      avgtimepre=0
      avgptime=0
      
      t1=time()
      if usepprocess:
          new=pop.massprecompute(pop.new)
          #avgptime+=ptime
          if new!=None:
              pop.new=new
      else:
          pop.precompute(pop.new)
      t2=time
      #avgtimepre+=(t2-t1)
    
      pop.addnew()
      
      
      
      if stop:
        # Stopping criteria met during initialisation
        if self.writetofile.get():
          frnakenstein.output_population(pop, logfile)
          print >> logfile, "Stopping criteria met after adding {0} initial sequences to population".format(len(pop.members))
          self.finpop=pop
        return
      # Evolve population for ngen generations
      avgtimec=0
      avgtimem=0
      avgtimered=0
      avgtimepre=0
      
      #These are all just variable to keep track of the times
      #each of the parts of the GA take in order to update the progress bar
      #and estimate the total amount of time left in the algorithm
      onegentime=1
      onegenmuttime=1
      onegenprec1time=1
      onegenprec2time=1
      onegenrectime=1
      onegenredtime=1
      
      for g in xrange(ngen):
           
        self.root.update()
        if self.StopFrank.get():
            self.resultprint("******User ended GA on Generation: "+str(g)+' of '+str(ngen)+'*******',clear=True)
            break
        percentdone=0
        # Create nmut mutated sequences
        strt=time()
        if usepprocess:
            # Create nmut mutated sequences
            t1=time()
            pop.pmutate(nmut)
            t2=time()
            onegenmuttime=(g*onegenmuttime+(t2-t1))/(g+1)
            percentdone+=onegenmuttime/(onegentime*ngen)
            self.frankprogress.set(float(g)/ngen+percentdone)
            #print 'time mutating: ',t2-t1
            avgtimem+=(t2-t1)
            
            # Create nxov sequences generated by crossover
            t1=time()
            pop.precombine(nxov)
            t2=time()
            onegenrectime=(g*onegenrectime+(t2-t1))/(g+1)
            percentdone+=onegenrectime/(onegentime*ngen)
            self.frankprogress.set(float(g)/ngen+percentdone)
            
            avgtimec+=(t2-t1)
            #print 'time recombining: ',t2-t1
        else:
            t1=time()
            pop.mutate(nmut)
            t2=time()
            onegenmuttime=(g*onegenmuttime+(t2-t1))/(g+1)
            percentdone+=onegenmuttime/(onegentime*ngen)
            self.frankprogress.set(float(g)/ngen+percentdone)
            avgtimem+=(t2-t1)
            #print 'time mutating: ',t2-t1
            t1=time()
            pop.recombine(nxov)
            t2=time()
            onegenrectime=(g*onegenrectime+(t2-t1))/(g+1)
            percentdone+=onegenrectime/(onegentime*ngen)
            self.frankprogress.set(float(g)/ngen+percentdone)
            avgtimec+=(t2-t1)
            #print 'time recombining: ',t2-t1
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
        onegenprec1time=(g*onegenprec1time+(t2-t1))/(g+1)
        percentdone+=onegenprec1time/(onegentime*ngen)
        self.frankprogress.set(float(g)/ngen+percentdone)
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
        onegenprec2time=(g*onegenprec2time+(t2-t1))/(g+1)
        percentdone+=onegenprec2time/(onegentime*ngen)
        self.frankprogress.set(float(g)/ngen+percentdone)
        avgtimepre+=(t2-t1)
        pop.addnew()
        
        
        # Reduce population to fixed size
        t1=time()
        pop.reduce(npop)
        t2=time()
        onegenredtime=(g*onegenredtime+(t2-t1))/(g+1)
        percentdone+=onegenredtime/(onegentime*ngen)
        self.frankprogress.set(float(g)/ngen+percentdone)
        avgtimered+=(t2-t1)
        
        
        if self.writetofile.get():
          # Output current status to log file
          print >> logfile, "\n\n=========== Population after generation {0} ==========".format(g)
          frnakenstein.output_population(pop, logfile)
    
        if report != None:
          report(pop, g)
    
        if stoppingcriteria != None and stoppingcriteria(pop.members.__iter__()):
          # Stop genetic algorithm
          if self.writetofile.get():
            print >> logfile, "Stopping criteria met after generation", g
          break
        
        # Promote random sequences to full members
        pop.addnew()
        stpt=time()
        
        onegentime=((g)*onegentime+(stpt-strt))/(g+1)
        timeleft=onegentime*(ngen-g+1)
        self.timeupdate(timeleft)
        
        self.frankprogress.set(float(g+1)/ngen)
    
      #These print out to the logfile
      if self.writetofile.get():
          print>>logfile, 'total time mutating: ', avgtimem
          print>>logfile,'average time mutating: ',avgtimem/ngen
          print>>logfile,'total time recombining: ', avgtimec
          print>>logfile, 'average time recombining: ',avgtimec/ngen
          print>>logfile, 'total time precomputing:',avgtimepre
      #These print out to the result box
      print 'total time reducing:',avgtimered 
      print 'total time mutating: ', avgtimem
      print 'average time mutating: ',avgtimem/ngen
      print 'total time recombining: ', avgtimec
      print 'average time recombining: ',avgtimec/ngen
      print 'total time precomputing:',avgtimepre
      #print 'total time in precompute task:',avgptime
      print 'total time reducing:',avgtimered 
      
      #Instead of having the ga return anything as it does in the command line version
      #since it is a bound method, we can just have it update the instance variable
      self.finpop=pop 
      
    #Emergency stop!
    def stopGA(self):
        self.StopFrank.set(True)
        self.frankrun.config(text="Run GA",command=self.runGA)
        
    #Calls self.runGA and resets the button to show stop  
    def hitrunGA(self):
        self.frankrun.config(text="Stop GA",command=self.stopGA)
        self._job=self.root.after(0, self.runGA)
    
    
    '''This is the central method of FrankGui class.  It is activated by pressing the RunGa
     button.  It goes through and gets the input from all of the components of the gui, and sets
     the corresponding variables accordingly.  These are all used to generate a list of options
     that get passed to the frnakenstein option parser.  It also parses the structure input out of the 
     main input box and sets through accordingly.  Finally it creates the Monster and Frnakenstein class,
     creates an instance of the Frnakenstein class (a population) and runs that population through the 
     ga method and outputs the results.
      '''     
    def runGA(self):
        self.frankrun.config(text="Stop GA",command=self.stopGA)
        self.StopFrank.set(False)
        tstrt=time()
        parser=frnakenstein.createParser()
        inputstring=self.strucinput.get(1.0, END).encode('ASCII')
        motifinput=self.motifinput.get(1.0,END).encode('ASCII')
        inputargs=[]
        
        self.PopSizeInt.set(self.addPopSizeEntry.get())
        try:
            size=self.PopSizeInt.get()
            inputargs.append('-s')
            inputargs.append(str(size))
            str(self.PopSizeInt.get()
                      )
        except ValueError:
            self.resultprint('Population size not valid')
            return
        
        try:
            self.RanNumInt.set(self.addranEntry.get())
            ran=self.RanNumInt.get()
            if ran!=0:
                inputargs.append('--addrandom')
                inputargs.append(str(ran))
        except ValueError:
            self.resultprint('Add Random number not valid')
            return
        
        if self.initialStr.get()=='random':
            inputargs.append('--randominit')
        else:
            inputargs.append('--designinit')
        
        inputargs.append(self.mutateVal.get())
        inputargs.append(self.combineVal.get())
        
        if self.containpseudo.get():
            inputargs.append('--predsoft')
            inputargs.append('pknotsRG')
        
        #These next four block all control input for the options that have the posibility of being combined
        #on the second tab:  position,xover,fitness, and objective.    
        try:
            if self.fitchoice.get()=='combined':
                if self.containpseudo.get():
                    raise ZeroDivisionError
                inputargs.append('--combinedfitness')
                for i,entry in enumerate(self.CombinedFitEntries):
                    
                    self.CombinedFitValues[entry].set(self.CombinedFitEntries[entry].get())
                    try:
                        weight=self.CombinedFitValues[i].get()
                    except ValueError:
                        self.resultprint(str(ValueError.message)+'Fit value not accepted.  Using weight of 0', clear=True)
                        continue
                    if weight==0:
                        continue
                    else:
                        inputargs.append(self.fitargkeys[i])
                        inputargs.append(str(weight))
            else:
                i=self.fitchoices.index(self.fitchoice.get())
                inputargs.append(self.fitargkeys[i])
        except ZeroDivisionError:
            self.resultprint('You cannot use these obj options with pseudo knots.  Using Wrong Objective instead')
            
        
        try:
            if self.objchoice.get()=='conjunctive' or self.objchoice.get()=='disjuntive':
                if self.containpseudo.get():
                    raise ZeroDivisionError
    
                inputargs.append(self.objargkeys[self.objchoices.index(self.objchoice.get())])
                
                for i,entry in enumerate(self.CombinedEntries):
                    try:
                        self.CombinedValues[entry].set(self.CombinedEntries[entry].get())
                        weight=self.CombinedValues[i].get()
                    except ValueError:
                        self.resultprint(str(ValueError.message)+'Objective value not accepted.  Using weight of 0', clear=True)
                        continue
                    if weight==0:
                        continue
                    else:
                        inputargs.append(self.objargkeys[i])
                        inputargs.append(str(weight))
            else:
                i=self.objchoices.index(self.objchoice.get())
                inputargs.append(self.objargkeys[i])
        except ZeroDivisionError:
            self.resultprint('You cannot use these Objective options with pseudo knots. Using Wrong Fitness instead')
        
        try:
            if self.recchoice.get()=='combined':
                if self.containpseudo.get():
                    raise ZeroDivisionError
    
                inputargs.append(self.xoverargkeys[self.recchoices.index(self.recchoice.get())])
                
                for i,entry in enumerate(self.CombinedRecEntries):
                    self.CombinedRecValues[entry].set(self.CombinedRecEntries[entry].get())
                    try:
                        weight=float(self.CombinedRecValues[i].get())
                    except ValueError:
                        self.resultprint(str(ValueError.message)+'Combined xover value not accepted.  Using weight of 0', clear=True)
                        continue
                    if weight==0:
                        continue
                    else:
                        inputargs.append(self.xoverargkeys[i])
                        inputargs.append(str(weight))
            else:
                i=self.recchoices.index(self.recchoice.get())
                inputargs.append(self.xoverargkeys[i])
        except ZeroDivisionError:
            self.resultprint('You cannot use these Recombination options with pseudo knots.  Using wrong selection instead')
            
            
        try:
            if self.mutchoice.get()=='combined':
                if self.containpseudo.get():
                    raise ZeroDivisionError
    
                inputargs.append(self.posargkeys[self.mutchoices.index(self.mutchoice.get())])
                for i,entry in enumerate(self.CombinedMutEntries):
                    self.CombinedMutValues[entry].set(self.CombinedMutEntries[entry].get())
                    try:
                        weight=float(self.CombinedMutValues[i].get())
                    except ValueError:
                        self.resultprint(str(ValueError.message)+'Combined mutation site selection value not accepted.  Using weight of 0', clear=True)
                        continue
                    if weight==0:
                        continue
                    else:
                        inputargs.append(self.posargkeys[i])
                        inputargs.append(str(weight))
            else:
                i=self.mutchoices.index(self.mutchoice.get())
                inputargs.append(self.posargkeys[i])
        except ZeroDivisionError:
            self.resultprint('You cannot use these mutation site options with pseudo knots.  Using wrong instead')   
        
        
        
        
        
        
        if self.useparallel.get():
            inputargs.append('-p')
            try:
                limit=self.useparallelLimit.get()
                inputargs.append('--limit')
                inputargs.append(str(limit))
            except ValueError:
                self.resultprint(str(ValueError)+"The processor limit is not valid")
        
        try:
            if self.softchoice.get()=='PPfold':
                if self.containpseudo.get():
                    raise ZeroDivisionError
                inputargs.append('--predsoft')
                inputargs.append('PPfold')
        except ZeroDivisionError:
            self.resultprint('You cannot use this software option with pseudo knots.\n Using pknotsRG instead')
            
        try:
            if self.softchoice.get()=='Both':
                if self.containpseudo.get():
                    raise ZeroDivisionError
                inputargs.append('--predsoft')
                inputargs.append('Both')
        except ZeroDivisionError:
            self.resultprint('You cannot use this software option with pseudo knots.\n Using pknotsRG instead')
            
        if self.softchoice.get()=='PknotsRG':
	  inputargs.append('--predsoft')
	  inputargs.append('pknotsRG')
        
        if self.keepsearching.get():
            inputargs.append('-c')
        if inputstring=='\n':
            self.resultprint('Please enter in a target')
            return
            
    
        
        
            
        #instead of getting the options from the command line, we simple parse them from
        #the list that we just created.  This means we can use all of the machninery from
        #the original frnakenstein
        print inputargs
        options=parser.parse_args(inputargs)
        
        #Check for the presence of a motif specification and then append all of the lines together
        inputlines=inputstring.split('\n')
        motiflines=motifinput.split('\n')
        
        #Motif can also be input in multiple lines just like the structure.
        if motiflines==[]:
            motif=False
        else:
            motif=""
            for line in motiflines:
                motif=motif+line
            motif=motif.strip('@')
            
        options.motif=motif  
         
        #This is the parsing of the input from the Input Text box
        #We want each of the structures to be seperated by a 
        #black line and the last one to either have nothing after it
        #or to have an @ signalling there are not more structures left
        #each structure can also have temperature after it so check to for temps
        ind=0
        t = []
        temp = None
        cont=True
        structleft=True
        wrkingstruct=''
        for line in inputlines:
            lnsplit=line.split()
            if line=='':
                continue
            if lnsplit[-1][-1].isdigit():
                wrkingstruct=wrkingstruct+(lnsplit[0].strip())
                try :
                    temp=float(lnsplit[-1])
                except ValueError:
                    self.resultprint(str(ValueError)+"Invalid Temperature")
                try:
                    frnakenstein.add_structure(options,wrkingstruct,temp)
                    wrkingstruct=''
                except ValueError:
                    self.resultprint(str(ValueError))
                    
            else:
                wrkingstruct=wrkingstruct+(lnsplit[0].strip())
                
        if temp==None:
            try:
                frnakenstein.add_structure(options,wrkingstruct,temp)
                wrkingstruct=''
            except ValueError:
                self.resultprint(str(ValueError))
    
         
        dfG=self.addGenEntry.get()
        if dfG==('-'):
            self.GenInt.set(len(options.target)*options.target.size())
            self.addGenEntry.delete(0,END)
            self.addGenEntry.insert(0,self.GenInt.get())
        try:
            options.g=self.GenInt.get()
        except ValueError:
            self.resultprint('Number of Generations not valid')
            return  
        
        frnakenstein.set_defaults(options)
       
        if "combinedxover" in options:
            selectcut = options.combinedxover[1]
        else:
            selectcut = options.xover[1]
  
        if "diversity" in options and options.diversity != 0:
            reduceclass = DiversityReduce(options.diversity)
        elif "diversity" in options and options.diversity == 0:
            reduceclass = ParetoReduce
        else:
            reduceclass = Reduce
            
        class Monster(frnakenstein.createMonster(options)):
            pass
        
        #I am actually not sure how this works, but it appears to allow for the Monster
        #class to be pickalizable because it is defined directly in the main code block of 
        #our program.  This allows for parallelization. 
        import __main__
        __main__.Monster=Monster
        #Just to see what all we have run it with.
        print options
        
        #Set up our popluation class with the necessary options
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
            
        #Here we copy the stdout and redirect it to a new file handle.  Thus anythign that prints
        #out within the GA doesn't got to the stdout, it gets printed out at the end of the Ga
        #on our results box.  There are three file handles that get created here, mainly:
        #so that we can control what gets outputted where and in what order
#        
#                        1) New out : This is the output from the output of the final population
#                        2) GA out:  this becomes the new stdout.  Includes all of the prints 
#                        statements within the GA or contained files
#                        3) fileout:  this is the handle that gets written into the log file 
#                        if it is selected.  
        
        stdoutcpy=sys.stdout
        #
        newout=StringIO()
        
        self.gaout=StringIO()
        
        if self.writetofile.get():
            self.fileout=StringIO()
        else:
            self.fileout=False
        report=None
        #create a population
        pop = Frankenstein(options.target, Monster, options.softtype, basedist, plimit=options.limit,motif=options.motif)
        #assign it to the instance variable
        self.finpop=pop
        
        
        
        sys.stdout=self.gaout
        self.ga(options.s, options.g, options.m, options.x, pop, logfile=self.fileout, stoppingcriteria=perfection, nrandom=options.addrandom, report=report, track=options.track, usepprocess=options.usepprocess)
        if self.StopFrank.get():
            self.resultprint("Final Population:")
            frnakenstein.output_population(self.finpop, file=newout)
            self.resultprint(newout.getvalue())
            return
        frnakenstein.output_population(self.finpop, file=newout)
        if self.writetofile.get():
            frnakenstein.output_population(self.finpop, file=self.fileout)
        self.resultprint(newout.getvalue(),clear=True)
        self.resultprint(self.gaout.getvalue())
        
        
        #Restore the stdout
        sys.stdout=stdoutcpy
        
        #write to file
        if self.writetofile.get():
            self.writeout()
        
        #update progress bar to finished and time to 0
        self.frankprogress.set(1)
        self.timeupdate(0)
        tfin=time()
        #Print the totaltime and reset the run GA button
        self.resultprint("Total Time: "+str(tfin-tstrt))
        self.frankrun.config(text="RunGA",command=self.hitrunGA )
        
    
    #method to create an open file dialog box and to get the file name from the user
        
    def OpenFile(self):
        self.writetofile.set(True)
        self.outputfile.set(tkFileDialog.asksaveasfilename(initialdir=os.getcwd()))
        self.outputfileEntry.config(state=NORMAL)
        self.outputfileEntry.delete(0,END)
        self.outputfileEntry.insert(0, os.path.basename(self.outputfile.get()))
        self.outputfileEntry.config(state='readonly')
     
    # This is where most of the instance variables for FrankGUi are setup (some are user
    # dependent and don't get created until the GA is run.  Most of the names are pretty   
    #self explanitory (but then again all programmers think that).
    def __init__(self):
        global allfranks
        allfranks.add(self)
        self.root = Tk()
        self.FrankTabBar=TabBarClass.TabBar(self.root,'Main Menu')
        self.FrankTabBar.config(bg=maintabcolor)
        
        
        self.root.title("Frnakenstein")
        self.root.geometry("875x675")
        self.root['padx']=1
        self.root['pady']=1
        self.root['bg']=rootcolor
        
        #On the first tab bar we want all of the input, outputs, and commands to run Frnakenstein
        self.maintab=TabBarClass.Tab(self.root,'Main Menu',color=color2)
        self.maintab.config(bg=maincolor)
        self.inputframe=GuiFrame(self.maintab,color='black')
        self.inputlabel=FrankLabel(self.inputframe,'Input the target stucture and press RunGA')
        self.inputlabel.config(bd=1,relief=RAISED,bg=color1)
        self.inputlabel.grid(stick=W)
        self.strucinput=Text(self.inputframe, width=50, height=10 )
        self.strucinput.config(highlightbackground='black')
        self.strucinput.config(bg=inputandoutputcolor)
        self.strucinput.grid(stick=E)
        self.inputscroll=Scrollbar(self.inputframe)
        self.inputscroll.config(troughcolor=color1)
        self.inputscroll.config(elementborder=1)
        self.inputscroll.config(highlightcolor=color1)
        self.inputscroll.config(relief=RAISED)
        self.strucinput.config(yscrollcommand=self.inputscroll.set)
        self.inputscroll.config(command=self.strucinput.yview)
        
        self.inputscroll.grid(row=1,column=1,sticky=N+S)
        
        self.targets=[]
        
        self.resultlabel=FrankLabel(self.inputframe, text="Results",bd=1, relief=RAISED,color=color1).grid(sticky=W)
        self.resultsbox=ReadOnlyText(self.inputframe, width =50,height=10,bg='#EDEDED')
        self.resultscroll=Scrollbar(self.inputframe)
        
        self.resultsbox.config(yscrollcommand=self.resultscroll.set)
        self.resultsbox.config(highlightbackground='black')
        self.resultsbox.config(bg=inputandoutputcolor)
        self.resultscroll.config(command=self.resultsbox.yview)
        
        
        self.resultsbox.grid(sticky=E)
        self.resultscroll.grid(row=3,column=1,stick=N+S)
        self.frankprogress=FrankProgress(self.inputframe)
        self.frankprogress.grid()
        self.inputframe.grid(row=0,column=0,columnspan=1,padx=10,pady=10,rowspan=6)
       
        self.FrankTabBar.add(self.maintab)
        
        #BEGINNING OF THE SECOND TAB
        #Objective and Fitness and Mutation and Recombination objects took up too much room so
        #moved them into another bar tab to access
        
        self.objfitmutrectab=TabBarClass.Tab(self.root,'Other Parameters',color=color2)
        self.objfitmutrectab.config(bg=color1)
        
        self.ObjFitFrame=GuiFrame(self.objfitmutrectab)
        self.ObjFitLabel=FrankLabel(self.ObjFitFrame,text="Objective and Fitness Options").grid(row=0,column=0,sticky=W)
        self.ObjFitFrame.grid(row=0,column=0,sticky=W,padx=3,pady=3,rowspan=4,columnspan=4)
        
        self.ObjFrame=GuiFrame(self.ObjFitFrame)
        self.objLabel=FrankLabel(self.ObjFrame,text='Objective:')
        self.objLabel.grid(sticky=W)
        self.objmenubut=FrankMenubutton(self.ObjFrame, width=20, height=2)
        self.objmenubut.pack_propagate(0)
        self.objchoices=['wrong','boltzmann','min boltzmann','product boltzmann','threshold boltzmann','relative boltzmann','ensemble','disjunctive','conjunctive']
        self.objchoice=StringVar()
        self.objmenubut.menu=FrankOptionMenu(self.objmenubut,self.objchoice,self.getObjective,'grey',*self.objchoices)
        self.objmenubut.menu.pack(expand=True,fill='both')
        self.objmenubut['menu']=self.objmenubut.menu
        self.objmenubut.grid(row=0,column=1,sticky=W)
        
        self.ObjFrame.grid(row=1,column=0,rowspan=2,padx=3,pady=3,sticky=NW)
        
        
        self.combineFrame=GuiFrame(self.ObjFrame,color=color1)
        self.CombinedValues=len(self.objchoices)*[DoubleVar()]
        self.CombinedEntries={}
        for i,choice in enumerate(self.objchoices):
            if i==7 or i==8:
                break
            label=FrankLabel(self.combineFrame,text=choice)
            label.grid(row=i%3, column=i/3*2)
            label.config(width=8)
            label.config(anchor=W)
            self.CombinedEntries[i]=FrankEntry(self.combineFrame, width=2)
            self.CombinedEntries[i].grid(row=i%3, column=i/3*2+1)
        self.combineFrame.grid(columnspan=3,padx=3,pady=3)
    
        self.FitFrame=GuiFrame(self.ObjFitFrame)
        self.FitLabel=FrankLabel(self.FitFrame,text='Fitness:')
        self.FitLabel.grid(sticky=W)
        self.fitmenubut=FrankMenubutton(self.FitFrame, width=20, height=2)
        self.fitmenubut.pack_propagate(0)
        self.fitchoices=['wrong','boltzmann','min boltzmann','product boltzmann','threshold boltzmann','relative boltzmann','ensemble','combined']
        self.fitchoice=StringVar()
        self.fitmenubut.menu=FrankOptionMenu(self.fitmenubut,self.fitchoice,self.getFitness,'grey',*self.fitchoices)
        self.fitmenubut.menu.pack(expand=True,fill='both')
        self.fitmenubut['menu']=self.fitmenubut.menu
        self.fitmenubut.grid(row=0,column=1,sticky=W)
        
        self.combineFitFrame=GuiFrame(self.FitFrame,color=color1)
        self.CombinedFitValues=len(self.fitchoices)*[DoubleVar()]
        self.CombinedFitEntries={}
        for i,choice in enumerate(self.fitchoices):
            if i==7:
                break
            label=FrankLabel(self.combineFitFrame,text=choice)
            label.grid(row=i%3, column=i/3*2)
            label.config(width=8)
            label.config(anchor=W)
            self.CombinedFitEntries[i]=FrankEntry(self.combineFitFrame, width=2)
            self.CombinedFitEntries[i].grid(row=i%3, column=i/3*2+1)
        self.combineFitFrame.grid(columnspan=3,padx=3,pady=3)
        self.FitFrame.grid(row=4,column=0,rowspan=2,columnspan=2,padx=3,pady=3)
        
        self.MutRecSiteLabel=FrankLabel(self.ObjFitFrame,text="Site Selection Options").grid(row=0,column=2)
        self.MutSiteFrame=GuiFrame(self.ObjFitFrame)
        self.MutSiteLabel=FrankLabel(self.MutSiteFrame,text='Mutation').grid(sticky=W)
        self.mutmenubut=FrankMenubutton(self.MutSiteFrame,width=15,height=2)
        self.mutmenubut.pack_propagate(0)
        self.mutchoices=['random','wrong','boltzmann','min boltzmann','log boltzmann','threshold boltzmann','relative boltzmann','combined']
        self.mutchoice=StringVar()
        self.mutmenubut.menu=FrankOptionMenu(self.mutmenubut,self.mutchoice,self.getMutationSite,'grey',*self.mutchoices)
        self.mutmenubut.menu.pack(expand=True,fill='both')
        self.mutmenubut['menu']=self.mutmenubut.menu
        self.mutmenubut.grid(row=0,column=1,sticky=W)
        self.MutSiteFrame.grid(row=1,column=2,columnspan=2,rowspan=2,padx=3,pady=3,sticky=N)
        
        self.combineMutFrame=GuiFrame(self.MutSiteFrame,color=color1)
        self.CombinedMutValues=len(self.mutchoices)*[DoubleVar()]
        self.CombinedMutEntries={}
        for i,choice in enumerate(self.mutchoices):
            if i==7:
                break
            label=FrankLabel(self.combineMutFrame,text=choice)
            label.grid(row=i%3, column=i/3*2)
            label.config(width=8)
            label.config(anchor=W)
            self.CombinedMutEntries[i]=FrankEntry(self.combineMutFrame, width=2)
            self.CombinedMutEntries[i].grid(row=i%3, column=i/3*2+1)
        self.combineMutFrame.grid(columnspan=3,padx=3,pady=3)
        
        
        
        
        
        self.RecSiteFrame=GuiFrame(self.ObjFitFrame)
        self.RecSiteLabel=FrankLabel(self.RecSiteFrame,text='Recombination').grid(sticky=W)
        self.recmenubut=FrankMenubutton(self.RecSiteFrame,width=15,height=2)
        self.recmenubut.config(highlightbackground='black')
        self.recmenubut.pack_propagate(0)
        self.recchoices=['random','wrong','boltzmann','min boltzmann','log boltzmann','threshold boltzmann','relative boltzmann','combined']
        self.recchoice=StringVar()
        self.recmenubut.menu=FrankOptionMenu(self.recmenubut,self.recchoice,self.getRecombinationSite,'grey',*self.recchoices)
        self.recmenubut.menu.pack(expand=True,fill='both')
        self.recmenubut['menu']=self.recmenubut.menu
        self.recmenubut.grid(row=0,column=1,sticky=W)
        self.RecSiteFrame.grid(row=3,column=2,rowspan=2,columnspan=2,padx=3,pady=3,sticky=N)
        
        self.combineRecFrame=GuiFrame(self.RecSiteFrame,color=color1)
        self.CombinedRecValues=len(self.recchoices)*[DoubleVar()]
        self.CombinedRecEntries={}
        for i,choice in enumerate(self.recchoices):
            if i==7:
                break
            label=FrankLabel(self.combineRecFrame,text=choice)
            label.grid(row=i%3, column=i/3*2)
            label.config(width=8)
            label.config(anchor=W)
            self.CombinedRecEntries[i]=FrankEntry(self.combineRecFrame, width=2)
            self.CombinedRecEntries[i].grid(row=i%3, column=i/3*2+1)
        self.combineRecFrame.grid(columnspan=3,padx=3,pady=3)
        
        #Base Pair Distribution Frames are in the other parameter tab as well
        self.BpDistributionframe=GuiFrame(self.objfitmutrectab)
        self.topBplabel=FrankLabel(self.BpDistributionframe,text='Select Nucleotide\ndistribution').grid(sticky=W)
        self.bpmenubut=FrankMenubutton(self.BpDistributionframe,width=14, height=2)
        self.bpmenubut.pack_propagate(0)
        self.distchoices=['unpaired','paired','basepair']
        self.bpdistchoice=StringVar()
        self.bpmenubut.menu=FrankOptionMenu(self.bpmenubut,self.bpdistchoice,self.getDistChoice,'grey',*self.distchoices)
        self.bpmenubut.menu.pack(expand=True,fill='both')
        self.bpmenubut['menu']=self.bpmenubut.menu
        self.unpairdist=[DoubleVar for i in xrange(4)]
        self.basepairdist=[DoubleVar for i in xrange(6)]
        self.pairdist=[DoubleVar for i in xrange(4)]
        self.bpmenubut.grid(sticky=W,padx=3,pady=3)
        
        self.distributionsframe=GuiFrame(self.BpDistributionframe)
        self.bpentrys=[FrankEntry(self.distributionsframe,width=3) for i in xrange(6)]
        self.bpdistlabels=[FrankLabel(self.distributionsframe) for i in xrange(6)]
        self.BpDistributionframe.grid(row=4,column=0,rowspan=3,padx=3,pady=3,sticky=W)
        
        self.motifframe=GuiFrame(self.objfitmutrectab)
        self.motiflabel=FrankLabel(self.motifframe,text='Input Motif').grid(row=0,sticky=W)
        self.motifinput=Text(self.motifframe, width=70, height=4)
        self.motifinput.config(highlightbackground='black')
        self.motifinput.config(bg=inputandoutputcolor)
        
        self.motifinput.grid(row=1)
        self.motif=StringVar()
        self.motifframe.grid(row=4,column=1)

        
        self.motifscroll=Scrollbar(self.motifframe)
        self.motifscroll.config(troughcolor=color1)
        self.motifscroll.config(elementborder=1)
        self.motifscroll.config(highlightcolor=color1)
        self.motifscroll.config(relief=RAISED)
        self.motifinput.config(yscrollcommand=self.inputscroll.set)
        self.motifscroll.config(command=self.motifinput.yview)
        
        self.motifscroll.grid(row=1,column=1,sticky=N+S)
    
        
        
        
        
        
        self.FrankTabBar.add(self.objfitmutrectab)
        
        
        
        #THIS IS BACK IN THE FIRST TAB
        self.GAparametersFrame=GuiFrame(self.maintab)
        
        self.PopGenRanframe=GuiFrame(self.GAparametersFrame,color=color1)
        self.GALabel=FrankLabel(self.PopGenRanframe,text="Genetic Algorithm\nParameters").grid(row=0, columnspan=2)
        self.addranLabel=FrankLabel(self.PopGenRanframe,text="Add Random:").grid(row=1,column=0)
        self.addranEntry=FrankEntry(self.PopGenRanframe, width=2)
        self.addranEntry.grid(row=1,column=1)
        self.RanNumInt=IntVar()
        self.addGenLabel=FrankLabel(self.PopGenRanframe,text="# Generations:").grid(row=2,column=0)
        self.addGenEntry=FrankEntry(self.PopGenRanframe, width=2) 
        self.addGenEntry.grid(row=2,column=1)
        self.GenInt=IntVar()
        self.addPopSizeLabel=FrankLabel(self.PopGenRanframe,text="Population Size:").grid(row=3,column=0)
        self.addPopSizeEntry=FrankEntry(self.PopGenRanframe, width=2)
        self.addPopSizeEntry.grid(row=3,column=1,padx=3,pady=3)
        self.PopSizeInt=IntVar()
        self.PopGenRanframe.grid(row=0,column=0,rowspan=3,columnspan=2,padx=3,pady=3,sticky=N+S+E+W)
        
        self.InitialiseFrame=GuiFrame(self.GAparametersFrame,color=color1)
        self.initialStr=StringVar()
        self.initialLabel=FrankLabel(self.InitialiseFrame,text="Initialization\nOptions")
        self.initialLabel.config(width='15')
        self.initialLabel.grid(columnspan=2)
        self.initial1=FrankRadiobutton(self.InitialiseFrame,variable=self.initialStr,value='random', text="Random", width=10)
        self.initial1.grid(columnspan=2)
        self.initial2=FrankRadiobutton(self.InitialiseFrame,variable=self.initialStr,value='design',text="Design",width=10)
        self.initial2.grid(columnspan=2)
        self.InitialiseFrame.grid(row=4,column=0,rowspan=3,columnspan=2,padx=3,pady=3,sticky=N+S+E+W)
       
      
        
        self.MutRecSeqLabel=FrankLabel(self.GAparametersFrame,text="Sequence Selection\nOptions")
        self.MutRecSeqLabel.grid(row=0,column=2)
        
        self.MutateOptionFrame=GuiFrame(self.GAparametersFrame,color=color1)
        self.mutateVal=StringVar()
        self.mutlabel=FrankLabel(self.MutateOptionFrame,text="Mutation")
        self.mutlabel.config(width='15')
        self.mutlabel.pack(side='top')
        self.mutate1=FrankRadiobutton(self.MutateOptionFrame,variable=self.mutateVal,value='--mutateall', text="All", width=10)
        self.mutate1.pack(side='top')
        self.mutate2=FrankRadiobutton(self.MutateOptionFrame,variable=self.mutateVal,value='--mutaterandom',text="Random",width=10)
        self.mutate2.pack(side='top')
        self.mutate3=FrankRadiobutton(self.MutateOptionFrame,variable=self.mutateVal,value='--mutatefit',text="Fitness",width=10)
        self.mutate3.pack(side='top')
        self.MutateOptionFrame.grid(row=0, rowspan=3, column=2,padx=3,pady=3,sticky=N+S+E+W)
    
        
        self.CombineOptionFrame=GuiFrame(self.GAparametersFrame,color=color1)
        self.combineVal=StringVar()
        self.combinelabel=FrankLabel(self.CombineOptionFrame,text="Recombination")
        self.combinelabel.config(width='15')
        self.combinelabel.pack(side='top')
        self.combine1=FrankRadiobutton(self.CombineOptionFrame,variable=self.combineVal,value='--randompairs', text="Random", width=11)
        self.combine1.pack(side='top')
        self.combine2=FrankRadiobutton(self.CombineOptionFrame,variable=self.combineVal,value='--weightedpairs',text="Weighted",width=11)
        self.combine2.pack(side='top')
        self.combine3=FrankRadiobutton(self.CombineOptionFrame,variable=self.combineVal,value='--combinedpairs',text="Combined",width=11)
        self.combine3.pack(side='top')
        self.CombineOptionFrame.grid(row=4, rowspan=3, column=2,padx=3,pady=3,sticky=N+S+E+W)
        
        
        
        
        
        #Software choice options
        self.softframe=GuiFrame(self.GAparametersFrame,color=color1)
        
        self.softlabel=FrankLabel(self.softframe,text="Select Folding\nSoftware").grid(columnspan=2)
        
        self.softmenubut=FrankMenubutton(self.softframe, width=14, height=2)
        self.softmenubut.pack_propagate(0)
        self.softchoice=StringVar()
        self.softwares=['RNAFold','PPfold','Both','PknotsRG']
        self.softmenubut.menu=FrankOptionMenu(self.softmenubut,self.softchoice,self.getSoftware,'grey',*self.softwares)
        self.softmenubut.menu.pack(expand=True,fill='both')
        self.softmenubut['menu']=self.softmenubut.menu                                                                                
        self.softmenubut.grid(columnspan=2,sticky=W)
        self.softframe.grid(row=0,column=4,padx=3,pady=3,rowspan=7,sticky=N+S+E+W)
        
        
        self.containpseudo=BooleanVar()
        self.containpseudobox=FrankCheckbutton(self.softframe, self.containpseudo, onvalue=True, offvalue=False,color=color2)  
        self.containpseudoLabel=FrankLabel(self.softframe,text="Pseudo knots: ")
        self.containpseudoLabel.grid(row=4,sticky=W)
        self.containpseudobox.grid(row=4,column=1)
        
        self.useparallel=BooleanVar()
        self.useparallelbox=FrankCheckbutton(self.softframe, self.useparallel, onvalue=True, offvalue=False,color=color2)  
        self.useparallelLabel=FrankLabel(self.softframe,text="Parallelize: ")
        self.useparallelLimit=IntVar()
        self.parallelLimitLabel=FrankLabel(self.softframe, text='Process Limit:')
        self.useparallellimitentry=FrankEntry(self.softframe, width=2,color=green3)
        self.useparallelLabel.grid(row=5,sticky=W)
        self.parallelLimitLabel.grid(row=6,column=0,sticky=W)
        self.useparallelbox.grid(row=5, column=1,sticky=S)
        self.useparallellimitentry.grid(row=6,column=1,sticky=E)
        
        self.keepsearching=BooleanVar()
        self.keepsearchingbox=FrankCheckbutton(self.softframe, self.keepsearching, onvalue=True, offvalue=False,color=color2)  
        self.keepsearchingLabel=FrankLabel(self.softframe,text="Continue Search : ")
        self.keepsearchingLabel.grid(row=7,sticky=W)
        self.keepsearchingbox.grid(row=7,column=1)
        

        self.GAparametersFrame.grid(row=2,column=1,columnspan=3,rowspan=4,padx=3,pady=3,sticky=NW)
        
        
        #Box keeping track of time
        self.timeframe=GuiFrame(self.maintab)
        self.defaultimeframe=GuiFrame(self.timeframe,color=color1)
        self.defaultimeframe.grid(padx=3,pady=3)
        self.timelabel=FrankLabel(self.defaultimeframe,text='Time Left: ')
        self.defaultsbutton=FrankButton(self.defaultimeframe,text='Set Defaults', command=self.setDefaults)
        self.defaultsbutton.config(width=15)
        self.defaultsbutton.grid(row=1,column=0,columnspan=2)
        self.timeentry=FrankEntry(self.defaultimeframe,width=8,color=green3,state=DISABLED)
        self.timelabel.grid(row=0)
        self.timeentry.grid(row=0,column=1,sticky=W)
        
        
        self.timeframe.grid(row=0,column=1,sticky=W)
        self.frankrun=FrankButton(self.maintab,text="Run GA",command=self.hitrunGA)
        self.StopFrank=BooleanVar()
        self.StopFrank.set(False)
        self.frankrun.config(width=50)
        self.frankrun.grid(row=1,column=1,columnspan=2,padx=3,pady=3,sticky=W)
      
        
        self.outputframe1=GuiFrame(self.maintab)
        self.outputframe2=GuiFrame(self.outputframe1,color=color1)
        self.outputfilelabel=FrankLabel(self.outputframe2,text="Output file:")
        self.outputfilelabel.grid(row=1)
        self.outputfileEntry=FrankEntry(self.outputframe2,width=15,color=green3,state='readonly')
        self.outputfileEntry.grid(row=1,column=1)
        self.writetofile=BooleanVar()
        self.outputfile=StringVar()
        self.openbutton=FrankButton(self.outputframe2,text='Open File',command=self.OpenFile)
        self.openbutton.grid(row=2,column=1,sticky=E+W,padx=3,pady=3)
        self.outputframe2.grid(padx=3,pady=3)
        self.outputframe1.grid(row=0,column=2,sticky=W,padx=3,pady=3)
        
        
        
        
        self.frankpicture=PhotoImage(file='frankensteintext.gif')
        self.frankpiclabel=Label(image=self.frankpicture)
        self.frankpiclabel.config(border=0)
        self.frankpiclabel.grid(row=1,column=0,columnspan=2,sticky=W)
    
          
        self.FrankTabBar.show()      
        self.setDefaults()
        self.root.mainloop()
        
        
#This whole program is only one line?
def main():
    frankscreen=FrankGui()
    
if __name__=="__main__":
    main()
    
