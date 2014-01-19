#!/usr/bin/env python
"""
Operations used to determine the transition state of a molecule using
kinetic isotope effects

Automatically checks if a Display is available, and, if so, enables the
draw_graphs to draw data stored in each class Experiment
"""

#######################################################
# TODO: Write in an option to automatically check TS structure if manual monitor not wished for
# Do this by using: sum xyz. Compute sum for every connection. User inputs rxn coordinates in list of tuples
# sum of xyz for all rxn coordinates must be >= 0.5 * (at least one of the negative frequencies)

# Put in get_frequencies after freq optimization + write function to get it

# Put in get_structure after freq optimization

########################################################
# Useful for converting """ """ into nice looking html's and pdf's
# pydoc -w ./ (writes a foo.html for each foo.py inside the current directory)
# epydoc --pdf filename will convert filename into pdf

########################################################
# Is there a way to check and see how much the script slows down Gaussian?

# Code breaks when just_minimize=True since depends on self.variables that have not yet been definied
# Need to change definition of newline error? It's being called in a lot of extraneous circumstances...
# Need to rework draw_graphs function

# Potential problem: allows user to be stupid and change algorithm midway through code

import string
import subprocess
import numpy
#from scipy import *
from scipy.optimize import fmin
import os
import itertools

# Checks if an X11 Display is available, and if so, enables the program to draw data reported
# in results.txt

if os.getenv('DISPLAY') is not None:
    import pylab
    from numpy import array

    def draw_isoeff(isoeff_output, inputs):
        """
        Plots kinetic isotope effects against atomic distances

        Keyword Arguments:
        
        * isoeff_output -- list of lists containing calculated kinetic isotope
        effects for each iteration of fmin
        * inputs -- list of lists containing distances between atoms for each
        iteration of fmin

        Assumes len(isoeff_output) == len(inputs), len(isoeff_output[i]) = K
        for every i, and len(inputs[j]) = C for every j

        Returns None
        
        Saves N graphs (N = len(inputs[i]) to Isoeff results.pdf, each
        depicting the effects of the change of distance to the kinetic
        isotope effects.
        """
        lists_output = zip(*isoeff_output)
        lists_input = zip(*inputs)
        #print lists_output
        #print lists_input
        for x, one_molecule_freq in enumerate(lists_output):
            for i, distance in enumerate(lists_input):
                pylab.figure(x)
                pylab.title('Distances between atoms vs Kinetic isotope effect ' + str(x + 1))
                pylab.xlabel('Distances Between Atoms')
                pylab.ylabel('Kinetic Isotope Effect')
                pylab.plot(distance, one_molecule_freq, 'o', label='Distance ' + str(i + 1))
            pylab.legend(loc='best')
        pylab.ion()
        pylab.draw()
        pylab.savefig('Isoeff results.pdf')
        pylab.close()
        return None

    def draw_minimize(min_output, inputs):
        """
        Plots output of minimization function against atomic distances

        Keyword Arguments:

        min_output -- list of lists containing output from fmin

        inputs -- list of lists containing distances between atoms for each iteration
        of fmin

        Assumes len(min_output) == len(inputs), len(min_output[i]) == 1 for every i,
        and len(inputs[j]) == K for every j

        Returns None
        
        Saves N graphs (N = len(inputs[i]) to Minimization Output.pdf, each depicting the
        effects of the change of ditance to the output of fmin
        """
        lists_input = zip(*inputs)
        print len(lists_input)
        for i, distance in enumerate(lists_input):
    ##        pylab.figure(0)
            pylab.title('Distances between atoms vs Value of minimized function')
            pylab.xlabel('Distances Between Atoms')
            pylab.ylabel('Minimized Function Value')
            pylab.plot(distance, min_output, 'o', label='Distance ' + str(i + 1))
            pylab.ylim([-0.05, 0.05])
        pylab.legend(loc='best')
        pylab.ion()
        pylab.draw()
        pylab.savefig('Minimization output.pdf')
        pylab.close()
        return None

    def draw_graphs(txt_file):
        """
        Saves data from txt_file into two different pdf files

        Keyword Arguments:

        txt_file -- string indicating file to recieve data

        Assumes the file contains only lists of numbers, ordered such that isotope effects
        are first, then atomic distances, then the output of the minimization function

        Returns None

        Saves graphs under Minimization output.pdf and Isoeff Results.pdf
        For specifications on graphs, see draw_isoeff and draw_minimize
        """
        lines = get_lines(txt_file)
        isoeff_output = eval(lines[2])
        inputs = eval(lines[3])
        minimized_output = []
        for line, input in zip(eval(lines[5]), inputs):
            line = [line]
            minimized_output.append(line)
        draw_isoeff(isoeff_output, inputs)
        draw_minimize(minimized_output, inputs)
        return None



def get_lines(txt_file):
    """Takes lines from txt_file and returns a list of strings"""
    text = open(txt_file, 'r')
    lines = text.readlines()
    text.close()
    return lines

def write_lines(txt_file, lines):
    """Takes a list of strings and rewrites the lines into txt_file"""
    text = open(txt_file, 'w')
    text.writelines(lines)
    text.close()
    return None

def get_structure(lines):
    """
    Finds the structure of an optimized molecule

    Keyword arguments:
    
    lines -- list of strings from frequency optimization output

    Returns a list of strings detailing the structure of the optimized molecule
    """
    structure = []
    for i, line in enumerate(lines):
        if "Input orientation:" in line:
            index = i
            break  
    line = lines[index]
    while "Distance" not in line:
        structure.append(lines[index])
        index += 1
        line = lines[index]
    return structure

def check_file(file):
    """
    Verifies that file is an output from GaussView

    Keyword arguments:
    
    file -- string indicating file to be checked

    Returns None

    Will raise GaussInputError if file is not in the recognized format from GaussView
    """
    lines = get_lines(file)
    checkpts = ['%mem=', '%nproc=', '%chk=', '# ']
    for line, checkpt in zip(lines, checkpts):
        if checkpt not in line:
            raise GaussInputError
    if lines[4].split() != []:
        raise GaussInputError
    if lines[6].split() != []:
        raise GaussInputError
##    print "Files checked."
    return None

def check_float(item):
    """Returns item as a float. Raises RuntimeError if float(item) raises a ValueError"""
    try:
        return float(item)
    except ValueError:
        raise RuntimeError (str(item) + " needs to be a float")
    
def check_int(item):
    """Returns item as an int. Raises RuntimeError if int(item) raises a ValueError"""
    try:
        return int(item)
    except ValueError:
        raise RuntimeError (str(item) + " needs to be an integer")

def check_mem(item):
    """Returns item. Asserts that item is a string and that "mw" ends the string"""
    assert type(item) == str
    if item[-2] + item[-1] != 'mw':
        raise RuntimeError ("You must have 'mw' at the end of memory")
    return item

def check_constant_list(item):
    """
    Checks the format of item to act as a molecular specification in Gaussian Input files

    Keyword Arguments:

    item -- a list of strings
    
    Returns item

    Raises ConstantInputError if "F" doesn't end every string in item
    """
    for line in item:
        if line.split()[-1] != 'F':
            raise ConstantInputError                
    return item

def temp_to_lines (template, fillers, line_range):
    """
    Converts template and fillers into a list that could be interpreted by Gaussian

    Keyword Arguments:

    template -- string with n %s in it

    fillers -- list of n strings

    line_range -- integer specifying how many lines in the final list should have
    % in front of them

    Returns a list of strings that may be interpreted by Gaussian

    Substitites %s in template with fillers, splits the template at "|"
    into a list of strings. Adds a % to the number of lines specified in line_range
    """
    new_template = template % fillers
    lines = new_template.split('|')
    for i in range(line_range + 1):
        lines[i] = '%' + lines[i]
    return lines



class Experiment (object):
    """
    A class to setup and run an Experiment to find a transition state structure
    """
    def __init__(self, Rfile, RfixedGaussian, TSfile, TSfixedGaussian, IsoeffFile, Data,\
                 ConstantDistances=None, Increment=None, Algorithm='# B3LYP/6-31G**', \
                 Solvent=None, TSOpt='TS', MaxCycle=10000, Eigentest=False, Memory='50mw', Nproc=4):
        """
        Initializes an Experiment

        Keyword Arguments:

        Rfile -- string indicating the file from GaussView of the substrate (R)

        RfixedGaussian -- list of strings containing fixed molecular specifications for the substrate

        TSfile -- string indicating the file from GaussView of the transition state (TS)

        TSfixedGaussian -- list of strings constaining fixed molecular specifications for the transition state

        IsoeffFile -- string indicating the file containing the input to isoeff98

        Data -- list of tuples of experimental data.
        (experimental kinetic isotope effect, experimental error)

        Keyword Arguments with Default:

        ConstantDistances -- list of strings containing the lines in TSfixedGaussian which
        should not be influenced by fmin. Default None

        Increment -- list of integers where minimization should begin. The length of Increment must equal
        the number of distances to be minimized. Default None

        Algorithm -- string formatted as the beginning of the Gaussian command line, but only
        including the method used for calculation. Default "# B3LYP/6-31G**"

        Solvent -- string following Gaussian formatting specifying a solvent in the calculations.
        Default None

        TSOpt -- string following Gaussian formatting specifying the optimization algorithm to
        use for finding a transition state. Default "TS"

        Maxcycle --  integer substituted for the variable MaxCycle in Gaussian operations.
        Specifies number of acceptable iterations for optimization calculations. Default 10000

        Eigentest -- boolean determining to test curvature in TS calculations. Eigentest=False is
        recommended only for those with large computing budgets. Default False

        Memory -- string specifying the amount of memory to be devoted to optimization.
        Default "50mw"

        Nproc -- integer specifying number of processors to be used in optimization. Default 4

        Defines inputs, isoeff_outputs, minimized_outputs, structures, and frequencies as empty
        lists. These are used to document previous tests using write_results
        """
            # First check all items which will never be None
        assert all(type(foo) == str for foo in [Rfile, TSfile, IsoeffFile, Algorithm, TSOpt, Memory])
            # Check Rfile
        check_file(Rfile)
        self.Rfile = Rfile
        assert all(type(foo) == list for foo in [RfixedGaussian, TSfixedGaussian, Data])
            # Verify RfixedGaussian is in the correct format
        assert all(type(line) == str for line in RfixedGaussian)
        self.RfixedGaussian = check_constant_list(RfixedGaussian)
            # Check TSfile
        check_file(TSfile)
        self.TSfile = TSfile
            # Verify TSfixedGaussian is in the correct format
        assert all(type(line) == str for line in TSfixedGaussian)
        self.TSfixedGaussian = check_constant_list(TSfixedGaussian)
            # Currently, IsoeffFile is assumed to be in the correct format
        self.IsoeffFile = IsoeffFile
            # Check Data tuples contain floats
        assert all(type(datum) == float and type(error) == float for datum, error in Data)
        self.Data = Data
            # If ConstantDistances exists, checks types of all subsets, ensures lines are in correct format,
            # and checks that each line in ConstantDistances exists in TSfixedGaussian
        if ConstantDistances != None:
            assert type(ConstantDistances) == list
            assert all(type(line) == str for line in ConstantDistances)
            ConstantDistances = check_constant_list(ConstantDistances)
            assert all(line in TSfixedGaussian for line in ConstantDistances)
        self.ConstantDistances = ConstantDistances
            # If Increment exists, turns it into a float if it's not already and if it can't be turned, raises error
        if Increment != None:
            Increment = check_float(Increment)
        self.Increment = Increment
            # If solvent exists, ensures solvent is type string
        if Solvent != None:
            assert type(Solvent) == str
            # The following will not equal none and were already tested or the testing is obvious
        self.Solvent = Solvent
        self.Algorithm = Algorithm
        self.TSOpt = TSOpt
        self.MaxCycle = check_int(MaxCycle)
        assert type(Eigentest) == bool
        self.Eigentest = Eigentest
        self.Memory = check_mem(Memory)
        self.Nproc = check_int(Nproc)
        self.inputs = []
        self.isoeff_output = []
        self.minimized_output = []
        self.current_distances = []

    def getStatus(self):
        """Returns a list of tuples: (keyword argument, value)"""
        status = []
        for key, value in zip(self.__dict__.keys(), self.__dict__.values()):
            status.append((key, value))
        return status

    # The following are statements to change any of the optional inputs to experiments.
    # To change required files, one must begin a new experiment
    def setConstantDistances (self, ConstantDistances):
        """
        Sets self.ConstantDistances to ConstantDistances

        Keyword Arguments:

        ConstantDistances -- list of strings

        Returns None

        Asserts that ConstantDistances is a list of strings and that every line in
        ConstantDistances is in TS.fixedGaussian
        """
        assert type(ConstantDistances) == list
        assert all(type(line) == str for line in ConstantDistances)
        assert all(line in self.TSfixedGaussian for line in ConstantDistances)
        self.ConstantDistances = check_constant_list(ConstantDistances)
    def setIncrement(self, Increment):
        """
        Sets self.Increment to Increment

        Keyword Arguments:

        Increment -- float

        Returns None

        Asserts that Increment is a float, or, if it is not, modifies Increment to a float
        """
        self.Increment = check_float(Increment)
    def setAlgorithm(self, Algorithm):
        """
        Sets self.Algorithm to Algorithm

        Keyword Arguments:

        Algorithm -- string formatted as the beginning of the Gaussian command line, but only
        including the method used for calculation.

        Returns None

        Asserts that Algorithm is a string
        """
        assert type(Algorithm) == str
        self.Algorithm = Algorithm
    def setSolvent(self, Solvent):
        """
        Sets self.Solvent to Solvent

        Keyword Arguments:

        Solvent -- string following Gaussian formatting specifying a solvent in the calculation

        Returns None

        Asserts that Solvent is a string
        """
        assert type(Solvent) == str
        self.Solvent = Solvent
    def setTSOpt(self, TSOpt):
        """
        Sets self.TSOpt to TSOpt

        Keyword Arguments:

        TSOpt -- string following Gaussian formatting specifying the optimization algorithm to
        use for finding a transition state.

        Returns None

        Asserts that TSOpt is a string
        """
        assert type(TSOpt) == str
        self.TSOpt = TSOpt
    def setMaxCycle(self, MaxCycle):
        """
        Sets self.MaxCycle to MaxCycle

        Keyword Arguments:

        MaxCycle -- integer substituted for the variable MaxCycle in Gaussian operations.
        Specifies number of acceptable iterations for optimization calculations.

        Returns None

        Formats MaxCycle to an int if possible
        """
        self.MaxCycle = check_int(MaxCycle)
    def setEigentest (self, Eigentest):
        """
        Sets self.Eigentest to Eigentest

        Keyword Arguments:

        Eigentest -- boolean determining to test curvature in TS calculations. Eigentest=False is
        recommended only for those with large computing budgets.

        Returns None

        Asserts Eigentest is a boolean
        """
        assert type(Eigentest) == bool
        self.Eigentest = Eigentest
    def setMemory(self, Memory):
        """
        Sets self.Memory to Memory

        Keyword Arguments:

        Memory -- string specifying the amount of memory to be devoted to optimization

        Returns None

        Asserts Memory is in Gaussian formatting
        """
        self.Memory = check_mem(Memory)
    def setNproc(self, Nproc):
        """
        Sets self.Nproc to Nproc

        Keyword Arguments:

        Nproc -- integer specifying number of processors to be used in optimization.

        Returns None

        Formats Nproc to an int if possible
        """
        self.Nproc = check_int(Nproc)

    # Creates file to pass to Gaussian
    # Returns the lines of the file
    def gauss_view_conversion (self, kind, StopInOptFreq=False):
        """
        Creates a file titled TS.gjf or R.gjf depending on the variable kind,
        writes the necessary Gaussian lines into the file,
        and returns the lines as a list of strings. The file is formatted to be a Gaussian
        input file

        Keyword Arguments:

        kind -- either "R" or "TS", specifing which type of molecule the gaussian file represents
        
        StopInOptFreq -- boolean. If False, Gaussian will skip the extra steps required to compute
        the Raman intensities during Hartree-Fock analytic frequency calculations,
        saving 10-30% in CPU time.
        If True, the program will stop after frequency calculations

        Returns lines, the list of lines written into the new file
        """
        ##Note: I am unsure about the validity of using Freq=NoRaman, Consult with Vipender
        assert type(StopInOptFreq) == bool
        if self.Eigentest == False:
            eigentest = 'noeigentest'
        else:
            eigentest = ''
        if self.Solvent == None:
            solvent = ''
        else:
            solvent = self.Solvent
        if StopInOptFreq == False:
            template = 'mem=%s\n|Chk=./' + kind + '.chk\n|nproc=%s\n|%s Opt=(%s, %sMaxCycle=%s, AddRed) %s Freq=NoRaman test\n' 
        else:
            template = 'mem=%s\n|Chk=./' + kind + '.chk\n|nproc=%s\n|%s Opt=(%s, %sMaxCycle=%s, AddRed) %s test\n'
        if kind == 'R':
            old_lines = get_lines(self.Rfile) 
            fixed_gaussian = self.RfixedGaussian
            opt_kind = ''
        elif kind == 'TS':
            old_lines = get_lines(self.TSfile)
            fixed_gaussian = self.TSfixedGaussian
            opt_kind = self.TSOpt + ', '
        else:
            raise RuntimeError('kind specified neither R nor TS')
        fillers = (self.Memory, self.Nproc, self.Algorithm, eigentest, opt_kind, self.MaxCycle, solvent)
        lines = temp_to_lines(template, fillers, 1)
        take_rest = True
        ## This mess manipulates the lines from the GaussView files to process into Gaussian
        ## Unfortunately, this requires that GaussView output files remain the same.
        for line in old_lines[4:9]:
            lines.append(line)
        for line in old_lines[9:]:
            if line.split() != []:
                if line.split()[0] == '1':
                    take_rest = False
            if take_rest:
                lines.append(line)
        try:
            for line in fixed_gaussian:
                lines.append(line)
        except NameError:
            pass
        ## Creates new file. There is certainly a more elegant way to do this.
        new_file = kind + '.gjf'
        subprocess.call(['cp', self.TSfile, new_file])
        write_lines(new_file, lines)
        return lines

    # Runs Gaussian on either default file or lines passed in with lines argument
    # Will raise an error if Gaussian faults or fails to complete
    # Returns the output of Gaussian
    def run_gaussian(self, kind, opt=True, lines=None):
        """
        Runs gaussian on lines or on the file TS.gjf or R.gjf, specified by kind.
        Formatted to run optimization jobs and finding vibrational frequency jobs.
        Writes the output of Gaussian to TS.log or R.log.

        Keyword Arguments:

        kind -- either "R" or "TS", specifing which type of molecule the Gaussian file represents
        
        opt -- boolean specifying type of job Gaussian will run. If True, Gaussian will optimize
        the molecule. If False, Gaussian will find the vibrational frequency. Default True
        
        lines -- list of strings normally contained in .gjf file. If not None, Gaussian will
        run the job on the list of strings. Default None.

        Returns: list of strings written to "TS.log" or "R.log"; the output of the Gaussian job.
        """
        # opt=True on the first iteration of Gaussian. opt=False on the second iteration. (see powerpoint presentation)
        log_file = kind + '.log'
        if lines != None:
            StdIn = lines
        elif opt:
            StdIn = open(kind + '.gjf', 'r')
        elif not opt:
            StdIn = open('freq' + kind + '.gjf', 'r')
        try:
            subprocess.call(['g03'], stdin = StdIn, stdout=open(log_file, 'w'))
        except IOError:
            raise GnotAvailable
        log_lines = get_lines(log_file)
        if not any("GradGrad" in line for line in log_lines):
            raise GnewlineError
        if "Normal termination of Gaussian" not in log_lines[-1]:
            raise GTimeoutError
        return log_lines

    # Writes a file to be passed to Gaussian for frequency calculations
    # Returns the lines of the file
    def write_freq_file(self, kind):
        """
        Creates a file for vibrational frequency calculations, freq[kind].gjf and prepares the directory to run
        vibrational frequency calculations. Copies [kind].chk, one of the files created by Gaussian
        optimization, to freq[kind].chk, which is used by Gaussian frequency calculations.

        Keyword Arguments:

        kind -- either "R" or "TS", specifying which type of molecule the Gaussian file represents

        Returns: lines written to freq[kind].gjf
        """
        chk_file = kind + '.chk'
        subprocess.call(['cp', chk_file, 'freq' + chk_file])
        log_file = kind + '.log'
        template = 'mem=%s\n|Chk=./freq' + kind + '.chk\n|nproc=%s\n|%s Freq=NoRaman Geom=AllCheck test\n'
        fillers = (self.Memory, self.Nproc, self.Algorithm)
        lines = temp_to_lines(template, fillers, 2)
        freq_file = 'freq' + kind + '.gjf'
        subprocess.call(['cp', log_file, freq_file])
        write_lines(freq_file, lines)
        return lines

    # Changes the log file to only have G98 in it
    def change_g03 (self, kind):
        """
        Formats TS.log or R.log, specified by kind, to be able to run on Isoeff98

        Keyword Arguments:

        kind -- either "R" or "TS", specifying which type of molecule the Gaussian file represents

        Returns: None
        """
        log_file = kind + '.log'
        lines = get_lines(log_file)
        new_lines = []
        for line in lines:
            new_line = string.replace(line, 'G03', 'G98')
            new_lines.append(new_line)
        write_lines(log_file, new_lines)
        return None

    # Runs Isoeff98 and returns a list of strings of frequency outputs. Raises MissingFileError
    # if TS or R.log file is missing. Raises IsoeffError if the program faults.
    def run_isoeff(self):
        """
        Runs Isoeff98. Assumes TS.log and R.log exist and are updated. Creates [self.IsoeffFile].iso.

        No keyword arguments

        Returns: List of strings, the frequency outputs
        """
        try:
            text1 = open('TS.log', 'r')
            text2 = open('R.log', 'r')
        except IOError:
            raise MissingFileError
        subprocess.call(['./isoeff98.exe', self.IsoeffFile]) 
        txt_file = self.IsoeffFile + '.iso'
        lines = get_lines(txt_file)
        kept_lines = []
        error = True
        for line in lines:
            if 'JOB NUMBER' in line: 
                error = False
            if 'EXECUTION OF PROGRAM TERMINATED ABNORMALLY' in line:
                error = True
            if 'ISOTOPE EFFECT(' in line:
                kept_lines.append(line.split()[3])
        if error:
            raise IsoeffError
        return kept_lines

    # Rewrites TS gauss file to be able to repass into Gaussian, including new_distances into calculation
    # rewrite_gauss_view should only be run on TS
    def rewrite_gauss_view(self, new_distances, StopInOptFreq=False):
        """
        Updates TS.gif to represent a new TS molecule. This molecule will have small changes in
        the distance between molecules. These distances are in new_distances

        Keyword Arguments:

        new_distances -- a list of floats indicating the new distances between each atom.
        The new distances are determined by fmin.
        
        StopInOptFreq -- a boolean. If False, Gaussian will skip the extra steps required to compute
        the Raman intensities during Hartree-Fock analytic frequency calculations,
        saving 10-30% in CPU time.
        If True, the program will stop after frequency calculations

        Returns: List of strings. The lines just written into TS.gjf
        """
        ## Note: check to make sure setting Freq=NoRaman is still relevant
        old_lines = get_lines('TS.gjf')
        if self.Eigentest == False:
            eigentest = 'noeigentest'
        else:
            eigentest = ''
        if self.Solvent == None:
            solvent = ''
        else:
            solvent = self.Solvent
        if StopInOptFreq == False:
            template = 'mem=%s\n|Chk=./TS.chk\n|nproc=%s\n|%s Opt=(%s, %s, MaxCycle=%s) %s Freq=NoRaman Geom=(AllCheck, ModRed)\n'
        else:
            template = 'mem=%s\n|Chk=./TS.chk\n|nproc=%s\n|%s Opt=(%s, %s, MaxCycle=%s) %s Geom=(AllCheck, ModRed)\n'
        fillers = (self.Memory, self.Nproc, self.Algorithm, eigentest, self.TSOpt, self.MaxCycle, solvent)
        lines = temp_to_lines(template, fillers, 2)
        lines.append('\n')
        ## This modifies the input of new_distances and prevents excessive Gaussian calculation.
        ## It is impractical to have a distance less than 0.1 between atoms or a distance greater
        ## than 5.6. Essentially implements a floor and ceiling function
        ## Assumes this will not change the calculations by fmin
        new_distances = [max(0.1, min(5.6, dist)) for dist in new_distances]
        self.current_distances = new_distances
        print "This is new_distances seen by rewrite_gauss_view", new_distances
        if self.ConstantDistances == None:
            constant_distances = []
        else:
            constant_distances = self.ConstantDistances
        fixed_command_lines = self.TSfixedGaussian
        ## Creates the new lines to be written into TS.gjf
        x = -1
        for line in fixed_command_lines:
            if all(line.split() != constant_distance.split() for constant_distance in constant_distances):
                print "This is line, and line.split()", line, line.split()
                x += 1
                new_line = line.split()[:-2]
                new_line.append(str(new_distances[x]) + ' F\n')
                print "This is new_line", new_line
                line = string.join(new_line)
            lines.append(line)
        write_lines('TS.gjf', lines)
        return lines

    # Writes current data into inp_file
    ## Currently this only writes temporary results, ie. the results will be rewritten each time write_results
    ## is called. Currently this problem is solved by using cumulative self.isoeff_output, self.inputs, and
    ## self.minimized_output
    def write_results(self, inp_file):
        """
        Writes the list of distances optimized for, their correlating result from Isoeff98,
        and their fmin output into inp_file.

        Keyword Arguments:

        inp_file -- string, the name of the file into which the data will be written

        Returns: None
        """
        assert type(inp_file) == str
        subprocess.call(['cp', 'R.gjf', inp_file])
        lines =['Temporary results\n', 'Isoeff Outputs, inputs\n', str(self.isoeff_output) \
                + '\n', str(self.inputs) + '\n', 'Output to be minimized\n', \
                str(self.minimized_output) + '\n']
        write_lines(inp_file, lines)
        ##if os.getenv('DISPLAY') is not None:
        ##    draw_graphs(self.minimized_output, self.inputs)
        return None

    # Will return the last set of distances that were written into a list of new_distances identical to the ones already passed into the function.
    # Used in calculate_frequency and initiate files
    ## This function was written before the implementation of the Class Experiment. There is probably a more elegant way to
    ## get the new distances which does not rely on file formatting

    def get_new_distances(self, kind):
        """
        Searches self.TSfixedGaussian or self.Rfixed Gaussian to find the distances to be modified
        by the experiment.

        Keyword Arguments:

        kind -- either "R" or "TS", specifying which type of molecule the Gaussian file represents

        Returns: list of floats - the new distances
        """
        print 'Getting new distances...'
        if self.ConstantDistances == None:
            constantdistances = []
        else:
            constantdistances = self.ConstantDistances
        def drop_line(line):
            return line in constantdistances
        ## itertools.dropwhile drops any string in self.TSfixedGaussian which is in constantdistances
        ## This ensures that only distances which are modified will be returned
        if kind == 'TS':
            new_distances = itertools.dropwhile(drop_line, self.TSfixedGaussian)
        else:
            new_distances = itertools.dropwhile(drop_line, self.RfixedGaussian)
        return [float(line.split()[-2]) for line in new_distances]

#########################################################################

                   #The beginning of cumulative functions# 

#########################################################################

    
    # Makes a kind file ready to be sent into Isoeff. Takes first_iteration as a variable to indicate
    # whether to make a file or rewrite a file

    # Much of this code is creative ways to generate proper print statements

    def calculate_frequency(self, kind, new_distances=None, stop_for_optimization=False, \
                            first_iteration=False):
        """
        Prepares R or TS files to be sent into Isoeff98.
        Runs gauss_view_conversion or rewrite_gauss_view, runs run_gaussian, will allow user to manually check structure
        using stop_for_optimization. Finally, runs change_g03 to prepare [kind].log to be sent into Isoeff98

        Keyword Arguments:

        kind -- either "R" or "TS", specifying which type of molecule the Gaussian file represents
        
        new_distances -- a list of floats - the distances between the atoms, to be placed in [kind].gjf
        Default None.
        
        stop_for_optimization -- a boolean. If True, the program will halt until the user has checked the
        molecule and approved it. The user also has the option to abort the program
        if the molecule is unsatisfactory. If True, the calculation will be much slower
        since Gaussian will be called twice. Default False.
        
        first_iteration -- a boolean. If True, the program will not run self.rewrite_gauss_view, instead
        self.gauss_view_conversion. A copy will also be created of [kind].log titled
        cp_has_g03_[kind].log. Default False.

        Returns: None
        """
        if new_distances == None:
            print "new_distances was None"
            new_distances = self.get_new_distances('TS')
        assert type(kind) == str
        print_statement = "Running Gaussian for structure optimization"
        if not stop_for_optimization:
            print_statement = print_statement + " and vibrational frequency calculation..."
        ## Prepare [kind].gjf to be run in Gaussian
        if first_iteration:
            print "Writing files for Gaussian structure optimization..."
            self.gauss_view_conversion(kind, StopInOptFreq=stop_for_optimization)
            print print_statement
        else:
            print "Writing files for Gaussian structure optimization using fmin generated distances..."
            self.rewrite_gauss_view(new_distances, StopInOptFreq=stop_for_optimization)
            print print_statement
        ## If stop_for_optimization=False, running Gaussian will create [kind].log with vibrational frequencies
        ## If not, running Gaussian will create [kind].log with molecular details
        self.run_gaussian(kind, opt=True)
        # Outputs to kind.log file
        if first_iteration:
            subprocess.call(['cp', kind + '.log', 'cp_has_g03_' + kind + '.log'])
        if stop_for_optimization:
            man_check = str(raw_input("Program paused until user manually checks structure. Press enter to continue. Press q to quit"))
            if man_check != 'q':
                print "Writing files for Gaussian vibrational frequency calculations..."
                self.write_freq_file(kind)
                print "Running Gaussian for vibrational frequency calculation..."
                ## Reruns Gaussian! Will be inefficient, but will allow user to manually check structure
                self.run_gaussian(kind, opt=False)
            else:
                print "Vibrational frequency calculations aborted"
                raise Abortion       
        self.change_g03(kind)
        return None



    # Runs calculate_frequency assuming the first step has been run (that only
    # TS should be rewritten and run)
    
    # Then runs Isoeff98.exe. Adds new_distances to self.inputs, adds distance vector
    # to self.minimized_output, and adds the calculated vibrational frequencies to
    # self.isoeff_output. Writes all three new lists to 'Temporary results' using
    # self.write_results
    
    # This is the function that is minimized. Returns a distance in N-dimensional space
    # where N is the number of atoms being examined and the distance is calculated
    # between the current "position" of isotopic vibrational effects and the expected
    # "position" of isotopic vibrational effects in self.Data[0]

    ## Note: This is repeated in the function only for help( ) uses
    def run_next_steps(self, new_distances, stop_for_optimization=False, stop_for_freq=False):
        """
        Runs calculate_frequency assuming that both R.log and TS.log exist and are updated.
        Assums the first step has been run (that only TS should be rewritten and run)
    
        Then runs Isoeff98.exe. Adds new_distances to self.inputs, adds distance vector
        to self.minimized_output, and adds the calculated vibrational frequencies to
        self.isoeff_output. Writes all three new lists to 'Temporary results' using
        self.write_results
    
        This is the function that is minimized.
        
        Keyword Arguments:

        new_distances -- a list of floats - the distances between the atoms
        
        stop_for_optimization -- a boolean. If True, will allow for manual checking of the molecule
        See calculate_frequency for details. Default False
        
        stop_for_freq -- a boolean. If True, the program will pause after calculating the
        vibrational frequencies, allowing the user to view them, before calling Isoeff98.
        Default False

        Returns: a distance in N-dimensional space where N is the number of atoms being examined
        and the distance is calculated between the current "position" of isotopic vibrational
        effects and the expected "position" of isotopic vibrational effects in self.Data[0]
        """
        print "New distance(s)" , new_distances
        # Creates list of lists incrementing each factor one at a time, if more than one new_distance
        # exists. Designed to decrease number of Gaussian errors thrown

        ## Should a dict to prevent tedious recalculation be inserted?
        ##    global io_dict
        ##    rounded_new_distances = round(distance, 6) for distance in new_distances]
        ##    if rounded_new_distances in io_dict:
        ##        return io_dict(rounded_new_distances)
        print "This is self.inputs", self.inputs
        if self.inputs != []:
            iterations = [self.inputs[-1][:i] + new_distances[i:] for i in \
                      reversed(xrange(len(self.inputs[-1]))) if self.inputs[-1][i] != new_distances[i]]
        else:
            iterations = [new_distances]
        print "Redefining iterations"
        iterations = [new_distances]
        for iteration in iterations:     
            print "This is iteration in run_next_steps", iteration
            self.inputs.append(iteration)        
            try:
                self.calculate_frequency('TS', new_distances=iteration, stop_for_optimization=stop_for_optimization)
            except GTimeoutError:
                print 'GTimeoutError raised. Trying again...'
                self.calculate_frequency('TS', new_distances=iteration, stop_for_optimization=stop_for_optimization)
            if stop_for_freq:
                raw_input("Program paused for user to check frequency calculations. Press enter to continue")
            print "Running isoeff98.exe..."
            calculated_data = self.run_isoeff()
            self.isoeff_output.append(calculated_data)
            print "Isoeff output", calculated_data 
            assert len(self.Data) == len(calculated_data), \
                   "Cross-check Data and IsoeffFile. Different numbers of atoms to check were specified."
            output = sum((float(exp_datum[0]) - float(calc_datum))**2 \
                       for calc_datum, exp_datum in zip(calculated_data, self.Data))
            self.minimized_output.append(output)
            self.write_results('Temporary results')
        return output


    def initiate_files(self, stop_TS_optimization=False, stop_R_optimization=False, stop_frequency=False,\
                        increment=False, rerun=False):
        """
        Prepares the experiment to be optimized using fmin.
        Calculates the vibrational frequency for TS and R, runs isoeff, sets self.current_distances to the
        distances calculated by Gaussian. If user decides to manually increment the distances between molecules,
        initiate_files will run run_next_steps iteratively until the increment value is reached.
        
        Keyword Arguments:

        stop_TS_optimization -- a boolean. If True, will allow for manual checking of the TS molecule.
        Default False
        
        stop_R_optimization -- a boolean. If True, will allow for manual checking of the R molecule.
        Default False
        
        stop_frequency -- a boolean. If True, the program will pause after calculating the
        vibrational frequencies, allowing the user to view them, before calling Isoeff98.
        Default False.
        
        increment -- a boolean. If True, the program will run run_next_steps until self.Increment has been reached.
        Default False.
        
        rerun -- a boolean. If True, the program will assume that TS has been run before. It is to be set to True
        if TS faulted from too low a value on MaxCycles, and the TS checkpoint file is still updated.
        Assumes that R has not been run before. Default False.

        Returns: None
        """
        print "Initiating TS calculations"   
        if rerun:
            self.calculate_frequency('TS', new_distances=self.current_distances, stop_for_optimizatation=stop_TS_optimization)
        else:
            self.calculate_frequency('TS', stop_for_optimization=stop_TS_optimization, first_iteration=True)
        print 'Initiating substrate calculations'
        self.calculate_frequency('R', stop_for_optimization=stop_R_optimization, first_iteration=True)
        if stop_frequency:
            raw_input("Program paused for user to check frequency calculations. Press enter to continue")
        print "Running isoeff..."
        calculated_data = self.run_isoeff()
        print 'First data calculated', calculated_data
        new_distances = self.get_new_distances("TS")
        print new_distances
        self.current_distances = new_distances
        ## TODO: increment has not been tested, as obvious by the infinite loops. 
        ## The value 0.2 has no significance. It was a number I predicted would work
        if increment:
            def change_distance(inp, distance, pos=True):
                if pos:
                    distance = distance + 0.2
                else:
                    distance = distance - 0.2
                revised_distances = self.Increment[:inp[0] - 1]
                revised_distances.append(distance)
                revised_distances.extend(new_distances[inp[0] + 1:])
                self.run_next_steps(revised_distances)
            if self.Increment != None:
                try:
                    for input_num, distance in zip(enumerate(self.Increment), new_distances):
                        if input_num[1] > distance:
                            while input_num > distance + 0.2:
                                change_distance(input_num, distance)
                        elif input_num[1] < distance:
                            while input_num < distance - 0.2:
                                change_distance(input_num, distance, pos=False)
                except IndexError:
                    pass
        return None
    

    def run_program(self, stop_beg_TS_optimization=False, stop_beg_R_optimization=False, stop_beg_freq=False, \
                    stop_iteration_optimization=False, stop_iteration_frequency=False, increment=False, \
                    rerun=False, just_minimize=False):
        """
        Runs the optimization program. Program will run initiate_files if just_minimize is False.
        Then the program will run fmin. fmin optimizes the distance vector returned by run_next_steps.
        The initial guess passed into fmin is self.current_distances.
        run_program thens reruns self.run_isoeff to find the final vibrational frequencies. It checks to ensure
        the error between the vibrational frequencies found experimentally with the ones just returned is
        within error defined in self.Data.

        Keyword Arguments:

        stop_beg_TS_optimization -- a boolean. If True, will allow for manual checking of the TS molecule before
        optimization begins. Default False
        
        stop_beg_R_optimization -- a boolean. If True, will allow for manual checking of the R molecule before
        optimization begins. Default False
        
        stop_beg_freq -- a boolean. If True, the program will pause after calculating the first
        vibrational frequencies, allowing the user to view them, before calling Isoeff98.
        Default False.
        
        stop_iteration_optimization -- a boolean. If True, will allow for manual checking of the TS molecule at
        every step of the optimization process. Default False.
        
        stop_iteration_frequency -- a boolean. If True, the program will pause after calculating the
        vibrational frequencies at every step of the optimization process,
        allowing the user to view them, before calling Isoeff98.
        Default False.
        
        increment -- a boolean. If True, the program will run run_next_steps until self.Increment has been reached.
        This will occur before optimization begins. Default False.
        
        rerun -- a boolean. If True, the program will assume that TS has been run before. It is to be set to True
        if TS faulted from too low a value on MaxCycles, and the TS checkpoint file is still updated.
        Assumes that R has not been run before. Default False.
        
        just_minimize -- a boolean. If True, the program will assume initiate_files has already been run. This
        assumes that TS.log, TS.chk, R.log, R.chk, and self.current_distances are valid and
        updated. This option will override rerun. Default False.

        Returns: a tuple. tuple[0] is a float - the final distance returned by fmin. ie. what was minimized
        tuple[1] is a list of floats - the errors between the experimental data and the calculated
        tuple[2] is a list of strings - "Within Error" if the specific error is within acceptable error
        "NOT WITHIN ERROR" if it is outside acceptable error
        """
        if not just_minimize:
            self.initiate_files(stop_TS_optimization=stop_beg_TS_optimization, stop_R_optimization=stop_beg_R_optimization, \
                           stop_frequency=stop_beg_freq, increment=increment, rerun=rerun)
        print self.current_distances
        ## I believe the args works correctly, but warrents test cases
        final_dist = fmin(self.run_next_steps, self.current_distances, \
                          args=(stop_iteration_optimization, stop_iteration_frequency))
        calculated_data = self.run_isoeff()
        errors = [abs(float(calc_datum) - float(exp_datum[0])) for calc_datum, exp_datum in zip(calculated_data, self.Data)]
        new_errors = []
        for error, exp_datum in zip(errors, self.Data):
            if error <= exp_datum[1]:
                error = "Within Error"
            else:
                error = "NOT WITHIN ERROR. Error was: " + str(error) + " measured error was: " + str(exp_datum[1])
            new_errors.append(error)
        return final_dist, errors, new_errors

#########################################################################

                   #Program meant for interaction with users# 

#########################################################################


    def find_transition_state(self, results_label, stop_beg_TS_optimization=False, stop_beg_R_optimization=False, \
                              stop_beg_freq=False, stop_iteration_optimization=False, stop_iteration_frequency=False, \
                              increment=False, rerun=False, just_minimize=False):
        """
        Program will find a likely transition state using calculated vibrational frequencies.
        Runs the minimization program, creates or updates the file results.txt to include the latest data.
        Will draw graphs using the draw_graphs function if the environment can draw graphs.
        If program is aborted, results.txt will still be updated to include the data gained before the keyboard
        interrupt and draw_graphs will still be called if possible.

        Keyword Arguments:
        
        results_label -- a string placed in front of the text in results.txt
        
        stop_beg_TS_optimization -- a boolean. If True, will allow for manual checking of the TS molecule before
        optimization begins. Default False
        
        stop_beg_R_optimization -- a boolean. If True, will allow for manual checking of the R molecule before
        optimization begins. Default False
        
        stop_beg_freq -- a boolean. If True, the program will pause after calculating the first
        vibrational frequencies, allowing the user to view them, before calling Isoeff98.
        Default False.
        
        stop_iteration_optimization -- a boolean. If True, will allow for manual checking of the TS molecule at
        every step of the optimization process. Default False.
        
        stop_iteration_frequency -- a boolean. If True, the program will pause after calculating the
        vibrational frequencies at every step of the optimization process,
        allowing the user to view them, before calling Isoeff98.
        Default False.
        
        increment -- a boolean. If True, the program will run run_next_steps until self.Increment has been reached.
        This will occur before optimization begins. Default False.
        
        rerun -- a boolean. If True, the program will assume that TS has been run before. It is to be set to True
        if TS faulted from too low a value on MaxCycles, and the TS checkpoint file is still updated.
        Assumes that R has not been run before. Default False.
        
        just_minimize -- a boolean. If True, the program will assume initiate_files has already been run. This
        assumes that TS.log, TS.chk, R.log, R.chk, and self.current_distances are valid and
        updated. This option will override rerun. Default False.
        

        Returns: None
        """
        try:
            distance, errors, new_errors = self.run_program(stop_beg_TS_optimization=stop_beg_TS_optimization,\
                                                            stop_beg_R_optimization=stop_beg_R_optimization, \
                                                            stop_beg_freq=stop_beg_freq, stop_iteration_optimization=stop_iteration_optimization,\
                                                            stop_iteration_frequency=stop_iteration_frequency, increment=increment, \
                                                            rerun=rerun, just_minimize=just_minimize)
            try:
                lines = get_lines('results.txt')
            except IOError:
                subprocess.call(['cp', 'R.gjf', 'results.txt'])
                lines = []
            lines.extend([results_label, str(distance) + '\n', str(errors) + '\n', 'Isoeff Outputs, inputs\n', \
                     str(self.isoeff_output) + '\n', str(self.inputs) + '\n', 'Output to be minimized\n', \
                          str(self.minimized_output) + '\n'])
            write_lines('results.txt', lines)
            if os.getenv('DISPLAY') is not None:
                draw_graphs('Temporary results')
            print "Final distance(s)", distance
            print "Error between experimental and calculated Frequencies", errors
            print "Are the errors found within experimental bounds", new_errors
        except KeyboardInterrupt:
            if len(self.inputs) != len(self.isoeff_output) and len(self.inputs) -1 == len(self.isoeff_output):
                self.inputs = self.inputs[:-1]
            lines = ['Outputs, inputs\n', str(self.isoeff_output) + '\n', str(self.inputs) + '\n', 'Output to be minimized\n', \
                     str(self.minimized_output) + '\n']
            write_lines('results.txt', lines)
            if os.getenv('DISPLAY') is not None:
                draw_graphs('results.txt')

        
##test = Experiment('Methylchloride_Substrate.gjf.txt', ['1 5 2.5 F'], 'file2', ['hiagain'], 'isofile', [1.0, 2.0, 3.0])
##test.gauss_view_conversion('Methylchloride_Substrate.gjf.txt', ['1 5 2.5 F'], 'TS')

##if __name__ == '__main__':
##    print "Hello world"
##    x = 1
##    print x + 2
##    test = Experiment('Methyl_Formate_R.gjf', ['4 5 1.9 F'], 'Methyl_Formate_TS.gjf',\
##                      ['4 5 1.9 F'], 'Isoeff_MethylFormate.txt', [(1.009, 0.001)])
##    #faulttest = Experiment('nonexistant', ['string'], 'nonexistant', ['string'], \
##    #                       'isoeff nonexistant', [(1, 2)])
##    print "getStatus"
##    print ''
##    for foo in test.getStatus():
##        print foo
##    print ''
##    print 'gauss_view_conversion, mult'
##    for line in test.gauss_view_conversion('R'):
##        print line
##    print '____________________________________'
##    print 'gauss_view_conversion, justopt'
##    for line in test.gauss_view_conversion('TS', StopInOptFreq=True):
##        print line
##    print '____________________________________'
##    for line in test.write_freq_file('R'):
##        print line
##    print '____________________________________'
##    for line in test.rewrite_gauss_view([0.5]):
##        print line
##    print "________"
##    for line in test.rewrite_gauss_view([0.5], StopInOptFreq=True):
##        print line
##lines = get_lines('get_structure_test.txt')
##for line in get_structure(lines):
##    print line
