import string
import subprocess
import numpy
#from scipy import *
from scipy.optimize import fmin
import os
import itertools

#TODO: if this is inported into the main file, will the "if" statement be run?

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
