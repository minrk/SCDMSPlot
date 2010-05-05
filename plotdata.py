#!/usr/bin/env python
"""
Script for parsing and plotting SCDMS deposition data

For use in command-line:
$> python /path/to/plotdata.py R*.txt

or within Python:
>>> plot_run("filename.txt") # to plot single file (plot_run(fname, True) to save the result)
>>> multiplot("file*.txt") # to plot each file that matches a pattern

by Min Ragan-Kelley (minrk@berkeley.edu)
for Matt Cherry (mcherry1@stanford.edu) 
"""

######### adjustable preferences #########

t_units = "min" # options are min,sec,scan (scan is default)
format = "pdf" #also does image formats, i.e. png,jpg, etc.
# lines to plot:
default_lines = None # if this is None, everything will be plotted
default_lines = [2,4]
# use default_lines = [2,28,123,45] etc. to change the default set of lines to be plottedh


# the labels:
#   keyed by weight, each value can either be a string: the label in the legend
#   or a length 2 or 3 tuple
#   if length is 2, first is the label, second is a matplotlib style string, i.e. '-' or 'g:' (see pylab.plot for full list)
#   if length is 3, first is the label, second is style, third is color
#           this allows for extra colors, like 'chartreuse' or '#abc123'
#           see matplotlib.colors.cnames for named colors
#           or use rgb hex values with '#aabbcc'
the_labels = {
2 : '$H_2$',
14: '$N$',
15: '$CH_3$',
16: '$O$',
18: '$H_2O$',
28: '$N_2$',
32: '$O_2$',
40: ('$Ar$','g--'),
44: '$CO_2$',
'TPressure': ('$Total$', '-', 'darkred')
}

# this is the list of colors to use for lines:
the_colors = ['b', 'g', 'r', 'c', 'm', 'k', 'orange','brown', 'y', 'pink', 'slategray']
the_styles = [':','-.','--', '-']
# see matplotlib.colors.cnames for named colors. '#RRGGBB' also works

######## end preferences  ########

import re
from glob import glob
from os import path
import pylab
from numpy import *

try:
    pylab.rcParams['axes.color_cycle'] = the_colors
except:
    # for mpl versions < 1.0
    import matplotlib
    matplotlib.axes.set_default_color_cycle(the_colors)

pylab.rcParams['lines.markeredgewidth'] = 0.2
pylab.rcParams['legend.fontsize']='small'
pylab.rcParams['legend.labelspacing']=0.


class LabelDict(dict):
    def __init__(self, dikt):
        self.dikt = dikt
        dict.__init__(self, dikt)
    
    def __getitem__(self,key):
        if self.dikt.has_key(key):
            return self.dikt[key]
        else:
            return '$%s$'%str(key)

the_labels = LabelDict(the_labels)

def lookup(val,dikt):
    for k,v in dikt.iteritems():
        if v == val:
            return k

def update_me():
    """updates the plotdata.py script from git"""
    import os,shutil,time,sys
    ans = raw_input("are you sure you want to overwrite this file? (y/n): ").lower()
    if "y" in ans and "n" not in ans:
        shutil.move(__file__, strip_extension(__file__)+'.prev.py')
        cmd = "curl http://github.com/minrk/SCDMSPlot/raw/master/plotdata.py -o '%s'"%__file__
        print "Performing: %s"%cmd
        sys.stdout.flush()
        # time.sleep(1)
        os.system(cmd)
        print "don't forget to 'run %s' to load the updates"%__file__
    else:
        print "cancelling"

def strip_extension(s):
    """remove the extension from a filename:
    file.txt -> file
    foo.bar.ex -> foo.bar
    foobar -> foobar"""
    sr = s[::-1]
    return s[:-sr.find('.') - 1 or None]

def parse(fname,keep=None):
    """parses a run file, i.e. RXXX.txt
    >>> parse('R123.txt')
    will parse the file, keeping all known labels
    
    >>> parse('R123.txt' keep=[2,18])
    will parse the file, but keep only masses 2 and 18
    >>> parse('R123.txt',keep=['H2','Ar'])
    will keep H2 and Ar, assuming they are in the weights dictionary
    """
    if keep is None:
        keep = default_lines
    times = []
    M = []
    with open(fname) as f:
        # skip the header
        s = f.readline()
        while s and not s.startswith("SOD"):
            s = f.readline()
        if not s: # reached EOF
            raise IOError("Could not parse %s, incorrect format?"%fname)
        sod_prefix = s.split("\\")[-1].strip()
        sod_prefix = sod_prefix.replace(".",'-')
        parent = path.dirname(fname)
        sod_fname = path.join(parent, sod_prefix+"V.txt")
        try:
            sod = parse_sod(sod_fname)
            print "  SOD: %s"%sod_fname
        except Exception, e:
            print "could not parse SOD %s"%sod_fname
            sod=None
            
        while s and not s.startswith("Scan#"):
            s = f.readline()
        header = s
        
        # Scan#	Time Into Run	MASS( 2)	MASS(18)	MASS(28)	MASS(32)	MASS(40)	MASS(44)	TPressure
        header = s.split('\t')
        # masses = map(int, [e.replace('Mass(','')])
        kept   = []
        weights = []
        for i,ms in enumerate(header[2:]): # ignore first 2, ecause they are time
            match = re.findall('[0-9]+',ms)
            if match:
                # get the number
                z = int(match[0])
            else:
                # string name, probably just TPressure
                z = ms.strip()
            label = the_labels.get(z, '')
            if not isinstance(label, str):
                label = label[0]
            if keep is None or isinstance(z,str) or z in keep or label in keep\
                        or label.replace('$','').replace('_','') in keep:
                weights.append(z)
                kept.append(i)
        
        kept.reverse()
        weights.reverse()
        
        if not s: # EOF
            raise IOError("Could not parse %s, incorrect format?"%fname)
        for lines in f.readlines():
            try:
                line = lines.split()
                ts = line[1]
                hr,m,sec = map(float, ts.split(':'))
                floats = map(float, line[2:])
                MM = [floats[i] for i in kept] # keep only requested
            except Exception, e:
                # print e, lines
                pass # skip over the end lines
            
            else:
                # save the data
                t = 3600*hr+60*m+sec
                times.append(t)
                M.append(MM)
    # M.reverse()
    return array(times),array(M), weights, sod

def feed(f):
    """skip over blank lines
    return first non-empty line"""
    s = f.readline()
    while s and not s.strip():
        s = f.readline()
    return s

def parse_sod(fname):
    """parses SOD file (i.e. feb2003-xxxV.txt) to get overlays and title"""
    wstart = rfstart = rfstop = wstop = 0
    # print fname
    runmatch = re.findall('[0-9]+V\.',fname)
    if runmatch:
        num = re.findall('[0-9]+',runmatch[0])[0]
    else:
        num = strip_extension(fname)
    with open(fname) as f:
        s = f.readline()
        while s and not re.findall("run\W*[0-9]+\.?[0-9]+",s.lower()):
            s = f.readline()
        datestr = s.split()[-1]
        header = feed(f).strip()
        # if 'batch' in s.lower
        while s and not re.findall("[0-9]\t",s):
            s = f.readline()
            if s.strip() and not re.match("^[0-9]\t",s):
                header += " "+s.strip()
        if not s: # reached EOF
            raise IOError("Could not parse %s, incorrect format?"%fname)
        
        firstline = s
        matches = re.findall("G[A-Z0-9]+",header)
        # if matches:
        title = "R%s (%s)"%(num, ", ".join(matches))
        # else:
        
        # header = header.replace('SCDMS','').strip()
        # now let's try to parse the header into a title...eesh
        # batch:
        # print datestr, header
        # return
        # matches = re.findall("batch\W*[0-9]+",header, flags=re.IGNORECASE)
        # print title
        # title = s[5:].strip() # drop the SCDMS prefix
        # # include next line for split titles:
        # title += " "+f.readline()
        # title = title.strip()
        
        
        
        for line in [firstline]+f.readlines():
            try:
                step,s = line.split('\t',1)
                step = int(step)
            except:
                pass
            else:
                s = s.lower().strip()
                if s.startswith("rf") or "etch" in s:
                     if rfstart and 'off' in s:
                         rfstop = step
                     else:
                         rfstart = step
                if rfstop and s.startswith("w"):
                    if 'dep' in s:
                        wstart = step
                    elif wstart and 'off' in s:
                        wstop = step
    if not (rfstart and rfstop and wstart and wstop):
        print "Could not fully parse SOD: %s, defaulting some values"%fname
    
    if not rfstop:
        rfstop = rfstart+50
    if not wstart:
        wstart = rfstart+100
    if not wstop:
        wstop = wstart+20
    return title,rfstart,rfstop,wstart,wstop
                    


def plot_run(fname,styles=None,keep=None,save=False,hold=False, cmap_plot=False,align_legend=True):
    """parse a file and plot it
    Also takes style, keep, and save keywords
    styles: str or list of str
        see pylab.plot for line styles, e.g. ':', '-o'
        default: solid line
    keep: a list of ints or strings
        filter for parsing the masses
        Will only plot lines for elements whose weight or name is in keep
        Always keeps the total pressure
        Default: plot all
        
        For instance: keep=['Ar','H2']
        or:           keep=[2,32]
    
    save: bool
        set whether to save or not at the end of the plot
        default: False
    
    cmap_plot: bool
        set whether to use a heat map to color lines, as opposed to the usual color cycling.
        default: False
    
    align_legend: bool
        set whether to sort the labels in the legend, such that the highest lines are the first labeled.
        If False, legend will be sorted by molecular weight.
        default: True
        
        
    >>> plot_run('R123.txt',keep=['Ar',2])
    >>> plot_run('R123.txt',':') # for all, but dotted
    """
    #
    ncolors = len(the_colors)
    
    if styles is None:
        styles = the_styles # accomodates 35 lines without repeating
    elif isinstance(styles,str):
        styles = [styles]
    styles=list(styles) # in case of arrays that won't grow when I do styles*2
    
    t,M,weights,sod = parse(fname,keep=keep) # get the data
    labels = []
    plotstyles = []
    colors = []
    
    # for w,tup in zip(weights, rich_labels):
    for w in weights:
        tup = the_labels[w]
        if isinstance(tup,str):
            l = tup
            s = styles[(hash(w)/len(the_colors))%len(styles)]
            c = the_colors[hash(w)%len(the_colors)]
            # print tup, w,s,c
        elif len(tup) == 2:
            l,s = tup
            c = None
        else:
            l,s,c = tup
        labels.append(l)
        plotstyles.append(s)
        colors.append(c)
        
    if 'min' in t_units.lower():
        t = t/60
    elif 'sec' in t_units.lower():
        pass
    else: # scan
        t = arange(0,len(t))
    
    
    if not hold:
        pylab.figure() # new figure
    fig = pylab.gcf()
    fig.subplotpars.right=0.82 # pad right for placing 
    
    if sod: # if we found a SOD file to parse
        times = array(sod[1:])-1 # adjust for index starting at 1 in SOD
        rfstart,rfstop,wstart,wstop = times
        t -= t[times[0]] # 
    
    # do the plot:
    # split in sets of ncolors for color repeating
    if not cmap_plot:
        if align_legend:
            MM = zeros((M.shape[0]+1,M.shape[1]))
            MM[1:] = M
            MM[0] = range(M.shape[1])
            lM = list(MM.transpose())
            sM = sorted(lM,key=lambda a: a[1:].max(),reverse=True)
            M2 = array(sM).transpose()
            M = M2[1:]
            
            neworder = map(int, M2[0])
            labels = [labels[i] for i in neworder]
            weights = [weights[i] for i in neworder]
            plotstyles = [plotstyles[i] for i in neworder]
            colors = [colors[i] for i in neworder]
        
        for s,c,line in zip(plotstyles,colors, M.transpose()):
            if c is not None:
                pylab.semilogy(t,line,s,color=c)
            else:
                pylab.semilogy(t,line,s)
                
    else:
        while ',' in styles:
            styles.remove(',') # drop the pixel
        # color mapped plot:
        colors = pylab.cm.ScalarMappable(cmap=pylab.cm.jet)
        colors.set_clim(0,M.shape[1])
        for i,line in enumerate(M.transpose()):
            pylab.semilogy(t,line,styles[i%len(styles)],color=colors.to_rgba(i))
    # for style in styles:
    
    if sod: # title and overlays from SOD file (if we have it)
        pylab.title(sod[0])
        overlay(t[rfstart],t[rfstop])
        overlay(t[wstart],t[wstop])
    else:
        pylab.title(fname)
    
    # draw the legend
    # leg = pylab.legend(labels,loc=(1.0,.25))
    leg = pylab.legend(labels,loc=(1.0,0))
    # leg.get_frame().set_alpha(0.75) # for legend transparency
    
    # shrink-to-fit and label x-axis:
    pylab.xlim(t.min(),t.max())
    
    if 'min' in t_units.lower():
        pylab.xlabel("t(min)")
    elif 'sec' in t_units.lower():
        pylab.xlabel("t(s)")
    else:
        pylab.xlabel("Scan (6s interval)")
    
    pylab.ylabel("Amps, Torr for Total")
    pylab.grid(True) # turn on the grid
    if save: # save to a file
        title = strip_extension(fname)
        pylab.savefig(title+'.'+format)

def overlay(left,right,color='grey'):
    """overlay a transparent rectangle between @left and $right, filling vertically."""
    ax = pylab.gca()
    y1,y2 = pylab.ylim()
    ax.fill_between([left,right],[y2*.99]*2,[y1*1.01]*2,color=color,alpha=0.25)


def multiplot(*matches,**kwargs):
    """for instance, multiplot('R*') or multiplot('R*.txt') or multiplot('*.txt','*.TXT')
    This method always saves the plot to a file with the same prefix as the source data."""
    # if not isinstance(matches, list):
        # matches=[matches]
    for match in matches:
      for fname in glob(match):
        print fname
        try:
            plot_run(fname,save=True,**kwargs)
        except Exception, e:
            print "skipping %s, plot failed due to:"%(fname)
            print e


if __name__ == '__main__':
    """when called from the command-line, just passes all arguments to multiplot"""
    import sys
    # for match in sys.argv[1:]:
        # print sys.argv[1:]
    multiplot(*sys.argv[1:])


