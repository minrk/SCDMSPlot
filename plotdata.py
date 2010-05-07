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
# lines to plot:  (Comment one of these out!)
default_lines = None # if this is None, everything will be plotted
#default_lines = [2,18,28,32,40,44] # This is the "standard suspects" default set

# the labels:
#   keyed by weight, each value can either be a string: the label in the legend
#   or a length 2 or 3 tuple
#   if length is 2, first is the label, second is a matplotlib style string, i.e. '-' or 'g:' (see pylab.plot for full list)
#   if length is 3, first is the label, second is style, third is color
#           this allows for extra colors, like 'chartreuse' or '#abc123'
#           see matplotlib.colors.cnames for named colors
#           or use rgb hex values with '#aabbcc'
the_labels = {
2 : ('$H_2$','-','crimson'),
18: ('$H_2O$','-','aqua'),
28: ('$N_2$','-','darkgreen'),
32: ('$O_2$','-','gold'),
40: ('$Ar$','-','darkred'),
44: ('$CO_2$','-','chocolate'),
'TPressure': ('$Total$','-','purple')
}

# this is the list of colors to use for lines:
the_colors = ['b', 'g', 'r', 'c', 'm', 'orange', 'y', 'pink', 'slategray']
the_styles = [':','-.','--']
# see matplotlib.colors.cnames for named colors. '#RRGGBB' also works

markersize = 4. # default is 6.

######## end preferences  ########

import re,os,shutil,sys

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

pylab.rcParams['lines.markersize'] = markersize
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
    print "Updating: ''%s"%__file__
    ans = raw_input("Are you sure you want to overwrite this file? (y/n): ").lower()
    
    cwd = path.abspath(os.curdir)
    gitdir = path.abspath(path.dirname(__file__))
    
    if "y" in ans and "n" not in ans:
        prevfile = strip_extension(__file__)+'.prev.py'
        shutil.copy(__file__, prevfile)
        if os.path.isdir(path.join(gitdir, '.git')):
            git=True
            os.chdir(gitdir)
            cmd = 'git pull'
        else:
            git=False
            cmd = "curl http://github.com/minrk/SCDMSPlot/raw/master/plotdata.py -o '%s'"%__file__
        print "Performing: %s"%cmd
        sys.stdout.flush()
        if not os.system(cmd):
            if not git:
                print "previous version located at %s"%prevfile
            else:
                os.chdir(cwd)
            print "don't forget to 'run %s' to load any updates"%__file__
        else:
            print "Download failed!"
    else:
        print "cancelling update"

def push_to_git():
    """pushes your changes to git"""
    ans = raw_input("Are you sure you want to commit your changes to git? (y/n): ").lower()
    
    cwd = path.abspath(os.curdir)
    gitdir = path.abspath(path.dirname(__file__))
    
    if "y" in ans and "n" not in ans:
        msg = raw_input("change summary: ")
        if not msg.strip():
            msg = 'interactive update'
            
        os.chdir(gitdir)
        cmd = "git add plotdata.py && git commit -m \"%s\" && git push"%msg
        print "Performing: %s"%cmd
        sys.stdout.flush()
        if not os.system(cmd):
            'git push failed'

        else:
            print "Upload failed!"
        os.chdir(cwd)
    else:
        print "cancelling update"
    

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
    if isinstance(keep, (str,int)):
        keep = [keep]
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
            # print "  SOD: %s"%sod_fname
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
            if keep is None or z in keep or label in keep\
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
                    

def setup_styles(weights, style):
    orig_style = style
    if style is None:
        style = the_styles # accomodates 35 lines without repeating
        wasNone = True
    elif isinstance(style,str):
        style = [style]
    
    labels = []
    plotstyles = []
    colors = []
    
    # for w,tup in zip(weights, rich_labels):
    for w in weights:
        tup = the_labels[w]
        if isinstance(tup,str):
            l = tup
            s = style[(hash(w)/len(the_colors))%len(style)]
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
        
    #
    if orig_style is not None:
        colors = [None]*len(plotstyles)
        while len(style) < len(plotstyles):
            style *= 2
        plotstyles = style
        colors = [None]*len(plotstyles)
    
    return labels,plotstyles,colors

def reorder(M,*also):
    MM = zeros((M.shape[0]+1,M.shape[1]))
    MM[1:] = M
    MM[0] = range(M.shape[1])
    lM = list(MM.transpose())
    sM = sorted(lM,key=lambda a: a[1:].max(),reverse=True)
    M2 = array(sM).transpose()
    M = M2[1:]
    reordered = [M]
    neworder = map(int, M2[0])
    for lis in also:
        reordered.append( [lis[i] for i in neworder ] )
    
    return reordered
    

def plot_run(fname,style=None,keep=None,save=False,hold=False, cmap_plot=False,align_legend=True,use_colors=True):
    """parse a file and plot it
    Also takes style, keep, and save keywords
    style: str or list of str
        see pylab.plot for line style, e.g. ':', '-o'
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
    
    t,M,weights,sod = parse(fname,keep=keep) # get the data
    labels = []
    plotstyles = []
    colors = []
    
    # for w,tup in zip(weights, rich_labels):
    labels,plotstyles,colors=setup_styles(weights, style)
    
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
            M,labels,weights,plotstyles,colors = reorder(M,labels,weights,plotstyles,colors)
        
        for s,c,line in zip(plotstyles,colors, M.transpose()):
            if c is not None:
                pylab.semilogy(t,line,s,color=c)
            else:
                pylab.semilogy(t,line,s)
                
    else:
        while ',' in style:
            style.remove(',') # drop the pixel
        # color mapped plot:
        colors = pylab.cm.ScalarMappable(cmap=pylab.cm.jet)
        colors.set_clim(0,M.shape[1])
        for i,line in enumerate(M.transpose()):
            pylab.semilogy(t,line,style[i%len(style)],color=colors.to_rgba(i))
    # for style in style:
    
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
    print M.min()
    realmin = (M+M.max()*(M==M.min())).min()
    print realmin
    pylab.ylim(ymin=10**floor(log(realmin)/log(10)))
    
    if 'min' in t_units.lower():
        pylab.xlabel("t(min)")
    elif 'sec' in t_units.lower():
        pylab.xlabel("t(s)")
    else:
        pylab.xlabel("Scan (6s interval)")
    
    
    
    pylab.ylabel("Amps, Torr for Total")
    pylab.grid(True) # turn on the grid
    if save: # save to a file
        if not isinstance(save, str):
             save = strip_extension(fname) + '.' +format
        pylab.savefig(save)

def timeslice(files, offset=0,radius=1, style=None,keep=None,save=False,hold=False, align_legend=True,normalize=False,return_data=False):
    """Parse several files, and plot a single timeslice offset from rfstart for each of them.
    
    files: str or list of strs
        either filename(s): ['R123.txt', 'R555.txt']
        or glob(s): ['R22*.txt', 'R3*']
        or: 'R*.txt'
        
    offset: int
        the slice to plot. units are scans from rfstart
        offset=10: 10 scans after rfstart
        offset=-20: 20 scans before rfstart
    
    radius: int
        the radius about @offset
        rfstart+offset +/- radius will be averaged
        default: 1
    
    Also takes style, keep, save, align_legend keywords
    
    style: str or list of str
        see pylab.plot for line style, e.g. ':', '-o'
        default: solid line
    keep: a list of ints or strings
        filter for parsing the masses
        Will only plot lines for elements whose weight or name is in keep
        Always keeps the total pressure
        Default: plot all
        
        For instance: keep=['Ar','H2']
        or:           keep=[2,32]
        or:           keep=['Ar', 2, '$H_2$']
    
    save: bool
        set whether to save or not at the end of the plot
        default: False
    
    align_legend: bool
        set whether to sort the labels in the legend, such that the highest lines are the first labeled.
        If False, legend will be sorted by molecular weight.
        default: True
        
    """
    the_files = []
    if isinstance(files, str):
        files = [files]
    # map(the_files.extend)
    for f in files:
        the_files.extend(glob(f))
    
    # keep=kwargs.get('keep',None)
    names = []
    lines = []
    for i,fname in enumerate(the_files):
        t,M,weights,sod = parse(fname,keep=keep) # get the data
        name = strip_extension(path.basename(fname))
        if sod: # if we found a SOD file to parse
            times = array(sod[1:])-1 # adjust for index starting at 1 in SOD
            rfstart,rfstop,wstart,wstop = times
            t -= t[times[0]] # 
        else:
            print "skipping %s"%name
            continue
        
        A = M[rfstart+offset-radius:rfstart+offset+radius+1].transpose()
        lines.append(map(mean, A))
        names.append(name)
    M = array(lines)
    
    labels,plotstyles,colors=setup_styles(weights, style)
    #
    if not hold:
        pylab.figure() # new figure
    fig = pylab.gcf()
    fig.subplotpars.right=0.82 # pad right for placing 
    
    if align_legend:
        M,labels, weights,plotstyles,colors = reorder(M,labels, weights,plotstyles,colors)
    
    x = range(len(M))
    for s,c,line in zip(plotstyles,colors, M.transpose()):
        if normalize:
            m = line.mean()
            line = line/(m+1e-20)
        try:
            if c is not None:
                pylab.semilogy(x, line,s+'o',color=c)
            else:
                pylab.semilogy(x,line,s+'o')
        except ValueError: # there was already a point-style specified?
            if c is not None:
                pylab.semilogy(x, line,s,color=c)
            else:
                pylab.semilogy(x,line,s)
    
    leg = pylab.legend(labels,loc=(1.0,0))
        
    title = "RF start + %i"%(offset)
    if radius:
            title += "$\pm$ %i"%(radius)
    pylab.title(title)
    pylab.xticks(x,names,fontsize='small')
    pylab.xlim(-1,len(names))
    pylab.ylabel("Amps, Torr for Total")
    pylab.grid(True) # turn on the grid
    if save: # save to a file
        if not isinstance(save, str):
             save = 'C-Dep Scatter.%i'%offset+'.'+format
        pylab.savefig(save)
    if return_data:
        return M.transpose()


def single_gas(files, keep=2, offsets=[20], radius=1, style=None):
    """runs timeslice for one gas at various times
    call it like:
    >>> single_gas("R*.txt", keep=4, offsets=[0,20,40,72], radius=5)
    there will be one line for each entry in @offsets
    """
    the_files = []
    if isinstance(files, str):
        files = [files]
    # map(the_files.extend)
    for f in files:
        the_files.extend(glob(f))
    
    if style is None:
        colors = 'r g b k c m y'.split()
        styles = ('-','--', '-.', ':')
        plotstyles = []
        for s in styles:
            for c in colors:
                plotstyles.append( c+s )
    else:
        if isinstance(style, str):
            style = list(style)
        plotstyles = style
        while len(plotstyles) < len(offsets):
            plotstyles.extend(style)
    
    pylab.figure()
    for off,s in zip(offsets,plotstyles):
        timeslice(files, style=s, keep=keep, offset=off, radius=radius, hold=True)
    pylab.legend(["rfstart$+%i\pm%i$"%(i,radius) for i in offsets ])
    label = the_labels[keep]
    if not isinstance(label, str):
        label = label[0]
    pylab.title(label)
    

def outliers(lines, tol=.25):
    """filters a set of lines, keeping only those with points with outliers beyond a tolerance"""
    keep = []
    for i,line in enumerate(lines):
        m = line.mean()
        diffs = abs((line-m)/m)
        if (diffs > tol).any():
            keep.append(i)
    return keep

def overlay(left,right,color='grey'):
    """overlay a transparent rectangle between @left and $right, filling vertically."""
    ax = pylab.gca()
    y1,y2 = pylab.ylim()
    ax.fill_between([left,right],[y2*.99]*2,[y1*1.01]*2,color=color,alpha=0.25)


def multiplot(*matches,**kwargs):
    """
    This is a wrapper for plot_run on multiple files. See plot_run for keyword arguments (kwargs).
    
    For instance: multiplot('R*') or multiplot('R*.txt') or multiplot('*.txt','*.TXT')
    
    This method always saves the plot to a file with the same prefix as the source data.
    
    """
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
    for arg in sys.argv[1:]:
        if 'update' in arg.lower():
            update_me()
    multiplot(*sys.argv[1:])


