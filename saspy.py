# python lib
'''
SASpy - ATSAS PLUGIN FOR PYMOL

(c) 2015-2016 A.PANJKOVICH FOR ATSAS TEAM AT EMBL-HAMBURG.
'''
import os
import sys
import re
import time
import shutil
import string
import Tkinter
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
import tempfile
import threading
import subprocess
import itertools


# pymol lib
try:
    from pymol import cmd
    from pymol.cgo import *
except ImportError:
    print 'Warning: pymol library cmd not found.'
    sys.exit(1)

# external lib
try:
    import Pmw
except ImportError:
    print 'Warning: failed to import Pmw. Exit ...'
    sys.exit(1)

## Plugin initialization

#global variables
currentDat = []
modelingRuns = 0
datViewer = Tkinter.StringVar()
cwd = Tkinter.StringVar()
cwd.set(os.getcwd())

from sys import platform
if "win32" == platform:
    datViewer.set("sasplotqt")
else:
    datViewer.set("sasplot") #on linux + mac

def __init__(self):
    """ SASpy - ATSAS Plugin for PyMOL
    """
    self.menuBar.addmenuitem('Plugin', 'command',
                             'SASpy', label = 'SASpy',
                             command = lambda s=self : SASpy(s))


class TemporaryDirectory:
    """Context Manager for working in a temporary directory"""

    def __init__(self, *args, **kwargs):
        self.temp_dir = tempfile.mkdtemp(*args, **kwargs)

    def __enter__(self):
        self.orig_dir = os.getcwd()
        os.chdir(self.temp_dir)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.orig_dir)
        # If there was an error, do not delete the temporary
        # directory, so that the user can examine its contents
        if exc_type is None:
            shutil.rmtree(self.temp_dir)

    def copy_in(self, src, dst=None):
        """Copy a file into the temporary directory

        Arguments:
        src -- Source file name (relative to the original working directory)
        dst -- Destination file name (relative to the temporary directory)
               If not present, same as the source file name
        """
        if dst is None:
            dst = src
        if os.path.isabs(dst):
            raise ValueError("Destination path should not be absolute")
        abs_src = os.path.join(self.orig_dir, src)
        abs_dst = os.path.join(self.temp_dir, dst)
        shutil.copy(abs_src, abs_dst)
        return abs_dst

    def move_out(self, src, dst=None):
        """Move a file out of the temporary directory

        Arguments:
        src -- Source file name (relative to the temporary directory)
        dst -- Destination file name (relative to the original directory)
               If not present, same as the source file name
        """
        if os.path.isabs(src):
            raise ValueError("Source path should not be absolute")
        if dst is None:
            dst = src
        abs_src = os.path.join(self.temp_dir, src)
        abs_dst = os.path.join(self.orig_dir, dst)
        shutil.move(abs_src, abs_dst)
        return abs_dst

    def move_out_numbered(self, src, prefix, suffix):
        """Move a file out of the temporary directory, without overwriting old files

        Chooses a new file name based on the given prefix and suffix and a unique number

        Arguments:
        src -- Source file name (relative to the temporary directory)
        prefix -- prefix for the destination filename
        suffix -- suffix for the destination filename
        """
        if os.path.isabs(src):
            raise ValueError("Source path should not be absolute")
        if suffix.startswith('.'):
            suffix = suffix[1:]
        abs_src = os.path.join(self.temp_dir, src)
        abs_prefix = os.path.join(self.orig_dir, prefix)
        abs_dst = "%s.%s" % (abs_prefix, suffix)
        if os.path.exists(abs_dst):
            for i in itertools.count(start=1):
                abs_dst = "%s_%d.%s" % (abs_prefix, i, suffix)
                if not os.path.exists(abs_dst):
                    break
        shutil.move(abs_src, abs_dst)
        return abs_dst


## class for GUI
class SASpy:

    def __init__(self, app):

        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Quit',
                                 #'Debug',
                                 'Refresh model list',
                                 '3. Execute'),
                                 title = 'SASpy - ATSAS Plugin for PyMOL',
                                 command = self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        self.atsasThreads  = []
        self.maxThreads    = 1
        self.procedure     = 'empty'
        self.saxsfn        = Tkinter.StringVar()
        self.crysolmode    = Tkinter.StringVar()
        self.crysolmode.set('predict')
        self.prefix        = Tkinter.StringVar()

        self.warnLabel = Tkinter.Label( self.dialog.interior(),
                                    anchor='center',
                                    fg ="red",
        text = 'No models found. Please open/load structures in PyMOL to proceed.')

        w = Tkinter.Label(self.dialog.interior(),
                          text = '\nSASpy - ATSAS Plugin for PyMOL\nVersion 1.3 - ATSAS 2.7.1\n\nEuropean Molecular Biology Laboratory\nHamburg Outstation, ATSAS Team, 2015-2016.\n',
                          background = 'white', foreground = 'blue')
        w.pack(expand = 1, fill = 'both', padx = 10, pady = 5)

        self.procLabel = Tkinter.Label( self.dialog.interior(),
                                        anchor='w',
                                        text = '1. Choose procedure:')
        self.procLabel.pack(fill='both', expand=True, padx=10, pady=5)

        #NOTEBOOK START
        self.notebook = Pmw.NoteBook(self.dialog.interior(),raisecommand=self.tabSelection)
        self.notebook.pack(fill = 'both', expand=2, padx=10, pady=10)

        # crysol tab
        crysolTab = self.createTab("crysol", 
         "Prediction of theoretical intensities and optionally fit\n"+
         "to experimental SAXS data. Please select at least one model\n"+
         "(and a SAXS .dat file for fit mode)."
        )
        self.crymodebut = Pmw.RadioSelect(crysolTab,
                                    buttontype='radiobutton',
                                    labelpos='w',
                                    label_text="Mode:",
                                    command=self.setCrysolMode,
                                    selectmode = 'single')
        self.crymodebut.grid(sticky='we', row=2, column=0, padx=5, pady=2)
        self.crymodebut.add('predict')
        self.crymodebut.add('fit')
        self.crymodebut.setvalue('predict')

        saxsfn_ent = Pmw.EntryField(crysolTab,
                                    label_text = 'SAXS .dat file:',
                                    labelpos='ws',
                                    entry_textvariable=self.saxsfn)
        saxsfn_but = Tkinter.Button(crysolTab, text = 'Browse...',
                                    command = self.getSAXSFile)
        saxsfn_ent.grid(sticky='we', row=3, column=0, padx=5, pady=2)
        saxsfn_but.grid(sticky='we', row=3, column=1, padx=5, pady=2)

        fn = []
        fn.append('crysol')
#        self.notebook.setnaturalsize(pageNames=fn)

        #alpraxin tab
        alpraxinTab = self.createTab("alpraxin", "Position a structure such that its principal intertia\nvectors are aligned with the coordinate axis.\nPlease select one model.")

        #supalm tab
        supalmTab = self.createTab('supalm', 'Superimposition of models and calculation of\nnormalized spatial discrepancy (NSD).\nPlease select two models.')

        #sasref tab
        sasreftab = self.createTab("sasref", "Quaternary structure modeling against solution scattering data.\nPlease select multiple models (rigid bodies) and a SAXS .dat file.")
        saxsfn_ent = Pmw.EntryField(sasreftab,
                                    label_text = 'SAXS .dat file:',
                                    labelpos='ws',
                                    entry_textvariable=self.saxsfn)
        saxsfn_but = Tkinter.Button(sasreftab, text = 'Browse...',
                                    command = self.getSAXSFile)
        saxsfn_ent.grid(sticky='we', row=3, column=0, padx=5, pady=5)
        saxsfn_but.grid(sticky='we', row=3, column=1, padx=5, pady=5)

        # sreflex tab
        sreflexTab = self.createTab("sreflex", "Model refinement based on SAXS data and normal mode analysis.\nPlease select models and a SAXS .dat file.")
        saxsfn_ent = Pmw.EntryField(sreflexTab,
                                    label_text = 'SAXS .dat file:',
                                    labelpos='ws',
                                    entry_textvariable=self.saxsfn)
        saxsfn_but = Tkinter.Button(sreflexTab, text = 'Browse...',
                                    command = self.getSAXSFile)
        saxsfn_ent.grid(sticky='we', row=3, column=0, padx=5, pady=5)
        saxsfn_but.grid(sticky='we', row=3, column=1, padx=5, pady=5)


        #dam display tab
        self.damColor = Tkinter.StringVar();
        self.damColor.set('white');
        self.damTrans = Tkinter.StringVar();
        self.damTrans.set('0.5');
        damdisplayTab = self.createTab("damdisplay", "Apply a predefined representation to a dummy-atom-model (DAM).\nPlease select one model.")
        damDisplayColorEntry = Pmw.EntryField(damdisplayTab,
                                    label_text = 'Color:',
                                    labelpos='ws',
                                    entry_textvariable=self.damColor)
        damDisplayColorEntry.grid(sticky='we', row=3, column=1, padx=5, pady=5)
        damDisplayTransEntry = Pmw.EntryField(damdisplayTab,
                                    label_text = 'Transparency:',
                                    labelpos='ws',
                                    entry_textvariable=self.damTrans)
        damDisplayTransEntry.grid(sticky='we', row=4, column=1, padx=5, pady=5)

        # config tab
        configTab = self.createTab("configure", "Settings available to configure SASpy:")
        # saxs viewer selection
        svi_ent = Pmw.EntryField(configTab,
                                    label_text = 'SAXS viewer:',
                                    labelpos='ws',
                                    entry_textvariable = datViewer)
        svi_but = Tkinter.Button(configTab, text = 'Select SAXS viewer',
                                    command = self.getSAXSViewer)
        svi_ent.grid(sticky='we', row=3, column=0, padx=5, pady=5)
        svi_but.grid(sticky='we', row=3, column=1, padx=5, pady=5)

        #working directory selection
        wd_ent = Pmw.EntryField(configTab,
                                    label_text = 'Current working dir:',
                                    labelpos='ws',
                                    entry_textvariable = cwd)
        #wd_ent.grid(sticky='w', row=2, column=0, columnspan=3, padx=5, pady=2)
        wd_ent.grid(sticky='w', row=2, column=0, padx=5, pady=2)
        scd_but = Tkinter.Button(configTab, text = 'Select working directory',
                                    command = self.setWorkingDirectory)
        scd_but.grid(sticky='we', row=2, column=1, padx=5, pady=2)

        #Model selection
        self.modsW = self.createModelSelectionWidget()
        self.modsW.pack(expand=1, fill='both', padx=10, pady=5)
        self.refreshModelSelectionWidget()

        self.notebook.setnaturalsize()

#GUI FUNCTIONS

    def errorWindow(self,title, msg):
        tkMessageBox.showerror(title,
                               "ERROR "+ msg,
                               parent=self.parent)
        return

    def createModelSelectionWidget(self):
        modsW = Pmw.RadioSelect(self.dialog.interior(),
                                    buttontype='button',
                                    labelpos='w',
                                    label_text="2. Model selection:",
                                    selectmode = 'multiple')
        mols = self.getListOfModels()
        for m in mols:
            modsW.add(m)
        return modsW

    def countSelectedModels(self):
        counter = 0
        ma = self.modsW.getcurselection()
        for m in ma:
            counter += 1
        return counter

    def setCrysolMode(self, mode):
        print "Setting crysol mode to "+mode
        self.crysolmode.set(mode)
        self.crymodebut.setvalue(mode)

    def setDatMode(self, mode):
        print "Setting open mode to " + mode
        global datmode
        datmode.set(mode)
        self.datmodebut.setvalue(mode)

    def submitJobAsThread(self, procType, models= []):
        viewer = datViewer.get()
        #check how many threads from this plugin are running
        running = self.checkAtsasThreads()
        if running == self.maxThreads:
            self.errorWindow("Max threads exceeded",
                "A process is already running, please wait for it to complete.")
            return
        #check if saxs file is available
        saxsfn = self.saxsfn.get()
        if False == os.path.isfile(saxsfn):
            self.errorWindow("FILE NOT FOUND",
                             "SAXS file \'"+saxsfn+"\' NOT FOUND.");
            return
 
        if 'sasref' == procType:
            t = threading.Thread(target = sasref, name = procType+'_thread', args = (saxsfn, models, viewer))
        elif 'sreflex' == procType:
            t = threading.Thread(target = sreflex, name = procType+'_thread', args = (saxsfn, models, viewer))
        t.setDaemon(1)
        t.start()
        self.atsasThreads.append(t)
        return

    def checkAtsasThreads(self):
        threadsAlive = 0
        for t in self.atsasThreads:
            a = t.is_alive()
            name = t.name
            if True == a:
                threadsAlive += 1
            else:
                idx = self.atsasThreads.index(t)
                t.join()
                del self.atsasThreads[idx]
                print "Just removed "+name+" from the list, with index "+repr(idx)
        return threadsAlive

    def submitSaspyJob(self, procType, models = []):
        if "alpraxin" == procType:
            alpraxin(models[0])
            return
        if "crysol" == procType:
            self.crysol(models)
            return
        if "damdisplay" == procType:
            damdisplay(models[0], self.damColor.get(), self.damTrans.get())
            return
        if "supalm" == procType:
            supalm(models[0], models[1])
            return
        #remaining job types (sreflex, sasref) 
        #should be submitted as a thread
        self.submitJobAsThread(procType, models)
        return

    def prepareJobAndSubmit(self):
        procType = self.procedure
        seln = self.countSelectedModels()

        if 'configure' == procType:
            return

        if 0 == seln:
            self.errorWindow("No model selected", "Please select models")
            return

        #some procedures need an exact number of models selected
        expect_dict = { #expected number of models
                     'alpraxin':1,
                     'damdisplay':1,
                     'supalm':2,
        #other procedures need a minimum number of selected models
                     'sreflex':11, #subtract ten to obtain min expected
                     'crysol':11,
                     'sasref':12,
                    }
        expn = expect_dict[procType]

        if 10 > expn: #procedure needs an exact number of selected models
            if expn != seln:
                self.errorWindow("Wrong number of models selected",
                "You selected "+ getPlural(seln) + ", but \'" + procType+ "\' expects "+getPlural(expn)+".\n")
                return
        else: #the procedure needs a minimum number of selected models
            expn = expn - 10
            if seln < expn:
                self.errorWindow("Wrong number of models selected",
                "You selected "+ getPlural(seln) + ", but \'" + procType +
                "\' expects at least " + getPlural(expn) +".\n")
                return

        self.submitSaspyJob(procType, self.modsW.getcurselection())
        return        

    def getListOfModels(self):
        #models can not contain the underscore character '_'
        #this is a Pmw limitation, but PyMOL does add such
        #characters often
        initialList = cmd.get_object_list()
        outputList = list();
        for m in initialList:
            if '_' in m:
                newName = m.translate(None,"_")
                cmd.set_name(m,newName)
                message("WARNING Renaming model \'"+m+ "\' to \'"+ newName+"\'")
                m = newName
            outputList.append(m)
        return outputList

    def refreshModelSelectionWidget(self):
        self.modsW.deleteall()
        mols = self.getListOfModels()
        for m in mols:
            if '_' in m:
                newName = m
                newName.replace('_', '-')
                cmd.set_name(m,newName)
                m = newName
                message("Renaming a model...")
            self.modsW.add(m)
        if 1 > len(self.getListOfModels()):
            self.warnLabel.pack(fill='both', expand=True, padx=10, pady=5)
        else:
            self.warnLabel.pack_forget()
        return

    def createTab(self, name = 'empty', description = 'empty'):
        page = self.notebook.add(name)
        tab_struc = Tkinter.LabelFrame(page, text = name)
        tab_struc.pack(fill='both', expand=True, padx=10, pady=10)
        desc = Tkinter.Label(tab_struc, justify="left", text = description, pady=2)
        desc.grid(sticky='w', row=0, column=0, columnspan=4, padx=10, pady=10)
        return tab_struc

    def openCurrentDatFile(self):
        global currentDat
#        message("About to open current dat file: " + repr(currentDat))
        if 0 == len(currentDat):
            tkMessageBox.showerror('No curve yet',
                                   'No SAXS intensities have been calculated yet',
                                    parent=self.parent)

        else:
            openDatFile(datViewer.get(), currentDat)
        return

    def getSAXSViewer(self):
        global datViewer
        file_name = tkFileDialog.askopenfilename(
            title='SAXS viewer', initialdir='',
            parent=self.parent)
        datViewer.set(file_name)
        return

    def setWorkingDirectory(self):
        newWorkDir = tkFileDialog.askdirectory(
            title='Set working directory', initialdir='',
            parent=self.parent)
        cmd.cd(newWorkDir)
        cwd.set(newWorkDir)
        message("Working directory changed to: " + newWorkDir);
        return

    def getSAXSFile(self):
        if 'crysol' == self.procedure:
            self.setCrysolMode('fit')
        file_name = tkFileDialog.askopenfilename(
            title='SAXS File', initialdir='',
            filetypes=[('saxs data files', '*.dat'), ('all files', '*')],
            parent=self.parent)
        self.saxsfn.set(file_name)
        return

    def tabSelection(self, pagename):
        #refresh each time a tab is selected
        cwd.set(os.getcwd()) #I don't know how to refresh this if the 
                             #user just calls 'cd'
        self.procedure = pagename
        return

    def crysol(self, selection, param=""):
        #wrapper for the different crysol modes
        if 1 < len(selection):
            message("CRYSOL will be executed for a complex")
            message("made of the following models: "+repr(selection))
        crymode = self.crysolmode.get()
        if 'predict' == crymode:
            return updateCurrentDat(predcrysol(selection))
        elif 'fit' == crymode:
            saxsfn = self.saxsfn.get()
            if False == os.path.isfile(saxsfn):
                self.errorWindow("FILE NOT FOUND",
                             "SAXS file \'"+saxsfn+"\' NOT FOUND.");
                return
            return updateCurrentDat(fitcrysol(saxsfn, selection))
            
    def execute(self, cmd):
        """ Run the cmd represented by the button clicked by user.
        """
        if cmd == 'OK':
            print 'is everything OK?'

        elif cmd == 'Refresh model list':
            self.refreshModelSelectionWidget()

#        elif cmd == 'Debug':

        elif cmd == '3. Execute':
            self.prepareJobAndSubmit()

        elif cmd == 'Quit':
            self.checkAtsasThreads()
            for p in self.atsasThreads:
                print "WARNING, a thread is still running: " + repr(p.name)

            message('Quit')
            if __name__ == '__main__':
                self.parent.destroy()
            else:
                self.dialog.withdraw()

        else:
            print 'Terminating SASpy Plugin...'
            self.dialog.withdraw()
            print 'Done.'

##################
# CLI Funtions

defprefix = 'saspy_wd'

def systemCommand(command, **kwargs):
    status = subprocess.call(command, **kwargs)
    if(0 != status):
        message("WARNING, something went wrong while executing:\n"
                + ' '.join(command))
    return status

def message(text):
    print "SASpy: "+text
    return

def getPlural(n):
    #get plurals right
    outstring = repr(n) +" model";
    if 1 != n: 
        outstring+='s'
    return outstring   

def destFile(folder, basename, suffix):
    #check if the destination filename already exists
    full = os.path.join(folder, basename + suffix)
    if(False == os.path.exists(full)):
        return full
    else: # generate a new filename, and check if it is available
        counter = 1;
        nf = os.path.join(folder, basename + "_" + repr(counter) + suffix)
        while(True == os.path.exists(nf)):
            counter += 1
            nf = os.path.join(folder, basename + "_" + repr(counter) + suffix)
        return nf

def writePdb(sel, prefix = ""):
    pdbfn = prefix + sel + ".pdb"
    npdbfn = pdbfn.replace(" or ", "");
    npdbfn = npdbfn.replace(" and ", "");
    npdbfn = npdbfn.translate(None, string.whitespace)
    cmd.save(npdbfn, sel)
    return npdbfn

#run crysol in predictive mode for a given selection
def predcrysol(models, prefix=defprefix, param = " "):
    #write all models into a single file  
    selection = " or ".join(models)
    with TemporaryDirectory() as tmpdir:
        #cmd.save(selection, pdbfn)
        pdbfn = writePdb(selection)
        systemCommand(["crysol"] + param.split() + [pdbfn])
        fid = pdbfn.replace(".pdb", "")
        df = tmpdir.move_out_numbered(fid+"00.int", fid, '.int')
    message( ".int file written to " + df)
    openSingleDatFile(datViewer.get(), df)
    return df

cmd.extend("predcrysol", predcrysol)

#run crysol in fit mode
def fitcrysol(SaxsDataFileName, models, prefix = defprefix, param = ""):
    if False == os.path.isfile(SaxsDataFileName):
        message("SAXS .dat file \'"+SaxsDataFileName+"\' not found")
        return
    fileFullPath = os.path.abspath(SaxsDataFileName);

    #write all models into a single file  
    selection = " or ".join(models)
    with TemporaryDirectory() as tmpdir:
        #cmd.save(selection, pdbfn)
        pdbfn = writePdb(selection)
        systemCommand(["crysol"] + param.split() + [pdbfn, fileFullPath])
        fid = pdbfn.replace(".pdb", "")
        df = tmpdir.move_out_numbered(fid+"00.fit", fid, '.fit')
    message( ".fit file written to " + df)
    openSingleDatFile(datViewer.get(), df)
    return df

cmd.extend("fitcrysol", fitcrysol)
        

#run sreflex
def sreflex(SaxsDataFileName, models,
            viewer='sasplot', prefix=defprefix):
    global modelingRuns
    global currentDat
    if False == os.path.isfile(SaxsDataFileName):
        message("SAXS .dat file \'"+SaxsDataFileName+"\' not found")
        return
    fileFullPath = os.path.abspath(SaxsDataFileName);
    cwd = os.getcwd()
    pdbs = []
    for m in models:
        pdbfn = writePdb(m)
        pdbs.append(pdbfn)
    s = ','
    coordsarg = s.join(pdbs)
    df = destFile("", prefix+"_sreflex","")
    systemCommand(["sreflex", "-p", df, fileFullPath, coordsarg])
    message( "sreflex finished." )
    modelingRuns += 1;
    #prepare files for load
    #open report and read entries there
    reportfn = os.path.join(df, "report.txt")
    if False == os.path.isfile(reportfn):
        message("SREFLEX report file \'"+reportfn+"\' not found, something went wrong")
        return
    rf = open(reportfn, 'r')
    currentDat = []
    for line in rf:
        sys.stdout.write(line)
        modelid = line.split()[0]
        if modelid.startswith('rc01') or modelid.startswith('uc01'):
            currentDat.append(df + "/fits/" + modelid + ".fit")
            cmd.load(df + "/models/" + modelid + ".pdb",
                    "sreflex" + repr(modelingRuns) + modelid)
    openDatFile(viewer, currentDat)
    return
cmd.extend("sreflex", sreflex)


def readTransformationMatrixFromPdbRemark(pdbfn):
    #read transformation matrix from output pdb
    #useful for alpraxin and supalm
    rf = open(pdbfn, 'r')
    read = 1
    a=[]
    c=4

    while(1):
        line=rf.readline()
        if re.match("(.*)Transformation(.*)", line):
            while(c):
                a.append(line[39:47])
                a.append(line[47:55])
                a.append(line[55:63])
                a.append(line[63:71])
                line=rf.readline()
                c=c-1;
            break 
    pymolOrder=[1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]
    p=[]; #output in pymol's format
    for i in pymolOrder:
        p.append(float (a[i-1]))

    return p

def alpraxin(sel):
    """run alpraxin and apply transformation matrix"""
    with TemporaryDirectory():
        pdbfn = writePdb(sel, "in_")
        outfn = sel + ".pdb"
        systemCommand(["alpraxin", "-o", outfn, pdbfn])
        tmat = readTransformationMatrixFromPdbRemark(outfn)
        cmd.transform_selection(sel, tmat)            

cmd.extend("alpraxin", alpraxin)

def supalm(template, toalign):
    """run supalm and apply transformation matrix"""
    if(toalign == template):
        message("ERROR Please choose different models for superimposition\n")
        return
    with TemporaryDirectory('supalm'):
        f1 = writePdb(template)
        f2 = writePdb(toalign, "in_")
        outfn = toalign + ".pdb"
        systemCommand(['supalm', '-o', outfn, '--prog2=crysol'] + [f1, f2])
        #reloading file approach, needed for ATSAS 2.7.1
        cmd.delete(toalign)
        cmd.load(outfn) 
        #the following will work once supalm is updated (ATSAS 2.7.2)
        #tmat = readTransformationMatrixFromPdbRemark(outfn)
        #cmd.transform_selection(toalign, tmat)

cmd.extend("supalm", supalm)

def openSingleDatFile(viewer, fn):
    openDatFile(viewer, [fn])

def openDatFile(viewer, fnlst = []):
    #message(repr(fnlst))
    for fn in fnlst:
        if not os.path.isfile(fn):
            message("ERROR Curve file \'" + fn + "\' not found.")
    viewerproc = subprocess.Popen([viewer] + fnlst)

def sasref(SaxsDataFileName, models = [], viewer='sasplot', param = " "):
    global modelingRuns
    if False == os.path.isfile(SaxsDataFileName):
        message("SAXS .dat file \'"+SaxsDataFileName+"\' not found")
        return
    fileFullPath = os.path.abspath(SaxsDataFileName)
    prefix='sasref'
 
    with TemporaryDirectory(prefix) as tmpdir:
        #sasref can not deal with long path/names
        tmpsaxsfn =  os.path.basename(SaxsDataFileName)+".dat"
        tmpdir.copy_in(fileFullPath, tmpsaxsfn);

        #create sasrefJob instance
        sasrefrun = sasrefJob()
        datarray = []
        datarray.append(tmpsaxsfn)

        subunits = []
        for m in models:
            pdbfn = writePdb(m)
            pd = dict()
            pd['filepath'] = pdbfn
            pd['symmetry'] = 'P1'
            pd['multiplicity'] = 1
            subunits.append(pd)
        modelingRuns += 1
        prefix = prefix + repr(modelingRuns)

        # It's possible to do this without creating a command file,
        # by using StringIO, but writing out a file here facilitates
        # debugging if anything goes wrong
        comfn = 'setup_sasrefcv.com'
        with open(comfn, 'w') as commandfile:
            com = sasrefrun.configure(datarray, 'angstrom', subunits,
                                      "", "P1", prefix )
            commandfile.write(com)
        with open(comfn, 'r') as commandfile:
            systemCommand(['sasrefcv'], stdin=commandfile)

        logfn = prefix + ".log"
        cf = prefix + "-1.fit"
        outpdb = prefix + ".pdb"
        cmd.load(outpdb)
        pdf = tmpdir.move_out_numbered(outpdb, prefix, '.pdb')
        df = tmpdir.move_out_numbered(cf, prefix, '.fit')

    message( "sasref finished")
    message( ".fit file written to " + df)
    message( ".pdb file written to " + pdf)
    message( "model loaded as " + prefix)
    openSingleDatFile(viewer, df)
    return df

cmd.extend("sasref", sasref)

def damdisplay(sel, color='white', transparency=0.5):
    #function to make dummy atom models look better
    #cmd.hide(representation="nonbonded", selection = sel);
    cmd.hide(representation="everything", selection = sel);
    cmd.color(color, selection = sel);
    cmd.set("transparency", transparency, selection = sel);
    cmd.show(representation="surface", selection = sel);

cmd.extend("damdisplay", damdisplay);

def updateCurrentDat(newDatFile):
    global currentDat
    currentDat=[]
    currentDat.append(newDatFile)

#class to execute sasrefcv
#stolen from sasflow
class sasrefJob(object):
  #
  # datfilename: experimental data
  # setup: a dictionary of subunits with symmetry and multiplicity info
  # contactsfilename: contact conditions (may be empty)
  #
  def configure(self, datfilenames, unit, subunits, contactsfilename, symmetry, prefix):
#    self.setPriority(task.TASK_PRIORITY_LOW)

    self.curve_con = self.curve_config(datfilenames, unit, prefix)
    self.subunit_con = self.subunit_config(subunits, prefix)
    self.xcorr_con = self.xcorr_config(datfilenames, subunits, prefix)
    self.contacts_con = contactsfilename

    #
    # NOTE: The default contacts keep things together, but severely
    #       increases runtime.
    #
    # if self.contacts_con == "":
    #  self.contacts_con = self.contacts_config(pdbfilenames, prefix)
    #

    ans = "U        ! User mode\n"                  \
          "%s       ! Output file prefix\n"         \
          "%s\n"                                    \
          "%s       ! Master symmetry\n"            \
          "%s       ! Curve config\n"               \
          "         ! Smearing parameters (none)\n" \
          "%s       ! Subunit config\n"             \
          "%s       ! Cross-correlation table\n"    \
          "%s       ! Contact conditions\n"         \
          "Unknown  ! Particle shape\n"

    #self.setStdin(ans % (prefix, prefix, self.symmetry,
    return (ans % (prefix, prefix, symmetry,
                                         self.curve_con,
                                         self.subunit_con,
                                         self.xcorr_con,
                                         self.contacts_con))
    #self.setProgramName("sasrefcv")

    #self.symmetry = symmetry
    #self.log      = prefix + ".log"
    #self.fit      = prefix + "-1.fit"
    #self.pdb      = prefix + ".pdb"


  def curve_config(self, datfilenames, unit, prefix):
    if unit == "angstrom":
      u = 1
    elif unit == "nanometer":
      u = 2

    config = "%d\n" % len(datfilenames)
    for datfilename in datfilenames:
      config += "%s -1.0 P1 %s 1.0 0 1.0 n\n" % (datfilename, u)

    cur_con = prefix + "_cur.con"
    return self.write_config(cur_con, config)

  def subunit_config(self, subunits, prefix):
    m = 0
    config = "%d\n"
    for subunit in subunits:
      n = int(subunit['multiplicity'])
      m += n
      for k in range(n):
        config += "%s y n %s\n" % ( subunit['filepath'],
                                    subunit['symmetry'])

    sub_con = prefix + "_sub.con"
    return self.write_config(sub_con, config % m)

  def xcorr_config(self, datfilenames, subunits, prefix):
    # total number of subunits
    m = 0
    for subunit in subunits:
      m += int(subunit['multiplicity'])

    config = ""
    for datfilename in datfilenames:
      # repeat '0.0' m-times and join, sparated by whitespace
      config += " ".join([ "0.0" ] * m) + "\n"

    xcorr_con = prefix + "_xcorr.con"
    return self.write_config(xcorr_con, config)

  def contacts_config(self, pdbfilenames, prefix):
    #
    # Default contacts conditions; subunits shouldn't be
    # terribly far apart (maxdist, in Angstrom)
    #
    maxdist = 10.0

    config = ""
    for i in range(len(pdbfilenames)):
      config += "dist %f\n" % maxdist
      for j in range(len(pdbfilenames)):
        if i <> j:
          config += "%d 1 0  %d 1 0\n" % (i + 1, j + 1)

    contacts_con = prefix + "_contacts.con"
    return self.write_config(contacts_con, config)

  def write_config(self, confilename, contents):
    with open(confilename, "w") as conf:
        conf.write(contents)
    return confilename
