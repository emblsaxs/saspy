# python lib
'''
SASpy - ATSAS PLUGIN FOR PYMOL

(c) 2015-2017 A.PANJKOVICH FOR ATSAS TEAM AT EMBL-HAMBURG.
'''
import os
import sys
import re
import time
import shutil
import string
import math
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
saspyVersion = "2.8.1"
currentDat = []
modelingRuns = 0
datViewer = Tkinter.StringVar()
cwd = Tkinter.StringVar()
cwd.set(os.getcwd())

from sys import platform
datViewer.set("sasplot") #on linux
if "win32" == platform:
    datViewer.set("sasplotqt")
if "darwin" == platform:
    datViewer.set("/Applications/ATSAS/sasplot.app")

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
        self.sasrefmode    = Tkinter.StringVar()
        self.sasrefmode.set('local')
        self.crysolmode    = Tkinter.StringVar()
        self.crysolmode.set('predict')
        self.prefix        = Tkinter.StringVar()

        self.warnLabel = Tkinter.Label( self.dialog.interior(),
                                    anchor='center',
                                    fg ="red",
        text = 'No models found. Please open/load structures in PyMOL to proceed.')

        description  = "SASpy - ATSAS Plugin for PyMOL\n"
        description += "ATSAS " + saspyVersion + "\n\n"
        description += "European Molecular Biology Laboratory\n"
        description += "Hamburg Outstation, ATSAS Team, 2015-2017.\n"


        w = Tkinter.Label(self.dialog.interior(),
                          text = description,
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
        #self.crymodebut.add('simulate')
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
        alpraxinTab = self.createTab("alpraxin", "Position a structure at the origin such that its principal\ninertia vectors are aligned with the coordinate axis.\nPlease select one or more models.")

        #supalm tab
        supalmTab = self.createTab('supalm', 'Superimposition of models and calculation of\nnormalized spatial discrepancy (NSD).\nPlease select two models.')

        #sasref tab
        sasreftab = self.createTab("sasref", "Quaternary structure modeling against solution scattering data.\nPlease select multiple models (rigid bodies) and a SAXS .dat file.\nRecommendation: execute alpraxin before refinement.")
        saxsfn_ent = Pmw.EntryField(sasreftab,
                                    label_text = 'SAXS .dat file:',
                                    labelpos='ws',
                                    entry_textvariable=self.saxsfn)
        saxsfn_but = Tkinter.Button(sasreftab, text = 'Browse...',
                                    command = self.getSAXSFile)
        saxsfn_ent.grid(sticky='we', row=3, column=0, padx=5, pady=5)
        saxsfn_but.grid(sticky='we', row=3, column=1, padx=5, pady=5)

        self.sasrefmodebut = Pmw.RadioSelect(sasreftab,
                                    buttontype='radiobutton',
                                    labelpos='w',
                                    label_text="Refinement mode:",
                                    command=self.setSasrefMode,
                                    selectmode = 'single')
        self.sasrefmodebut.grid(sticky='we', row=2, column=0, padx=5, pady=2)
        self.sasrefmodebut.add('local')
        self.sasrefmodebut.add('global')
        self.sasrefmodebut.setvalue('local')


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

        self.ATSAS_sanityCheck()

#GUI FUNCTIONS

    def errorWindow(self, title, msg):
        tkMessageBox.showerror(title,
                               "ERROR " + msg, 
                               parent=self.parent)
        return

    def notificationWindow(self, title, msg):
        tkMessageBox.showinfo(title, msg,
                               parent=self.parent)
        return

    def ATSAS_sanityCheck(self):
        msg = checkAtsasBin()
        if "OK" != msg:
            message(msg)
            self.errorWindow("ERROR", msg)
            self.execute("Quit") 
            return
         
        msg = checkAtsasVersion()
        if "OK" != msg:
            message(msg)
            self.errorWindow("ERROR", msg)
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

    def setSasrefMode(self, mode):
        print "Setting sasref mode to "+mode
        self.sasrefmode.set(mode)
        self.sasrefmodebut.setvalue(mode)

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
            t = threading.Thread(target = sasref, name = procType+'_thread', args = (saxsfn, models, self.sasrefmode.get(), viewer))
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
            alpraxin(models)
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
                     'damdisplay':1,
                     'supalm':2,
        #other procedures need a minimum number of selected models
                     'alpraxin':11,
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
        #dots are also removed, as they confuse crysol
        initialList = cmd.get_object_list()
        outputList = list();
        for m in initialList:
            if '_' in m:
                newName = m.translate(None,"_")
                cmd.set_name(m,newName)
                message("WARNING Renaming model \'"+m+ "\' to \'"+ newName+"\'")
                m = newName
            if '.' in m:
                newName = m.translate(None,".")
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
            title='SAXS File', initialdir=cwd,
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
        if 'simulate' == crymode:
            return updateCurrentDat(simulateScattering(selection))
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

#parse crysol log file
def parseCrysolLog (logFileName):
    '''Parse Crysol log file, obtain Chi2 and Rg'''
    #will not parse crysol_summary.txt, but the .log file 
    #created for each individual run

    chi2 = 9999;
    Rg = 9999;

    position = -1
    counter = 0
    with open(logFileName, 'r') as rf:
        for line in rf:
            counter += 1
            if re.match("(.*)Fitting parameters(.*)", line):
                print "line number: " + repr(counter)
                position = counter + 2
            if counter == position:
                if line[66:73] != "*******":
                    chi2 = float(line[66:73])
            if re.match("(.*)Rg from the slope of net intensity(.*)", line):
                Rg = float(line[59:65])
    rf.close()
    return {'chi2':chi2, 'Rg':Rg}

def simulateScattering(models, prefix=defprefix, param = " "):
    '''Use CRYSOL and ADDERRORS to simulate scattering''' 
    #write all models into a single file  
    selection = " or ".join(models)
    Rg = -9999
    df = 'unknown'
    with TemporaryDirectory() as tmpdir:
        pdbfn = writePdb(selection)
        systemCommand(["crysol", "-ns", "800"] + param.split() + [pdbfn])
        fid = pdbfn.replace(".pdb", "")
        result = parseCrysolLog(fid+"00.log")
        Rg = result['Rg']
        tmpint = fid + "00.int"
        tmpout = fid + ".dat"
        systemCommand(["adderrors", tmpint, "-o", tmpout])
        df = tmpdir.move_out_numbered(tmpout, fid, '.dat')

    message("CRYSOL Theoretical Rg = " + repr(Rg))
    message( ".dat file written to " + df)
    openSingleDatFile(datViewer.get(), df)
    return df

cmd.extend("simulateScattering", simulateScattering)

#run crysol in predictive mode for a given selection
def predcrysol(models, prefix=defprefix, param = " "):
    #write all models into a single file  
    selection = " or ".join(models)
    Rg = -9999
    df = 'unknown'
    with TemporaryDirectory() as tmpdir:
        #cmd.save(selection, pdbfn)
        pdbfn = writePdb(selection)
        systemCommand(["crysol"] + param.split() + [pdbfn])
        fid = pdbfn.replace(".pdb", "")
        result = parseCrysolLog(fid+"00.log")
        Rg = result['Rg']
        df = tmpdir.move_out_numbered(fid+"00.int", fid, '.int')

    message("CRYSOL Theoretical Rg = " + repr(Rg))
    message( ".int file written to " + df)
    openSingleDatFile(datViewer.get(), df)
    return df

cmd.extend("predcrysol", predcrysol)

def checkAtsasBin():
    '''Check if ATSAS binaries are available on PATH''' 
    try: 
        status = subprocess.check_output(["crysol", "-v"])
    except OSError:
        return "\nATSAS executables not found in PATH.\n\nPlease install ATSAS.\n"
    return "OK"

def checkAtsasVersion():
    '''Check if the installed ATSAS version matches what SASpy expects'''
    try:
        output = subprocess.check_output(["crysol", "-v"],
                                    stderr=subprocess.STDOUT)
    except OSError:
        return "\nATSAS executables not found in PATH.\n\nPlease install ATSAS."
    except subprocess.CalledProcessError as exc:
        msg = "\nError running `crysol -v` to check the ATSAS version.\n\n"
        msg += "Check that ATSAS is properly installed.\n\n"
        msg += "Output of the failed command:\n"
        msg += exc.output
        return msg

    versionMatch = re.match(r"crysol, ATSAS v?(\S+)\s+", output)
    if versionMatch is None:
        msg = "\nCould not parse the ATSAS version from `crysol -v` output:\n\n"
        msg += output
        return msg

    installedAtsasVersion = versionMatch.group(1)
    msg = "Found binaries for ATSAS version " + installedAtsasVersion
    message(msg)
    if installedAtsasVersion != saspyVersion:
        msg = "\nSASpy and ATSAS versions do not match.\n\n"
        msg += "Currently installed:\n"
        msg += "SASpy " + saspyVersion + "\n"
        msg += "ATSAS " + installedAtsasVersion + "\n\n"
        msg += "SASpy and ATSAS versions should be equivalent.\nPlease update accordingly."
        return msg
    return "OK"

cmd.extend("checkAtsasVersion", checkAtsasVersion)
 
#run crysol in fit mode
def fitcrysol(SaxsDataFileName, models, prefix = defprefix, param = ""):
    if False == os.path.isfile(SaxsDataFileName):
        message("SAXS .dat file \'"+SaxsDataFileName+"\' not found")
        return
    fileFullPath = os.path.abspath(SaxsDataFileName);

    Rg = -9999
    chi2 = -9999
    
    #write all models into a single file  
    selection = " or ".join(models)
    with TemporaryDirectory() as tmpdir:
        #cmd.save(selection, pdbfn)
        pdbfn = writePdb(selection)
        systemCommand(["crysol"] + param.split() + [pdbfn, fileFullPath])
        fid = pdbfn.replace(".pdb", "")
        result = parseCrysolLog(fid+"00.log")
        Rg = result['Rg']
        chi2 = result['chi2']
        df = tmpdir.move_out_numbered(fid+"00.fit", fid, '.fit')

        #if there is more than one model, we are evaluating a complex
        #in this case we should provide the coordinates of the complex
        #to the user.
        if 1 < len(models):
            pdbfn=tmpdir.move_out_numbered(pdbfn, fid, '.pdb')
            message( ".pdb file written to " + pdbfn)

    message( ".fit file written to " + df)
    message("CRYSOL Theoretical Rg = " + repr(Rg))
    message("CRYSOL Chi-square = " + repr(chi2))
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
    currentDat = []
    with open(reportfn, 'r') as rf:
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

def readNSDFromSupalmPdb(pdbfn):
    """parse NSD value from SUPALM output PDB file
    useful for Windows users without console access"""
    nsd = 9999;
    with open(pdbfn, 'r') as rf:
        for line in rf:
            if re.match("(.*)Final distance(.*)", line):
                nsd = float(line[38:49])
                break
    return nsd

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
                a.append(line[39:51])
                a.append(line[51:63])
                a.append(line[63:75])
                a.append(line[75:87])
                line=rf.readline()
                c=c-1;
            break 
    rf.close()
    pymolOrder=[1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16]
    p=[]; #output in pymol's format
    for i in pymolOrder:
        p.append(float (a[i-1]))
    return p

def alpraxin(models):
    """run alpraxin and apply transformation matrix"""
    sel = " or ".join(models)
    with TemporaryDirectory():
        pdbfn = writePdb(sel, "in_")
        outfn = sel + ".pdb"
        outfn = outfn.replace(" ", "")
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
        sargs = ['supalm', '-o', outfn]
        sargs.append('--prog2=crysol')
        sargs.append('--enantiomorphs=N')
        sargs.append(f1)
        sargs.append(f2)
        systemCommand(sargs)
        tmat = readTransformationMatrixFromPdbRemark(outfn)
        nsd = readNSDFromSupalmPdb(outfn)
        cmd.transform_selection(toalign, tmat)
        message("SUPALM NSD = " + repr(nsd))
cmd.extend("supalm", supalm)

def openSingleDatFile(viewer, fn):
    openDatFile(viewer, [fn])

def openDatFile(viewer, fnlst = []):
    #message(repr(fnlst))
    for fn in fnlst:
        if not os.path.isfile(fn):
            message("ERROR Curve file \'" + fn + "\' not found.")
    if('darwin' == platform):
        viewerproc = subprocess.Popen(["open"] + ["-a"] + [viewer] + fnlst)       
    else:
        viewerproc = subprocess.Popen([viewer] + fnlst) 
 
def damdisplay(sel, color='white', transparency=0.5):
    '''set visualization of dummy atom models'''
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

def allToRefAlign(ref):
    for i in cmd.get_object_list():
        cmd.align(i, ref);

cmd.extend("allToRefAlign", allToRefAlign);

def parseEulerAngles(filename):
    '''Parse Euler angles and translations 
    per subunit from SASREF output PDB'''

    #temporal storage of values, one set per subunit
    collection = []
    move = []

    with open(filename, 'r') as rf:
        for line in rf:
            if re.match("REMARK Old center positioned at(.*)", line):
                move.append(float(line[32:39]))
                move.append(float(line[40:47]))
                move.append(float(line[48:55]))
            if re.match("REMARK Rotated by Euler angles(.*)", line):
                move.append(float(line[32:39]))
                move.append(float(line[40:47]))
                move.append(float(line[48:55]))
            if re.match("REMARK New center positioned at(.*)", line):
                move.append(float(line[32:39]))
                move.append(float(line[40:47]))
                move.append(float(line[48:55]))
            if re.match("TER", line):
                collection.append(move)
                move = []
    return collection

def anglesToTTTMat(movement):
    '''Return a PyMOL transformation matrix from a set of 
    SASREF translation vectors and Euler angles'''
    alpha = math.radians(movement[3])
    beta  = math.radians(movement[4])
    gamma = math.radians(movement[5])
    
    sinalp = - math.sin(alpha)
    cosalp =   math.cos(alpha)
    sinbet = - math.sin(beta)
    cosbet =   math.cos(beta)
    singam = - math.sin(gamma)
    cosgam =   math.cos(gamma)

    #transforming to PyMOL TTT as explained in
    #http://www.pymolwiki.org/index.php/Transform_selection

    output = []
    output.append(cosalp*cosbet*cosgam - sinalp*singam)
    output.append(-cosalp*cosbet*singam - sinalp*cosgam)
    output.append(cosalp*sinbet)
    output.append(movement[6])
    output.append(sinalp*cosbet*cosgam + cosalp*singam)
    output.append(-sinalp*cosbet*singam + cosalp*cosgam)
    output.append(sinalp*sinbet)
    output.append(movement[7])
    output.append(-sinbet*cosgam)
    output.append(sinbet*singam)
    output.append(cosbet)
    output.append(movement[8])
    output.append(-movement[0])
    output.append(-movement[1])
    output.append(-movement[2])
    output.append(1.0)
    return output

def sasref(SaxsDataFileName, models = [], mode = 'local', viewer='sasplot'):
    '''Execute SASREF and apply obtained transformations to subunits'''

    global modelingRuns

    #parameter configuration
    #local refinment (to reproduce MASSHA behaviour):
    confp = {'spst':'1.0', # spatial step, SASREF default is 5.0 Angstrom
             'anst':'5.0', #angular step, SASREF default is 20 degrees
             'init':'1.0', # Initial annealing temperature (def= 10) 
             'sche':'0.9', # Annealing schedule factor
             'iter':'500', # Max # of iterations at each T (def = 10000)
             'maxs':' 50', # Max # of successes at each T (def = 1000)
             'mins':' 10', # Min # of successes to continue (def = 100)
             'maxa':' 10', # Max # of annealing steps (def = 100)
             'msol':'  1', # Max # of solutions to store (def = 1)
             'shft':' '}   #xyz init shift, blank for local
                           #0.0 for global search

    if 'global' == mode: #settings for default global search
        confp['spst'] =  '5.0'
        confp['anst'] = '20.0'
        confp['init'] = '10.0'
        confp['sche'] =  '0.9'
        confp['iter'] = '10000'
        confp['maxs'] = '1000'
        confp['mins'] = '100'
        confp['maxa'] = '100'
        confp['msol'] = '1'
        confp['shft'] = '0.0'

    if False == os.path.isfile(SaxsDataFileName):
        message("SAXS .dat file \'"+SaxsDataFileName+"\' not found")
        return 113
    fileFullPath = os.path.abspath(SaxsDataFileName)

    prefix = 'sasref'
    modelingRuns += 1
    prefix = prefix + repr(modelingRuns)

    numberOfSubunits = len(models);

    with TemporaryDirectory(prefix) as tmpdir:
        #sasref can not deal with long path/names
        tmpsaxsfn =  os.path.basename(SaxsDataFileName)
        tmpdir.copy_in(fileFullPath, tmpsaxsfn);

        sc = ""
        sc += "Expert  ! Configuration mode\n"
        sc += prefix + " ! logfilename\n"
        sc += prefix + " ! projectDescription\n"
        sc += "        ! initRandomSeed\n"
        sc += repr(numberOfSubunits) + "     ! totalNumberOfSubunits\n"
        sc += "P1      ! Symmetry\n"
        sc += "1       ! totalNumberOfScatteringCurves\n"
        sc += "N       ! KratkyGeometry\n"
        sc += "1,2     ! Input first & last subunits in 1-st construct\n"
        sc += tmpsaxsfn+ " ! Enter file name, 1-st experimental data\n"
        sc += "        ! Angular units input file\n"
        sc += "1.0     ! Fitting range in fractions of Smax\n"

        count = 1
        for m in models:
            pdbtmpfn = writePdb(m)
            pdbtmpfn = os.path.basename(pdbtmpfn)
            fid = os.path.splitext(pdbtmpfn)[0]
            #print "abs: " + m + " and temp: "+pdbtmpfn + "while id: "+fid;
            #tmpdir.copy_in(m, pdbtmpfn);
            #compute amplitudes
            systemCommand(["crysol"] + ["-p"] + [fid] + [pdbtmpfn])
            print "computed alm for : "+ fid +"\n"
            sc += fid+".alm"+"! subunit " +repr(count)+" amplitudes\n"
            sc += "0.0     ! Initial rotation by alpha\n"
            sc += "0.0     ! Initial rotation by beta\n"
            sc += "0.0     ! Initial rotation by gamma\n"
            sc += confp['shft']+"       ! Initial shift along X\n"
            sc += confp['shft']+"       ! Initial shift along Y\n"
            sc += confp['shft']+"       ! Initial shift along Z\n"
            sc += "N       ! Movements limitations of subunit:  N/F/X/Y/Z/D\n"
            sc += confp['spst'] + " ! Spatial step in Angstrom\n"
            sc += confp['anst'] + " ! Angular step in degrees\n"


        #continue with rest of commands
        sc += "        ! Cross penalty weight\n" 
        sc += "        ! Disconnectivity penalty weight\n"
        sc += "        ! Docking penalty weight\n"
        sc += "        ! File name, contacts conditions, CR for none\n"
        sc += "U       ! Expected particle shape: <P>rolate, <O>blate or <U>nknown\n"
        sc += "        ! Shift penalty weight \n"
        sc += confp['init'] + "     ! Initial annealing temperature (def= 10)   \n"
        sc += confp['sche'] + "     ! Annealing schedule factor\n"
        sc += confp['iter'] +  "    ! Max # of iterations at each T (def = 10000)\n"
        sc += confp['maxs'] + "     ! Max # of successes at each T (def = 1000)\n"
        sc += confp['mins'] +"      ! Min # of successes to continue (def = 100)\n"
        sc += confp['maxa'] +"      ! Max # of annealing steps (def = 100)\n"
        sc += confp['msol']+"       ! Max # of solutions to store (def = 1)\n"
        comfn = 'setup_sasref.com'
 
        with open(comfn, 'w') as commandfile:
            commandfile.write(sc)
        with open(comfn, 'r') as commandfile:
            systemCommand(['sasref'], stdin=commandfile)

        outpdb = prefix + ".pdb" 
        #read and apply movements 
        moves = parseEulerAngles(outpdb);
        idx = 0
        for mov in moves:
            tmat = anglesToTTTMat(mov)
            cmd.transform_selection(models[idx], tmat)            
            idx = idx+1

        outpdbfn = tmpdir.move_out_numbered(prefix + ".pdb", prefix, '.pdb')
        message( ".pdb file written to " + outpdbfn)
        cf = tmpdir.move_out_numbered(prefix + "-1.fit", prefix, '.fit')
        message( ".fit file written to " + cf)
        openSingleDatFile(viewer, cf)
    return 



