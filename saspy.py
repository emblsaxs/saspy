# python lib
'''
SASpy - ATSAS PLUGIN FOR PYMOL

(c) 2015 A.PANJKOVICH FOR ATSAS TEAM AT EMBL-HAMBURG.
'''
import os
import sys
import time
import shutil
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
datmode = Tkinter.StringVar()
datmode.set('immediately')

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
                                 '3. Execute', 
                                 '4. View SAXS curve'),
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
                          text = '\nSASpy - ATSAS Plugin for PyMOL\nVersion 1.0 - ATSAS 2.7.0\n\nEuropean Molecular Biology Laboratory\nHamburg Outstation, ATSAS Team, 2015.\n',
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
        crysolTab = self.createTab("crysol", '''Prediction of theoretical intensities and
optionally fit to experimental SAXS data.
Please select one model (and a SAXS .dat file for fit mode).'''
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
        self.notebook.setnaturalsize(pageNames=fn)
            
        #alpraxin tab
        alpraxinTab = self.createTab("alpraxin", "Position a structure such that its principal intertia\nvectors are aligned with the coordinate axis.\nPlease select one model.")
    
        #supcomb tab
        supcombTab = self.createTab('supcomb', 'Superimposition of models and calculation of\nnormalized spatial discrepancy (NSD).\nPlease select two models.')
        #join models tab
        joinTab = self.createTab("createComplex", "Create a complex from selected models.\nPlease provide a name and select models to join.")
        self.newModelName = Tkinter.StringVar();
        newModelNameEntry = Pmw.EntryField(joinTab,
                                    label_text = 'New complex name:', 
                                    labelpos='ws',
                                    entry_textvariable=self.newModelName)
        newModelNameEntry.grid(sticky='we', row=3, column=0, padx=5, pady=5)

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


        # config tab
        configTab = self.createTab("configure", "Settings available to configure SASpy:")
        svi_ent = Pmw.EntryField(configTab,
                                    label_text = 'SAXS viewer:', 
                                    labelpos='ws',
                                    entry_textvariable = datViewer)
        svi_but = Tkinter.Button(configTab, text = 'Select...',
                                    command = self.getSAXSViewer)
        svi_ent.grid(sticky='we', row=3, column=0, padx=5, pady=5)
        svi_but.grid(sticky='we', row=3, column=1, padx=5, pady=5)

        self.datmodebut = Pmw.RadioSelect(configTab,
                                    buttontype='radiobutton',
                                    labelpos='w',
                                    label_text="Open results:",
                                    command=self.setDatMode,
                                    selectmode = 'single')
        self.datmodebut.grid(sticky='we', row=2, column=0, padx=5, pady=2)
        self.datmodebut.add('immediately')
        self.datmodebut.add('on request')
        self.datmodebut.setvalue('immediately')

        #Model selection
        self.modsW = self.createModelSelectionWidget()
        self.modsW.pack(expand=1, fill='both', padx=10, pady=5)        
        self.refreshModelSelectionWidget()


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
        
       
    def submitSingleModelJob(self, procType, selection):
        global currentDat
        if "crysol" == procType:
            crymode = self.crysolmode.get()
            print 'crymode is: ' + crymode
            if 'predict' == crymode:
                updateCurrentDat(predcrysol(selection))
            elif 'fit' == crymode:
                saxsfn = self.saxsfn.get()
                if False == os.path.isfile(saxsfn):
                    self.errorWindow("FILE NOT FOUND",
                                 "SAXS file \'"+saxsfn+"\' NOT FOUND.");
                    return
                updateCurrentDat(fitcrysol(saxsfn, selection))

        elif "alpraxin" == procType:
            alpraxin(selection)
            return
        else:
            self.errorWindow("Not enough models selected",
                "Please select more models for routine \'"+procType+"\'")
            return
        return;

    def submitJobAsThread(self, procType, models= [], saxsfn = ""):
        viewer = datViewer.get()
        #check how many threads from this plugin are running
        running = self.checkAtsasThreads()
        if running == self.maxThreads:
            self.errorWindow("Max threads exceeded",
                "A process is already running, please wait for it to complete.")
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
        
    def submitMultiModelJob(self, procType, models = []):
        if "createComplex" == procType:
            self.createComplex(self.newModelName.get(), models)
        if "sasref" == procType or "sreflex" == procType:
            saxsfn = self.saxsfn.get()
            if False == os.path.isfile(saxsfn):
                self.errorWindow("FILE NOT FOUND",
                                 "SAXS file \'"+saxsfn+"\' NOT FOUND.");
                return
            if "sreflex" == procType:
                self.submitJobAsThread(procType, models, saxsfn)
            if "sasref" == procType:
                self.submitJobAsThread(procType, models, saxsfn)
        if "supcomb" == procType:
            supcomb(models[0], models[1])
        return 

    def prepareJobForSubmit(self):
        procType = self.procedure
        nModels = self.countSelectedModels()

        if 'configure' == procType:
            return
        expn_dict = {
                     'crysol':1,
                     'alpraxin':1,
                     'supcomb':2,
                     'createComplex':99,
                     'sreflex':99,
                     'sasref':99
                     }

        expected_n = expn_dict[procType]
        if expected_n != 99:
            if expected_n != nModels:
                self.errorWindow("Wrong number of models selected", 
                "You selected "+ repr(nModels) +" models, but \'" + procType+"\' expects " + repr(expected_n) + " model(s).\n")
                return
            
       
        if 0 == nModels:
            self.errorWindow("No model selected", "Please select models")
            return
        if 1 == nModels:
            if 'sreflex' == procType:
                self.submitMultiModelJob(procType, self.modsW.getcurselection())
            else:
                selection = self.modsW.getcurselection()[0]
                self.submitSingleModelJob(procType, selection)
        if 1 < nModels:
            self.submitMultiModelJob(procType, self.modsW.getcurselection())
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
        if 0 == len(currentDat):
            tkMessageBox.showerror('No curve yet',
                                   'No SAXS intensities have been calculated yet',
                                    parent=self.parent)
               
        else:
            openDatFile(datViewer.get(), currentDat) 
        return

    def createComplex(self, newName, models = []):
        mols = self.getListOfModels()
        if newName in mols:
            self.errorWindow("Model name in use", "Please choose a different name for the new model")
            return
        #create working dir
        with TemporaryDirectory('createComplex'):
            #write every model to file, then read back in as single model
            with open(newName+".pdb", 'w') as of:
                for p in models:
                    pdbfn = writePdb(p)
                    with open(pdbfn, 'r') as rf:
                        for line in rf:
                            if not line.startswith("END"): 
                                of.write(line)
            cmd.load(newName+".pdb")

    def getSAXSViewer(self):
        global datViewer
        file_name = tkFileDialog.askopenfilename(
            title='SAXS viewer', initialdir='',
            parent=self.parent)
        datViewer.set(file_name)
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
        self.procedure = pagename
        return
 
    def execute(self, cmd):
        """ Run the cmd represented by the button clicked by user.
        """        
        if cmd == 'OK':
            print 'is everything OK?'

        elif cmd == 'Refresh model list':
            self.refreshModelSelectionWidget()

#        elif cmd == 'Debug':

        elif cmd == '3. Execute':
            self.prepareJobForSubmit()
      
        elif cmd == '4. View SAXS curve':
            self.openCurrentDatFile()            

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
    cmd.save(pdbfn, sel)
    return pdbfn

#run crysol in predictive mode for a given selection
def predcrysol(sel, prefix=defprefix, param = " "):
    with TemporaryDirectory() as tmpdir:
        pdbfn = writePdb(sel)
        cf = sel + "00.int"
        systemCommand(["crysol"] + param.split() + [pdbfn])
        df = tmpdir.move_out_numbered(cf, sel, ".int")
    message( ".int file written to " + df)
    if 'immediately' == datmode.get():
        openSingleDatFile(datViewer.get(), df)
    return df

cmd.extend("predcrysol", predcrysol)

#run crysol in fit mode
def fitcrysol(SaxsDataFileName, sel, prefix = defprefix, param = ""):
    if False == os.path.isfile(SaxsDataFileName):
        message("SAXS .dat file \'"+SaxsDataFileName+"\' not found")
        return
    fileFullPath = os.path.abspath(SaxsDataFileName);
    with TemporaryDirectory() as tmpdir:
        pdbfn = writePdb(sel)
        systemCommand(["crysol"] + param.split() + [pdbfn, fileFullPath])
        cf = sel + "00.fit"
        df = tmpdir.move_out_numbered(cf, sel, '.fit')
    message( ".fit file written to " + df)
    if 'immediately' == datmode.get():
        openSingleDatFile(datViewer.get(), df)
    return df

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


def alpraxin(sel):
    """run alpraxin and load new orientation"""
    with TemporaryDirectory():
        pdbfn = writePdb(sel, "in_")
        outfn = sel + ".pdb"
        systemCommand(["alpraxin", "-o", outfn, pdbfn])
        cmd.delete(sel)
        cmd.load(outfn)

cmd.extend("alpraxin", alpraxin)


def openSingleDatFile(viewer, fn):
    openDatFile(viewer, [fn])


def openDatFile(viewer, fnlst = []):
    for fn in fnlst:
        if not os.path.isfile(fn):
            message("ERROR Curve file \'" + fn + "\' not found.")
    viewerproc = subprocess.Popen([viewer] + fnlst)


def supcomb(template, toalign, output="empty", prefix=defprefix, param=""):
    """run supcomb and load new orientation"""
    if(toalign == template):
        message("ERROR Please choose different models for superimposition\n")
        return
    with TemporaryDirectory():
        f1 = writePdb(template)
        f2 = writePdb(toalign, "in_")
        outfn = toalign + ".pdb"
        if "" == param:
            param = ['-m', 'F', '--superposition=B']
        else:
            param = param.split()
        systemCommand(['supcomb', '-o', outfn] + param + [f1, f2])
        cmd.delete(toalign)
        cmd.load(outfn)

cmd.extend("supcomb", supcomb)


def sasref(SaxsDataFileName, models = [], viewer='sasplot', param = " "):
    global modelingRuns
    if False == os.path.isfile(SaxsDataFileName):
        message("SAXS .dat file \'"+SaxsDataFileName+"\' not found")
        return
    fileFullPath = os.path.abspath(SaxsDataFileName)
    prefix='sasref'
    with TemporaryDirectory(prefix) as tmpdir:
        #create sasrefJob instance
        sasrefrun = sasrefJob() 
        datarray = []
        datarray.append(SaxsDataFileName)

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
        df = tmpdir.move_out_numbered(cf, prefix, '.fit')
    message( "sasref finished")
    message( ".fit file written to " + df)
    message( "model loaded as " + prefix)
    openSingleDatFile(viewer, df)
    return df

cmd.extend("sasref", sasref)

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
