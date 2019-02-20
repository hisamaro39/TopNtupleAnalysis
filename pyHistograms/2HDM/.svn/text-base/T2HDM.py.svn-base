#!/usr/bin/python

import numpy as NUMPY
import subprocess
from contextlib import contextmanager
import os
import sys
from pipes import quote
from pprint import pprint
import imp
import ROOT
from ROOT import *

class t2HDM:
   """The 2HDM definitions"""
   def __init__(self, nameX="A", mX=500, typeX=2, sba=1, mintanb=0.3, maxtanb=0.3):
      self.nameX   = nameX
      self.mX      = mX
      self.typeX   = typeX
      self.sba     = sba
      self.mintanb = mintanb
      self.maxtanb = maxtanb
      self.cuts  = "m"+nameX+"=="+str(mX)+" && sba=="+str(sba)+" && tanb>="+str(mintanb)+" && tanb<="+str(maxtanb)
      self.cuts += " && TMath::ATan(tanb)>0. && TMath::ATan(tanb)<TMath::Pi()/2."
      self.cuts += " && TMath::Abs(cba)<=1."
      self.cuts += " && type=="+str(typeX)
      # self.cuts += " && (status&7)==0" 
      self.cuts += " && (status&3)==0"
      self.AlphaS  = 0.1 #lhapdf.mkAlphaS("NNPDF30_nlo_as_0118")
      self.workpath = "/afs/cern.ch/user/h/hod/data/lib2HDM/"
      self.mgpath   = "/afs/cern.ch/user/h/hod/MadGraph/MG5_aMC_v2_4_3/"
      self.diagramsSM = [ "P0_gg_ttx",
                          "P0_uux_ttx",
                          "P1_uux_ttxg",
                          "P1_gux_ttxux",
                          "P1_gu_ttxu",
                          "P1_gg_ttxg",
                          "P2_gu_ttxgu",
                          "P2_uxux_ttxuxux",
                          "P2_uxcx_ttxuxcx",
                          "P2_uux_ttxuux",
                          "P2_uux_ttxgg",
                          "P2_uux_ttxccx",
                          "P2_uu_ttxuu",
                          "P2_ucx_ttxucx",
                          "P2_uc_ttxuc",
                          "P2_gux_ttxgux",
                          "P2_gg_ttxgg",
                          "P2_gg_ttxuux"
                        ]
      self.diagramsXX = self.diagramsSM +  [ "P0_bbx_ttx",
                                             "P1_bbx_ttxg",
                                             "P1_gbx_ttxbx",
                                             "P1_gb_ttxb",
                                             "P2_gb_ttxgb",
                                             "P2_bxbx_ttxbxbx",
                                             "P2_uxbx_ttxuxbx",
                                             "P2_bbx_ttxbbx",
                                             "P2_bbx_ttxgg",
                                             "P2_uux_ttxbbx",
                                             "P2_bbx_ttxuux",
                                             "P2_bb_ttxbb",
                                             "P2_ubx_ttxubx",
                                             "P2_ub_ttxub",
                                             "P2_gbx_ttxgbx",
                                             "P2_gg_ttxbbx",
                                             #"P2_uxb_ttxuxb"
                                            ]

      self.ptest4 = [[148.7495, 0.0, 0.0, 148.7495],
                     [346.258375, 0.0, 0.0, -346.258375],
                     [285.4675, -99.5479765625, 60.29248046875, -197.2335625],
                     [209.539609375, 99.5479765625, -60.29248046875, -0.2753184814453125]]
      self.ptest5 = [[264.3046875, 0.0, 0.0, 264.3046875],
                     [328.02553125, 0.0, 0.0, -328.02553125],
                     [242.8330625, 18.229767578125, 58.28096484375, 157.198890625],
                     [199.76946875, -59.23303515625, -29.940345703125, -79.73228125],
                     [149.727625, 41.00326953125, -28.34062109375, -141.18746875]]
      self.ptest6 = [[596.8451875, 0.0, 0.0, 596.8451875],
                     [245.78221875, 0.0, 0.0, -245.78221875],
                     [433.40378125, -140.238125, -88.751625, 361.48684375],
                     [238.4015625, 154.838421875, 24.286923828125, 54.7362890625],
                     [56.67294140625, -6.11514990234375, 34.39092578125, 44.62837109375],
                     [114.1488359375, -8.485150390625, 30.0737734375, -109.7885546875]]
   def constrainWidth(wmin,wmax):
      self.cuts += " && (width_"+nameX+"/m"+nameX+">"+str(wmin)+" && width_"+nameX+"/m"+nameX+"<"+str(wmax)+")"

   def getProcName(self,proc):
      procname = proc
      procname = procname.replace("P0_","")
      procname = procname.replace("P1_","")
      procname = procname.replace("P2_","")
      procname = procname.replace("_","")
      return procname

##########################
### load the default model
model = t2HDM() ##########
##########################



class fermions:
   def __init__(self):
      self.u = 1
      self.d = 2
      self.s = 3
      self.c = 4
      self.b = 5
      self.t = 6
      self.e = 11
      self.mu = 13
      self.tau = 15


MC = 1.42
MB = 4.7
MT = 172.5
MM = 0.10566
MTAU = 1.777

### the data structure
parameters = []

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def invert_momenta(p):
   """ fortran/C-python do not order table in the same order"""
   new_p = []
   for i in range(len(p[0])):  new_p.append([0]*len(p))
   for i, onep in enumerate(p):
      for j, x in enumerate(onep):
         new_p[j][i] = x
   return new_p


def couplings(typeX,nameX,fermion,tanb,sba,cba):
   Fermions = fermions()
   g = 0.
   if(typeX==1):
      if(nameX=="h"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = sba+cba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = sba+cba/tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = sba+cba/tanb
      if(nameX=="H"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = cba-sba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = cba-sba/tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = cba-sba/tanb
      if(nameX=="A"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = +1./tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = -1./tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = -1./tanb
   if(typeX==2):
      if(nameX=="h"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = sba+cba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = sba-cba*tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = sba-cba*tanb
      if(nameX=="H"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = cba-sba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = cba+sba*tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = cba+sba*tanb
      if(nameX=="A"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = +1./tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = tanb
   return g


def setParameters(nameX,mX,cuts="",typeX=2,sba=1):
   f = TFile("/afs/cern.ch/user/h/hod/data/thdm_grid_v166.root","READ")
   t = f.Get("thdm")
   b_tb  = NUMPY.zeros(1, dtype=float)
   b_sba = NUMPY.zeros(1, dtype=float)
   b_cba = NUMPY.zeros(1, dtype=float)
   b_wA  = NUMPY.zeros(1, dtype=float)
   b_wH  = NUMPY.zeros(1, dtype=float)
   b_mA  = NUMPY.zeros(1, dtype=float)
   b_mH  = NUMPY.zeros(1, dtype=float)
   t.SetBranchAddress("tanb",b_tb);
   t.SetBranchAddress("sba",b_sba);
   t.SetBranchAddress("cba",b_cba);
   t.SetBranchAddress("width_A",b_wA);
   t.SetBranchAddress("width_H",b_wH);
   t.SetBranchAddress("mA",b_mA);
   t.SetBranchAddress("mH",b_mH);

   print "Before cuts = "+str(t.GetEntries())
   print "cuts: "+cuts
   t.Draw(">>elist",cuts,"entrylist");
   elist = gDirectory.Get("elist");
   t.SetEntryList(elist);
   Nlist = elist.GetN()
   print "After cuts = "+str(Nlist)

   Fermions = fermions()
   global parameters
   for i in range(0,Nlist):
      n = elist.Next()
      t.GetEntry(n)
      tanb = b_tb[0]
      sba  = b_sba[0]
      cba  = b_cba[0]
      wA   = b_wA[0]
      wH   = b_wH[0]
      mA   = b_mA[0]
      mH   = b_mH[0]
      YMT   = couplings(typeX,nameX,Fermions.t,tanb,sba,cba)*MT
      YMB   = couplings(typeX,nameX,Fermions.b,tanb,sba,cba)*MB
      YMC   = couplings(typeX,nameX,Fermions.c,tanb,sba,cba)*MC
      YMM   = couplings(typeX,nameX,Fermions.mu,tanb,sba,cba)*MM
      YMTAU = couplings(typeX,nameX,Fermions.tau,tanb,sba,cba)*MTAU
      print "["+str(i)+"] tanb="+'%.6f' % tanb+" sba="+'%.6f' % sba+" cba="+'%.6f' % cba+" wA="+'%.6f' % wA+" wH="+'%.6f' % wH+" YMT="+'%.6f' % YMT+" YMB="+'%.6f' % YMB+" YMC="+'%.6f' % YMC+" YMM="+'%.6f' % YMM+" YMTAU="+'%.6f' % YMTAU
      adict = {}
      adict = {'tanb':tanb, 'sba':sba, 'cba':cba, 'mA':mA, 'wA':wA, 'mH':mH, 'wH':wH, 'YMT':YMT, 'YMB':YMB, 'YMC':YMC, 'YMM':YMM, 'YMTAU':YMTAU}
      parameters.append(adict.copy())
   print "N="+str(len(parameters))+" parameters set !"


def getItanb(tanb):
   for i in xrange(len(parameters)):
      if(parameters[i]['tanb']==tanb): return i
   print "Cannot find index of tanb=%g in parameters list of disctionaries" % tanb
   return -1


def compileSM(nameX,mX,sba):	
   X       = "matrix/"+nameX+"/"+str(mX)+"/"+str(sba)+"/"
   command = "./bin/mg5_aMC run/proc.SM.dat"
   procdirbase = "pp-SM-ttjets/SubProcesses/"
   procdirs = []
   for proc in model.diagramsSM:
      procdirs.append(procdirbase+proc+"/")
   matxdir = model.workpath+X
   libdir  = matxdir+"SM/"

   global parameters

   p = subprocess.Popen("rm -rf "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   p = subprocess.Popen("mkdir -p "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   for proc in model.diagramsSM:
      p = subprocess.Popen("mkdir -p "+libdir+"/"+proc, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

   # enter the directory like this:
   with cd(model.mgpath):
      # make sure the old directory is removed
      p = subprocess.Popen("rm -rf "+procdirbase, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      # execute the generation of the process
      p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      p = subprocess.Popen("/bin/cp -f "+procdirbase+"makefile "+procdirbase+"makefile_ORIG", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      notMade = True
      for procdir in procdirs:
         # go to make the library
         with cd(procdir):
            proc = procdir
            proc = proc.replace(procdirbase,"")
            proc = proc.replace("/","")
            procname = model.getProcName(proc)
            if(notMade):
               p = subprocess.Popen("make", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
               out, err = p.communicate()
               notMade = False
         
            ### cahnge the makefile, make and copy
            p = subprocess.Popen("/bin/cp -f ../makefile_ORIG ../makefile", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen('sed -i -e "s/MENUM)py/MENUM)SM'+procname+'py/g" '+model.mgpath+procdir+'makefile', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("make matrix2SM"+procname+"py.so", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("cp matrix2SM"+procname+"py.so "+libdir+proc+"/", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("cp ../../Cards/*.dat "+libdir+proc+"/", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("rm -f "+libdir+proc+"/*_default.dat", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()


def compileX(index,nameX,mX,sba):
   X       = "matrix/"+nameX+"/"+str(mX)+"/"+str(sba)+"/"
   command = "./bin/mg5_aMC run/proc."+nameX+".dat"
   procdirbase = "pp-"+nameX+"-ttjets/SubProcesses/"
   S = ""
   if(nameX=="A"): S = "h"
   if(nameX=="H"): S = "h1"
   procdirs = []
   for proc in model.diagramsXX:
      procdirs.append(procdirbase+proc+"_no_"+S+"/")
   matxdir = model.workpath+X
   libdir  = matxdir+str(index)+"/"

   global parameters

   if(index==0):
      p = subprocess.Popen("rm -rf "+matxdir+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

   p = subprocess.Popen("mkdir -p "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   for proc in model.diagramsXX:
      p = subprocess.Popen("mkdir -p "+libdir+"/"+proc+"_no_"+S, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

   # enter the directory like this:
   with cd(model.mgpath):
      # make sure the old directory is removed
      p = subprocess.Popen("rm -rf "+procdirbase, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      # modify the main parameters card
      ifname = "models/HEFTff/parameters.py_ORIG"
      ofname = ifname.replace("_ORIG","")
      p = subprocess.Popen("/bin/cp -f "+ifname+" "+ofname, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
      out, err = p.communicate()

      replacements = {}
      if(nameX=="H"):
         replacements.update({' = 111.,':          ' = 100000000000.,'})
         replacements.update({' = 0.11111,':       ' = 100000000000.,'})
         replacements.update({' = 120.,':          ' = '+str(parameters[index].get("mH"))+','})
         replacements.update({' = 0.00575308848,': ' = '+str(parameters[index].get("wH"))+','})
      if(nameX=="A"):
         replacements.update({' = 111.,':          ' = '+str(parameters[index].get("mA"))+','})
         replacements.update({' = 0.11111,':       ' = '+str(parameters[index].get("wA"))+','})
         replacements.update({' = 120.,':          ' = 100000000000.,'})
         replacements.update({' = 0.00575308848,': ' = 100000000000.,'})
      replacements.update({' = 1.27,':          ' = '+str(parameters[index].get("YMC"))+','})
      replacements.update({' = 4.2,':           ' = '+str(parameters[index].get("YMB"))+','})
      replacements.update({' = 164.5,':         ' = '+str(parameters[index].get("YMT"))+','})
      replacements.update({' = 0.10567,':       ' = '+str(parameters[index].get("YMM"))+','})
      replacements.update({' = 1.778,':         ' = '+str(parameters[index].get("YMTAU"))+','})
      replacements.update({' = 172.,':          ' = 172.5,'})
      for sold in replacements.keys():
         snew = replacements[sold]
         p = subprocess.Popen('sed -i -- "s/'+sold+'/'+snew+'/g" '+ofname, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)     
         out, err = p.communicate()

      # execute the generation of the process
      p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      p = subprocess.Popen("/bin/cp -f "+procdirbase+"makefile "+procdirbase+"makefile_ORIG", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      notMade = True
      for procdir in procdirs:
         # go to make the library
         with cd(procdir):
            proc = procdir
            proc = proc.replace(procdirbase,"")
            proc = proc.replace("/","")
            procname = model.getProcName(proc)
            if(notMade):
               p = subprocess.Popen("make", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
               out, err = p.communicate()
               notMade = False
	
            ### cahnge the makefile, make and copy
            p = subprocess.Popen("/bin/cp -f ../makefile_ORIG ../makefile", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen('sed -i -e "s/MENUM)py/MENUM)'+nameX+str(index)+procname+'py/g" '+model.mgpath+procdir+'makefile', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("make matrix2"+nameX+str(index)+procname+"py.so", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("cp matrix2"+nameX+str(index)+procname+"py.so "+libdir+proc+"/", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("cp ../../Cards/*.dat "+libdir+proc+"/", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            p = subprocess.Popen("rm -f "+libdir+proc+"/*_default.dat", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()

modules = {}
def setModules(nameX,libs="All",index=-1):
   libmatrix = model.workpath+"matrix/"+nameX+"/"+str(model.mX)+"/"+str(model.sba)+"/"
   nX = len(parameters)
   S = ""
   if(nameX=="A"): S = "h"
   if(nameX=="H"): S = "h1"
   if(libs=="All" or libs=="AllX"):
      for i in range (0,nX):
         if(index!=-1 and i!=index): continue
         for proc in model.diagramsXX:
            proc += "_no_"+S
            procdir = proc
            proc = model.getProcName(proc)
            name = 'matrix2'+nameX+str(i)+proc+'py'
            sindex = str(i)
            print "changing dir to: "+libmatrix+sindex+"/"+procdir+"/"
            with cd(libmatrix+sindex+"/"+procdir+"/"):
               print "in "+os.getcwd()+", trying to import ",name
               module_info = imp.find_module(name,[libmatrix+str(i)+"/"+procdir+"/"])
               modules.update({name:imp.load_module(name, *module_info)})
               modules[name].initialise(libmatrix+str(i)+"/"+procdir+"/param_card.dat")
               print "Successfully initialised ",name
   if(libs=="X" and index!=-1):
      for proc in model.diagramsXX:
         proc += "_no_"+S
         procdir = proc
         proc = model.getProcName(proc)
         name = 'matrix2'+nameX+str(index)+proc+'py'
         sindex = str(index)
         with cd(libmatrix+sindex+"/"+procdir+"/"):
            print "in "+os.getcwd()+", trying to import ",name
            module_info = imp.find_module(name,[libmatrix+str(index)+"/"+procdir+"/"])
            modules.update({name:imp.load_module(name, *module_info)})
            modules[name].initialise(libmatrix+str(index)+"/"+procdir+"/param_card.dat")
            print "Successfully initialised ",name
   if(libs=="All" or libs=="SM"):
      for proc in model.diagramsSM:
         procdir = proc
         proc = model.getProcName(proc)
         name = 'matrix2SM'+proc+'py'
         with cd(libmatrix+"SM/"+procdir+"/"):
            print "in "+os.getcwd()+", trying to import "+name
            module_info = imp.find_module(name,[libmatrix+"SM/"+procdir+"/"])
            modules.update({name:imp.load_module(name, *module_info)})
            modules[name].initialise(libmatrix+"SM/"+procdir+"/param_card.dat")
            print "Successfully initialised ",name


def testImport(nameX,mX,sba,index=-1):
   if(index<0): setModules(nameX,"SM",index)
   else:        setModules(nameX,"X",index)
   S = ""
   if(nameX=="A"): S = "h"
   if(nameX=="H"): S = "h1"
   alphaS = 0.1 # just a random value
   nhel = 0     # sum over all helicity
   me2 = -1
   me2 = {}
   if(index<0):
      for proc in model.diagramsSM:
         p = model.ptest4
         if("P0" in proc):   p = model.ptest4
         elif("P1" in proc): p = model.ptest5
         elif("P2" in proc): p = model.ptest6
         else:
            print "ERROR - unknown SM process:", proc
            quit()
         procname = model.getProcName(proc)
         P=invert_momenta(p)
         me2.update({proc:modules['matrix2SM'+procname+'py'].get_me(P,alphaS,nhel)})                   ### calculate the SM ME^2
   else:
      for proc in model.diagramsXX:
        p = model.ptest4
        if("P0" in proc):   p = model.ptest4
        elif("P1" in proc): p = model.ptest5
        elif("P2" in proc): p = model.ptest6
        else:
           print "ERROR - unknown SM process:", proc
           quit()
        procname = model.getProcName(proc)+"no"+S
        P=invert_momenta(p)
        me2.update({proc+"_no_"+S:modules['matrix2'+nameX+str(index)+procname+'py'].get_me(P,alphaS,nhel)}) ### calculate the X ME^2
   return me2


def makeSM(nameX,mX,sba,test=False):
   compileSM(nameX,mX,sba)
   if(test):
      me2 = testImport(nameX,mX,sba)
      print "Done making library SM -> test ME^2="+str(me2)
   else:
      print "Done making library SM"


def make2HDM(nameX,mX,sba,test=False):
   for i in range(0,len(parameters)):
      compileX(i,nameX,mX,sba)
      if(test):
        me2 = testImport(nameX,mX,sba,i)
        print "Done making library "+str(i)+" -> test ME^2="+str(me2)
      else:
        print "Done making library "+str(i)


def makeAll(nameX,mX,sba,test=False):	
   make2HDM(nameX,mX,sba,test)
   makeSM(nameX,mX,sba,test)


def test2HDM():
   setParameters(model.nameX,model.mX,model.cuts,model.typeX,model.sba)
   print parameters
   makeAll(model.nameX,model.mX,model.sba,True)
