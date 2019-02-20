#!/usr/bin/env python
import os, sys
from ROOT import *

gROOT.Reset()

gROOT.SetStyle('Plain')
gStyle.SetPaintTextFormat('4.2f')
gStyle.SetPalette(1)
gStyle.SetTextSize(3)

doprint = ".png,.eps,.pdf"

def saveCanvas(canvas, name, File, OUTDIR):
        File.cd()
        a = 1;
        canvas.SetCanvasSize(int(700*a),int(500*a));
        canvas.SetWindowSize(int(700*a),int(500*a));
        canvas.Update();
        canvas.Write(name)
        os.system("mkdir -p "+OUTDIR)
        if doprint!="": 
                for ext in doprint.split(","):
                        canvas.Print(OUTDIR+"/"+name+ext)

	return

def MergeFiles(oname, list_dir):
	com = "hadd -f "+oname
	for Dir in list_dir:
		com += " "+Dir
	#print com	
	os.system(com)
	
def Rebin_2D(h1, ngx, ngy):

	nbinsx = h1.GetXaxis().GetNbins()
	nbinsy = h1.GetYaxis().GetNbins()
	xmin   = h1.GetXaxis().GetXmin()
	xmax   = h1.GetXaxis().GetXmax()
	ymin   = h1.GetYaxis().GetXmin()
	ymax   = h1.GetYaxis().GetXmax()
	nx = nbinsx/ngx;
	ny = nbinsy/ngy;
	print nx,xmin,xmax,ny,ymin,ymax
	h1.SetBins(nx,xmin,xmax,ny,ymin,ymax)
	
	#loop on all bins to reset contents and errors
	for biny in range(nbinsy):
		for binx in range(nbinsx):
			ibin = h1.GetBin(binx,biny)
			h1.SetBinContent(ibin,0)
	
	#loop on all bins and refill
	for biny in range(nbinsy):
		by  = h1.GetYaxis().GetBinCenter(biny)
		iy  = h1.GetYaxis().FindBin(by)
		for binx in range(nbinsx):
			bx = h1.GetXaxis().GetBinCenter(binx)
			ix  = h1.GetXaxis().FindBin(bx)
			bin = h1.GetBin(binx,biny)
			ibin= h1.GetBin(ix,iy)
			cu  = h1.GetBinContent(bin)
			h1.AddBinContent(ibin,cu)

	return h1

def effRates(inputDir):
	
	merge = 1
	extraPlots = 1
	verbose = 0
	
	channels  = []
	channels += [('resolved','e' )]
	#channels += [('resolved','mu')]
	#channels += [('boosted', 'e' )]
	#channels += [('boosted', 'mu')]

	outfile = TFile('eff_ttbar.root','RECREATE')
	#outfile = TFile('eff_wjets_wev.root','RECREATE')
	
	for ichan in channels:

		#if(ichan[0]=='resolved'):	regime = "r";
		#elif(ichan[0]=='boosted'):	regime = "b";
				
		hadd_in=""
		hadd_in = inputDir+ichan[0]+"_"+ichan[1]+"_ttbar*root"
		#hadd_in = inputDir+regime+ichan[1]+"_nom_wev.root"
		
		ttbarFile = ichan[0]+"_ttbar_"+ichan[1]+".root "
		
		print "hadd -f " + ttbarFile + hadd_in
		os.system("hadd -f " + ttbarFile + hadd_in)
		
		inputfile = TFile(ttbarFile,'READ')
		
		
		parametrization1D =  ["lepPt", "minDeltaR", "cosDPhi", "MET", "d0sig", "Dz0sin", "mwt", "mwt_met", "nJets", "nTrkBtagJets","ptvarcone","topoetcone","topoetcone_effBins"]
		#parametrization1D += ["lepPhi", "MET_phi", "DeltaPhi", "minDeltaR_tjet_effBins","minDeltaR_tjet", "mwt_effBins", "mwt_met_effBins", "Dz0sin"]
		
		for ipar in parametrization1D:
		
			histo_t = inputfile.Get("eff_"+ipar).Clone()		   
			histo_l = inputfile.Get("eff_"+ipar+"_Loose").Clone()
			histo_l.Add(histo_t,+1)
						
			#histo_t.Scale(lumi)
			if extraPlots:	histo_t.Write("hP_eff_"+ipar+"_tight_"+ichan[0]+'_'+ichan[1])
		
			#histo_l.Scale(lumi)
			if extraPlots:	histo_l.Write("hP_eff_"+ipar+"_loose_"+ichan[0]+'_'+ichan[1])
						
			h_ratio = histo_t
			h_ratio.Divide(histo_t, histo_l, 1.0, 1.0, "B")
				
			c = TCanvas("c_eff_"+ipar+"_"+ichan[0]+"_"+ichan[1])
                	maxi = 0
			mini = 0
			h_ratio.Draw("e")
			h_ratio.SetStats(0)
			h_ratio.GetYaxis().SetRangeUser(0,1)
			h_ratio.SetMarkerStyle(20)
			h_ratio.SetMarkerSize(0.8)		 
			c.SetGridy(1)
			c.SetLogx(0)
			if ipar in ["lepEta", "cos_metPhi_lepPhi"]:
				c.SetLogx(0)
			c.SetLogy(0)
			c.Update()
               		c.Modified()
			outfile.cd()
			h_ratio.Write("realRate_"+ipar+"_"+ichan[0]+"_"+ichan[1])	               			
			saveCanvas(c, "h_real_"+ipar+"_"+ichan[0]+"_"+ichan[1], outfile, "real_plots")
		
		parametrization2D = ["LepPt_DR", "mwt_met_map", "mwt_met_map_lowDR", "mwt_met_map_medDR", "mwt_met_map_highDR", "lepPt_topoetcone"]
		for ipar in parametrization2D:

			if(ipar.find("LepPt_DR")!=-1):
				histo_t = inputfile.Get("eff_"+ipar).Clone() 		   
				histo_l = inputfile.Get("eff_"+ipar+"_Loose").Clone()
				histo_l.Add(histo_t,+1)	
				name = 'pTdr'		

			else:
				histo_t = inputfile.Get("eff_"+ipar).Clone() 		   
				histo_l = inputfile.Get("eff_"+ipar+"_Loose").Clone()
				histo_l.Add(histo_t,+1)
				name = ipar
			
						
			# ---> 2D fake rates: leptPt vs leptEta
			if extraPlots:	histo_t.Write("hP_2D_"+ipar+"_tight_"+ichan[0]+'_'+ichan[1])
		
			if extraPlots:	histo_l.Write("hP_2D_"+ipar+"_loose_"+ichan[0]+'_'+ichan[1])
		
			h2D_ratio = histo_t.Clone()
			h2D_ratio.Divide(histo_t, histo_l, 1.0, 1.0, "B")
		
			gStyle.SetPaintTextFormat("1.2f")
			c2 = TCanvas("c_2D_"+ipar+"_"+ichan[0]+"_"+ichan[1])
			c2.cd()
			h2D_ratio.Draw()
			h2D_ratio.Draw("colz TEXT e1 SAME")
			h2D_ratio.SetStats(0)
			
			#h2D_ratio.GetZaxis().SetRangeUser(0.7, 1.0)
			h2D_ratio.GetXaxis().SetRangeUser(30, 300)
			h2D_ratio.GetYaxis().SetRangeUser(0.4, 7)
			#if(ichan[0]=='boosted'):	h2D_ratio.GetYaxis().SetRangeUser(0., 1.5)
			#else:				h2D_ratio.GetYaxis().SetRangeUser(0., 5.0)
			
			gPad.RedrawAxis()	
			c2.SetLogx(1)
			#c2.SetLogy(1)
			c2.Update()
               		c2.Modified()	
		
			l = TLatex()
			l.SetNDC()
			l.SetTextFont(72)
			l.SetTextSize(0.03)
			l.SetTextColor(kBlack)
			#l.DrawLatex(0.1,0.92, "#intLdt ="+`lumi/1000.`[:3]+" fb^{-1}")	
			h2D_ratio.Write('eff_'+name+'_'+ichan[0]+'_'+ichan[1])
			saveCanvas(c2, "2Dh_"+name+"_"+ichan[0]+"_"+ichan[1], outfile, "real_plots")
		

	return
	
	
def fakeRates(inputDir, lumi):

	merge = 1
	extraPlots = 1
	verbose = 0
		
	channels  = []
	
	
	channels += [('resolved','e' )]
	#channels += [('resolved','mu')]
	#channels += [('boosted', 'e' )]
	#channels += [('boosted', 'mu')]
		
	iPad = 0
	
	outfile = TFile("fake.root","RECREATE")
		
	for ichan in channels:
	
		hadd_in= []		
		
		# --> Merging the bkg
		for ibkg in ["ttbar", "Wev", "Wmv", "Wtv", "Zjets", "st", "VV"]:			
			
			if ibkg=="Zjets" :
				for File in os.popen("ls "+inputDir+ichan[0]+"_"+ichan[1]+"_Z*root ").readlines():
					hadd_in.append(File[:-1])
			elif ibkg=="st":
				for File in os.popen("ls "+inputDir+ichan[0]+"_"+ichan[1]+"_st.root ").readlines():
					hadd_in.append(File[:-1])	
			elif ibkg=="VV":
				for File in os.popen("ls "+inputDir+ichan[0]+"_"+ichan[1]+"_VV.root ").readlines():
					hadd_in.append(File[:-1])
                        elif ibkg=="Wtv":
                                for File in os.popen("ls "+inputDir+ichan[0]+"_"+ichan[1]+"_Wtv.root ").readlines():
                                        hadd_in.append(File[:-1])
                        else:
                                hadd_intmp= []
                                for File in os.popen("ls "+inputDir+ichan[0]+"_"+ichan[1]+"_"+ibkg+"*.root ").readlines():
                                        hadd_intmp.append(File[:-1])                            
                                hadd_Path = inputDir+ichan[0]+"_"+ichan[1]+"_"+ibkg+".root"
                                MergeFiles(hadd_Path, hadd_intmp)
                                for File in os.popen("ls "+inputDir+ichan[0]+"_"+ichan[1]+"_"+ibkg+".root ").readlines():
                                        hadd_in.append(File[:-1])

		hadd_bkgPath = ichan[0]+"_BKG_"+ichan[1]+".root"
		if merge:	MergeFiles(hadd_bkgPath, hadd_in)
		
		# --> Merging the datasets
		hadd_in=[]		
		print "ls "+inputDir+ichan[0]+"_"+ichan[1]+"_DT*root "
		
		for File in os.popen("ls "+inputDir+ichan[0]+"_"+ichan[1]+"_DT*root ").readlines():
			hadd_in.append(File[:-1])
		
		hadd_dataPath =  ichan[0]+"_DATA_"+ichan[1]+".root"
		if merge:	MergeFiles(hadd_dataPath, hadd_in)
		
		# --> Check QCD at low MWT
		
		bkgFile  = TFile(hadd_bkgPath,'READ')
		dataFile = TFile(hadd_dataPath,'READ')
		
		varList = [] 
		varList += ['mwt', 'lepPt_effBins','lepPt','minDeltaR_effBins','minDeltaR','minDeltaR_tjet_effBins','minDeltaR_tjet','closJetPt']
		varList += ['MET', 'lepPhi', 'MET_phi', 'DeltaPhi', 'cos_metPhi_lepPhi', 'cos_metPhi_lepPhi_effBins', 'lepEta','nTrkBtagJets','nJets', 'closJetJVT', 'mwt_met', 'Dz0sin', 'd0sig']
		
		varList += ["ptvarcone", "topoetcone", "topoetcone_effBins"]
		#varList += ['mwt_effBins1', 'mwt_effBins2', 'mwt_effBins3', 'mwt_effBins4', 'mwt_effBins5', 'mwt_effBins6']
		
		#varList += ['chi2']
		#varList += ['ljet_m', 'ljet_pt', 'ljet_tau32', 'ljet_tau32_wta']
		
		for var in varList:
			
			print var		
					
			if (var.find("chi2")!=-1):	
				if (ichan[0]=='b'):	
					continue		
		
			#tight
			c_t = TCanvas("c_"+ichan[0]+'_'+var+'_tight_'+ichan[1])
			leg = TLegend(0.6, 0.6, 0.8, 0.88, ichan[0]+" "+ichan[1]+" channel (tight):")
                	leg.SetTextSize(0.028)
                	leg.SetFillColor(10)
                	leg.SetBorderSize(0)
		
			h_data = dataFile.Get("fake_"+var).Clone()
			if (var.find("mwt")!=-1):
				h_data.GetXaxis().SetRangeUser(0, 200)
			if (var.find("MET")!=-1):
				h_data.GetXaxis().SetRangeUser(0, 200)
			if (var.find("d0sig")!=-1):
				#h_data.Rebin(2)
				h_data.GetXaxis().SetRangeUser(-100, 100)	
			h_data.Draw("e")
			h_data.SetMarkerStyle(20)
			h_data.SetStats(0)
			leg.AddEntry(h_data, "Data: "+`h_data.Integral()`[:7], "P")
			if(verbose):	print ichan, "tight:"
			if(verbose):	print 'data &' , `h_data.Integral()`[:7], '\\\ \\hline'
			
			pallete   = [2, 3, 3, 3, 6, 6, 6, 7, 11]
			styleLine = [1, 1 ,2, 3, 1, 2, 3, 1, 1]
			name 	  = ["ttbar", "We#nu", "W#mu#nu", "W#tau#nu", "Zee", "Z#mu#mu", "Z#tau#tau", "st", "Diboson"]
			h_ibkg = []
			i=0
			
			#bkgList = ["ttbar", "Wev", "Wmv", "Wtv", "Zee", "Zmm", "Ztt", "st", "VV"]
                        bkgList = ["ttbar", "Wev", "Wmv", "Wtv", "Zmm", "Ztt", "st", "VV"]
						
			for ibkg in bkgList:
				
				tmp_bkgFile  = TFile(inputDir+ichan[0]+"_"+ichan[1]+"_"+ibkg+".root",'READ')	
				if(verbose):	print ibkg, i, var
				#if not tmp_bkgFile.IsOpen():	continue
								
				h_ibkg.append(tmp_bkgFile.Get("fake_"+var).Clone())
				h_ibkg[i].SetDirectory(0)
				h_ibkg[i].Scale(lumi)
				
				
				#print h_ibkg[i].Integral()
				
				if h_ibkg[i].Integral()>0:	
					if h_ibkg[i].Integral()>0.0001: #for better visualization
						if(verbose):	print ibkg, '&' , `h_ibkg[i].Integral()`[:7], '\\\ \\hline'	
						h_ibkg[i].Draw("histo SAME")
						if (var.find("mwt")!=-1):
							h_ibkg[i].GetXaxis().SetRangeUser(0, 200)
						if (var.find("MET")!=-1):
							h_ibkg[i].GetXaxis().SetRangeUser(0, 200)
						if (var.find("d0sig")!=-1):
							#h_ibkg[i].Rebin(2)
							h_ibkg[i].GetXaxis().SetRangeUser(-100, 100)	
						h_ibkg[i].SetStats(0)
					h_ibkg[i].SetLineColor(pallete[i]) 
					h_ibkg[i].SetLineStyle(styleLine[i])
					h_ibkg[i].SetLineWidth(3)					
					#print h_ibkg[i].Integral(), i
					leg.AddEntry(h_ibkg[i],  name[i]+": "+`h_ibkg[i].Integral()`[:7], "L")
				i += 1
			
			h_bkg = bkgFile.Get("fake_"+var).Clone()
			h_bkg.Scale(lumi)
			h_bkg.Draw("histo SAME")
			if (var.find("mwt")!=-1):
				h_bkg.GetXaxis().SetRangeUser(0, 200)
			if (var.find("MET")!=-1):
				h_bkg.GetXaxis().SetRangeUser(0, 200)
			if (var.find("d0sig")!=-1):
				#h_bkg.Rebin(2)
				h_bkg.GetXaxis().SetRangeUser(-100, 100)
				
			h_bkg.SetLineColor(4) 
			h_bkg.SetLineWidth(3)
			h_bkg.SetStats(0)
						
			leg.AddEntry(h_bkg,  "total bkg: "+`h_bkg.Integral()`[:7], "L")
			leg.Draw()
			if(verbose):	print 'total bkg &' , `h_bkg.Integral()`[:7], '\\\ \\hline'
			
			
			if h_data.GetMaximum()>h_bkg.GetMaximum():	h_data.SetMaximum(h_data.GetMaximum()*1.3)
			else: 						h_data.SetMaximum(h_bkg.GetMaximum()*1.3)
		
			c_t.SetLogy(0)
			if (var.find("lepPt")!=-1):	c_t.SetLogx(1)
			c_t.SetGridy(1)
			c_t.Update()
                	c_t.Modified()  	    
			saveCanvas(c_t, 'h_'+ichan[0]+'_'+var+'_tight_'+ichan[1], outfile, 'fake_plots')
			c_t.SetLogy(1)
			c_t.Modified()  
			saveCanvas(c_t, 'hLOG_'+ichan[0]+'_'+var+'_tight_'+ichan[1], outfile, 'fake_plots')
			
			#loose
			c_l = TCanvas("c_"+ichan[0]+'_'+var+'_loose_'+ichan[1])
			leg = TLegend(0.6, 0.6, 0.8, 0.88, ichan[0]+" "+ichan[1]+" channel (anti Tight):")
                	leg.SetTextSize(0.028)
	                leg.SetFillColor(10)
        	        leg.SetBorderSize(0)
		
			h_data = dataFile.Get('fake_'+var+'_Loose').Clone()
			h_data.Draw("e")
			if (var.find("mwt")!=-1):
				h_data.GetXaxis().SetRangeUser(0, 200)
			if (var.find("MET")!=-1):
				h_data.GetXaxis().SetRangeUser(0, 200)
			if (var.find("d0sig")!=-1):
				#h_data.Rebin(2)
				h_data.GetXaxis().SetRangeUser(-100, 100)	
				
			h_data.SetMarkerStyle(20)
			h_data.SetStats(0)
			leg.AddEntry(h_data, "Data: "+`h_data.Integral()`[:7], "P")
			if(verbose):	print ichan, "Loose:"
			if(verbose):	print 'data &' , `h_data.Integral()`[:7], '\\\ \\hline'
			
			pallete   = [2, 3, 3, 3, 6, 6, 6, 7, 11]
			styleLine = [1, 1 ,2, 3, 1, 2, 3, 1, 1]
			name 	  = ["ttbar", "We#nu", "W#mu#nu", "W#tau#nu", "Zee", "Z#mu#mu", "Z#tau#tau", "st", "Diboson"]
			h_ibkg = []
			i=0
			#for ibkg in ["ttbar", "Wev", "Wmv", "Wtv", "Zee", "Zmm", "Ztt", "st", "VV"]:
                        for ibkg in ["ttbar", "Wev", "Wmv", "Wtv", "Zmm", "Ztt", "st", "VV"]:
				tmp_bkgFile  = TFile(inputDir+ichan[0]+"_"+ichan[1]+"_"+ibkg+".root",'READ')
				h_ibkg.append(tmp_bkgFile.Get("fake_"+var+'_Loose').Clone())
				h_ibkg[i].SetDirectory(0)
				h_ibkg[i].Scale(lumi)
				if h_ibkg[i].Integral()>0:	
					if h_ibkg[i].Integral()>0.0001:	
						h_ibkg[i].Draw("histo SAME")
						if (var.find("mwt")!=-1):
							h_ibkg[i].GetXaxis().SetRangeUser(0, 200)
						if (var.find("MET")!=-1):
							h_ibkg[i].GetXaxis().SetRangeUser(0, 200)
						if (var.find("d0sig")!=-1):
							#h_ibkg[i].Rebin(2)
							h_ibkg[i].GetXaxis().SetRangeUser(-100, 100)	
						if(verbose):	print ibkg, '&' , `h_ibkg[i].Integral()`[:7], '\\\ \\hline'
						h_ibkg[i].SetStats(0)
					h_ibkg[i].SetLineColor(pallete[i]) 
					h_ibkg[i].SetLineStyle(styleLine[i])
					h_ibkg[i].SetLineWidth(3)					
					#print h_ibkg[i].Integral(), i
					leg.AddEntry(h_ibkg[i],  name[i]+": "+`h_ibkg[i].Integral()`[:7], "L")
				tmp_bkgFile.Close()	
				i += 1
						
			h_bkg = bkgFile.Get('fake_'+var+'_Loose').Clone()
			h_bkg.Scale(lumi)
			h_bkg.Draw("histo SAME")
			if (var.find("mwt")!=-1):
				h_bkg.GetXaxis().SetRangeUser(0, 200)
			if (var.find("MET")!=-1):
				h_bkg.GetXaxis().SetRangeUser(0, 200)
			if (var.find("d0sig")!=-1):
				#h_bkg.Rebin(2)
				h_bkg.GetXaxis().SetRangeUser(-100, 100)
				
			h_bkg.SetLineColor(4) 
			h_bkg.SetLineWidth(3)
			h_bkg.SetStats(0)		
		
                	if h_data.GetMaximum()>h_bkg.GetMaximum():	h_data.SetMaximum(h_data.GetMaximum()*1.3)
			else: 						h_data.SetMaximum(h_bkg.GetMaximum()*1.3)	
			
			leg.AddEntry(h_bkg,  "bkg: "+`h_bkg.Integral()`[:7], "L")
			leg.Draw()
			if(verbose):	print 'total bkg &' , `h_bkg.Integral()`[:7], '\\\ \\hline'
			c_l.SetLogy(0)
			if (var.find("lepPt")!=-1):	c_l.SetLogx(1)
			
			c_l.SetGridy(1)
			c_l.Update()
        	        c_l.Modified()  	    
			saveCanvas(c_l, 'h_'+ichan[0]+'_'+var+'_loose_'+ichan[1], outfile, 'fake_plots')
			c_l.SetLogy(1)
			c_l.Modified() 
			saveCanvas(c_l, 'hLOG_'+ichan[0]+'_'+var+'_loose_'+ichan[1], outfile, 'fake_plots')
		# --> Fake rate
		parametrization1D = []
		parametrization1D =  ["lepPt_effBins", "lepPt", "lepEta", "minDeltaR_effBins", "minDeltaR", "closJetPt_effBins", "cos_metPhi_lepPhi", "cos_metPhi_lepPhi_effBins","MET", "d0sig"]
		#parametrization1D += ["minDeltaR_highEta", "lepPt_highEta", "cos_metPhi_lepPhi_highEta", "minDeltaR_lowEta", "lepPt_lowEta", "cos_metPhi_lepPhi_lowEta"]
		parametrization1D += ["ptvarcone", "topoetcone", "topoetcone_effBins"]
		#parametrization1D += ["lepPt_highmWt", "lepPt_medmWt", "lepPt_lowmWt", "cos_metPhi_lepPhi_highmWt", "cos_metPhi_lepPhi_medmWt", "cos_metPhi_lepPhi_lowmWt"]
		#parametrization1D += ["lepPhi", "MET_phi", "DeltaPhi", "minDeltaR_tjet_effBins","minDeltaR_tjet", "mwt_effBins", "mwt_met_effBins", "Dz0sin", "nJets", "nTrkBtagJets"]
		#parametrization1D += ["lepPt_lowDR", "lepPt_highDR", "minDeltaR_lowDR", "minDeltaR_highDR", "cos_metPhi_lepPhi_lowDR", "cos_metPhi_lepPhi_highDR"]
		parametrization1D += ["mwt_effBins", "MET_effBins"]
#		parametrization1D += ["mwt_met_highDR", "mwt_met_lowDR", "met_highDR", "met_lowDR"]
               # parametrization1D += ["met_highDR"]
		
		for ipar in parametrization1D:
		
			histo_t = []
			histo_l = []
			
			iPad=0
			for item in ["BKG","DATA"]:	
			
				if item=="BKG":		inputfile = bkgFile
				elif item=="DATA":	inputfile = dataFile

				histo_t.append( inputfile.Get("fake_"+ipar).Clone() )			   
				histo_l.append( inputfile.Get("fake_"+ipar+"_Loose").Clone() )
				histo_l[iPad].Add(histo_t[iPad],+1)
				
				#hdummy = histo_l[iPad].Clone()
				#histo_l[iPad].Rebin(hdummy.GetXaxis().GetNbins())
				
				#hdummy = histo_t[iPad].Clone()
				#histo_t[iPad].Rebin(hdummy.GetXaxis().GetNbins())
				
				iPad+=1
				
			histo_t[0].Scale(lumi)
			h1 = histo_t[1]
			h1.Add(histo_t[0],-1)
			if extraPlots:	h1.Write("hP_"+ipar+"_tight_"+ichan[0]+'_'+ichan[1])
		
			histo_l[0].Scale(lumi)
			h2 = histo_l[1]
			h2.Add(histo_l[0],-1)		
			if extraPlots:	h2.Write("hP_"+ipar+"_loose_"+ichan[0]+'_'+ichan[1])
						
			h_ratio = h1
			h_ratio.Divide(h1, h2, 1.0, 1.0, "B")
				
			c = TCanvas("c_"+"fake_"+ipar+"_"+ichan[0]+"_"+ichan[1])
                	maxi = 0
			mini = 0
			h_ratio.Draw("e")
			h_ratio.SetStats(0)
			h_ratio.GetYaxis().SetRangeUser(0,1)
			h_ratio.SetMarkerStyle(20)
			h_ratio.SetMarkerSize(0.8)		 
			c.SetGridy(1)
			c.SetLogx(0)
			if ipar in ["lepEta", "cos_metPhi_lepPhi"]:
				c.SetLogx(0)
			c.SetLogy(0)
			c.Update()
               		c.Modified()
			h_ratio.Write("fakeRate_"+ipar+"_"+ichan[0]+"_"+ichan[1])	               			
			saveCanvas(c, "h_fake_"+ipar+"_"+ichan[0]+"_"+ichan[1], outfile, "fake_plots")
		
		parametrization2D = []
		#parametrization2D =  ["lepPt_lepEta", "lepPt_closJetPt", "lepPt_closJetPt_lowDR", "lepPt_closJetPt_medDR", "lepPt_closJetPt_highDR"]
		#parametrization2D += ["lepPt_cosDPhi", "lepPt_cosDPhi_lowDR", "lepPt_cosDPhi_medDR", "lepPt_cosDPhi_highDR"]
		#parametrization2D += ["lepPt_cosDPhi_lowEta", "lepPt_cosDPhi_highEta", "lepPt_closJetPt_lowEta", "lepPt_closJetPt_highEta"]
		#parametrization2D += ["lepPt_cosDPhi_lowmWt", "lepPt_cosDPhi_medmWt", "lepPt_cosDPhi_highmWt"]
		#parametrization2D += ["lepPt_met", "lepPt_met_lowDR", "lepPt_met_highDR"]
		parametrization2D += ["lepPt_minDeltaR","lepPt_topoetcone"]
		#parametrization2D += ["minDeltaR_met_highLepPt", "minDeltaR_met_lowLepPt", "lepPt_closJetPt"]
		parametrization2D += ['mwt_met_map', 'mwt_met_map_lowDR', 'mwt_met_map_medDR', 'mwt_met_map_highDR']
		for ipar in parametrization2D:
		
			histo_t = []
			histo_l = []
			
			iPad=0
			for item in ["BKG","DATA"]:	
			
				if item=="BKG":		inputfile = bkgFile
				elif item=="DATA":	inputfile = dataFile

				histo_t.append( inputfile.Get("fake_"+ipar).Clone() )			   
				histo_l.append( inputfile.Get("fake_"+ipar+"_Loose").Clone() )
				histo_l[iPad].Add(histo_t[iPad],+1)
				
				iPad+=1
			
			# ---> 2D fake rates: leptPt vs leptEta
			histo_t[0].Scale(lumi)
			h1_2D = histo_t[1].Clone()
			h1_2D.Add(histo_t[0],-1)	
			if extraPlots:	h1_2D.Write("hP_2D_"+ipar+"_tight_"+ichan[0]+'_'+ichan[1])
		
			histo_l[0].Scale(lumi)		
			h2_2D = histo_l[1].Clone()
			h2_2D.Add(histo_l[0],-1) 
			if extraPlots:	h2_2D.Write("hP_2D_"+ipar+"_loose_"+ichan[0]+'_'+ichan[1])
		
			if extraPlots:
				c = TCanvas("c_"+"hP_2D_"+ipar+"_t_"+ichan[0]+"_"+ichan[1])
				ht_DTMC = histo_t[1].Clone()
				ht_DTMC.Divide(histo_t[0].Clone())
				if ipar.find("mwt_met_map")!=-1:	ht_DTMC.Draw("colz")
				else:					ht_DTMC.Draw("colz TEXT e")
				ht_DTMC.SetStats(0)
				c.SetLogx(1)
				ht_DTMC.Write("hP_2D_DTMCr_"+ipar+"_tight_"+ichan[0]+"_"+ichan[1])
				saveCanvas(c, "c_hP_2D_DTMCr_"+ipar+"_tight_"+ichan[0]+"_"+ichan[1], outfile, 'fake_plots')
			
				c = TCanvas("c_"+"hP_2D_"+ipar+"_l_"+ichan[0]+"_"+ichan[1])
				hl_DTMC = histo_l[1].Clone()
				hl_DTMC.Divide(histo_l[0].Clone())
				if ipar.find("mwt_met_map")!=-1:	hl_DTMC.Draw("colz")
				else:					hl_DTMC.Draw("colz TEXT e")
				hl_DTMC.SetStats(0)
				c.SetLogx(1)
				hl_DTMC.Write("hP_2D_DTMCr_"+ipar+"_loose_"+ichan[0]+"_"+ichan[1])
				saveCanvas(c, "c_hP_2D_DTMCr_"+ipar+"_loose_"+ichan[0]+"_"+ichan[1], outfile, 'fake_plots')
		
			h2D_ratio = h1_2D.Clone()
			h2D_ratio.Divide(h1_2D, h2_2D, 1.0, 1.0, "B")
		
			gStyle.SetPaintTextFormat("1.2f")
			c2 = TCanvas("c_"+"2D_"+ipar+"_"+ichan[0]+"_"+ichan[1])
			c2.cd()
			h2D_ratio.Draw()
			h2D_ratio.Draw("colz TEXT e1 SAME")
			h2D_ratio.SetStats(0)
			gPad.RedrawAxis()	
			c2.SetLogx(0)
						 
			if ipar.find("met")!=-1:
				h2D_ratio.GetXaxis().SetRangeUser(20, 100)
				h2D_ratio.GetYaxis().SetRangeUser(0, 100)
			
			if ipar.find("lepPt")!=-1:	
			#if ipar in ["lepPt_cosDPhi", "lepPt_cosDPhi_lowDR", "lepPt_cosDPhi_highDR", "lepPt_cosDPhi_medDR", "lepPt_minDeltaR"]:	
				#h2D_ratio.GetXaxis().SetRangeUser(25, 100)
				#h2D_ratio.GetYaxis().SetRangeUser(0,7)
				c2.SetLogx(1)
				c2.SetLogy(0)
				
			c2.Update()
               		c2.Modified()	
		
			l = TLatex()
			l.SetNDC()
			l.SetTextFont(72)
			l.SetTextSize(0.03)
			l.SetTextColor(kBlack)
			l.DrawLatex(0.1,0.92, "#intLdt ="+`lumi/1000.`[:3]+" fb^{-1}")	
			h2D_ratio.Write("2Dfake_"+ipar+"_"+ichan[0]+"_"+ichan[1])
			saveCanvas(c2, "2Dh_"+ipar+"_"+ichan[0]+"_"+ichan[1], outfile, "fake_plots")	
	
	bkgFile.Close()
	dataFile.Close()	
	outfile.Close()
	return



#--------------------------------#
#       Multijet estimation      #
#--------------------------------#

#Produce eff rate plots

#inputDir = '/AtlasDisk/users/romano/fakeStudies/2.3.41/LPCTools/ProduceMiniTuple/030816_ePreTag_v1.0_2j_realRates/'
inputDir = '/AtlasDisk/users/sanmay/TTBar/AnalysisTop-2.4.16/LPCTools_New/ProduceMiniTuple/COMB_REALRATES/'

if 0:
	effRates(inputDir)


#inputDir = '/AtlasDisk/users/romano/fakeStudies/2.3.41/LPCTools/ProduceMiniTuple/250616_muPreTag_v1.0_2j_fakeRates/'
#inputDir = '/AtlasDisk/users/romano/fakeStudies/2.3.41/LPCTools/ProduceMiniTuple/250616_muTag_v1.0_2j_fakeRates/'

#inputDir = '/AtlasDisk/users/romano/fakeStudies/2.3.41/LPCTools/ProduceMiniTuple/030816_ePreTag_v1.0_2j_fakeRates/'
inputDir = '/AtlasDisk/users/sanmay/TTBar/AnalysisTop-2.4.16/LPCTools_New/ProduceMiniTuple/080816_eTag_v1.1_2j_fakeRates_2016_1btag/'

#lumi = 3193.68 #pb-1
#lumi =  5807.5 #pb-1
#lumi =  3200+3500 #pb-1
lumi = 3200+11571


if 1:	
	fakeRates(inputDir, lumi)


