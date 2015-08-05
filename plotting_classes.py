# import plottingUtils as putils
import signalUtils as sutils
import plotDetails as pdets
import numpy as np
import ROOT as r
import math as ma

###---------------------------------------------------------------------------------------------###
###---------------------------------------------------------------------------------------------###
# set some global root variables for the instance imported above

r.gROOT.SetBatch(r.kTRUE)
r.gStyle.SetOptStat(0)
r.TH1.SetDefaultSumw2(1)
r.TH2.SetDefaultSumw2(1)

###---------------------------------------------------------------------------------------------###
###---------------------------------------------------------------------------------------------###

def safe_divide(num = None, denom = None):
    if denom > 0.:
        return float(num/denom)
    else:
        return 0.

def safe_sqrt(num = None):
    if num > 0.:
        return ma.sqrt(num)
    else:
        return 0.

def syst_picker(pos = (), neg = ()):

    class syst(object):
        def __init__(self, data = ()):
            self.val = data[0]
            self.err = data[1]
            self.rerr = safe_divide(data[1], data[0])
        def __str__(self):
            return "syst_obj: %.2f +/- %2f (%.4fper)" % (self.val, self.err, self.rerr*100.)

    # assign largest central value to 'a'
    if pos[0] >= neg[0]:
        a = syst(pos)
        b = syst(neg)
    else:
        a = syst(neg)
        b = syst(pos)

    # if rel err on a is bigger than twice rel err on b
    if a.rerr > 2*b.rerr:
        # and if lower bound on a is less than b, use b
        if (a.val-a.err < b.val):
            print "DODGEY!"
            return b.val, b.err

    return a.val, a.err


class effMap(object):
    '''container for efficiency map'''
    def __init__(self, nom = None, denom = None, nom_err = None, denom_err = None, noweight = False):
        self._hist = None
        self._errHist = None
        self._relErrHist = None
        self._nom = nom
        self._nomerr = nom_err
        self._denom = denom
        self._denomerr = denom_err
        self._noweight = noweight #if True, then calculate errors from noweight yields provided
        self._mean = 0.
        self._max = 0.
        self._min = 1.
        self._rms = 0.

        # if only a nom is passed (i.e. eff externally calculated)
        if not denom:
            self._hist = nom
            # if a err hist is also passed
            if nom_err:
                self._errHist = nom_err

        self.process()

        ###
        # 1. Add check that both hists are self-consistent

    def process(self):
        '''make eff hist and calculate values'''

        if self._denom:
            self._hist = self._nom.Clone()
            self._hist.Divide(self._denom)
        self._nBins = self._hist.GetNbinsX() * self._hist.GetNbinsY() + 1000
        vals = []
        for i in range(1, self._nBins):
            val = self._hist.GetBinContent(i)
            if val > 0:
                vals.append(val)
        
        val_array = np.array(vals)
        
        if len(vals):
            self._rms   = np.var(val_array)
            self._mean  = np.mean(val_array)
            self._min   = np.amin(val_array)
            self._max   = np.amax(val_array)
        else:
            self._rms   = 0.
            self._mean  = 0.
            self._min   = 0.
            self._max   = 0.
        
        if self._nom and self._denom:
            # create a histogram of eff errors
            self.makeErrorHist()

        self.makeRelErrorHist()

    def makeErrorHist(self):
        '''caclulate bin by bin stat absolute error'''

        # clone a hist to fill with errors
        self._errHist = self._nom.Clone()

        # loop all bins
        for i in range(1, self._nBins):

            # skip points that aren't in the scan
            # note this skips points with no num, a denom, and assoc errors...
            if self._nom.GetBinContent(i) == 0:
                self._errHist.SetBinContent(i, 0.)
                continue

            # if no existing error hists are passed, calc from scratch
            if not self._nomerr and not self._denomerr:
                # print ">>> No error histograms passed to effMap."
                # if denom hist, sum stat errors
                if self._denom:
                    nom_val     = self._nom.GetBinContent(i)
                    denom_val   = self._denom.GetBinContent(i)

                    # if nom_val > 0. and denom_val > 0.:
                        
                    #     # calculate the statistical error for each yield
                    #     nom_err     = 1./ma.sqrt(nom_val)
                    #     denom_err   = 1./ma.sqrt(denom_val)

                    #     nom_ratio = nom_err/nom_val
                    #     denom_ratio = denom_err/denom_val
                    
                    #     eff_val = nom_val/denom_val
                    #     eff_err = eff_val*ma.sqrt(nom_ratio+denom_ratio) # should be squared sum
                    
                    # else:
                    #     # set to zero if either hist has no entries
                    #     eff_err = 0.

                    nom_err     = safe_divide(1., safe_sqrt(nom_val))
                    denom_err   = safe_divide(1., safe_sqrt(denom_val))

                    nom_ratio   = safe_divide(nom_err, nom_val)
                    denom_ratio = safe_divide(denom_err, denom_val)

                    eff_val = safe_divide(nom_val, denom_val)
                    eff_err = eff_val*safe_sqrt(ma.pow(nom_ratio, 2) + ma.pow(denom_ratio, 2))

                else:
                    # if no denom hist, just use error cald'd from nom hist
                    eff_val = self._nom.GetBinContent(i)
                    # BUG - is this correct?
                    eff_err = safe_divide(1., safe_sqrt(eff_val))
                    
            else:
                # in here if error histograms were passed

                if self._noweight:
                    # if self._noweight is True, then the noweight yields have been passed as error
                    # histograms
                    nom_val_err     = self._nomerr.GetBinContent(i)
                    denom_val_err   = self._denomerr.GetBinContent(i)

                    # error is calculated from each as the pure MC stat uncert: err = 1/sqrt(N)
                    nom_err = safe_divide(1., safe_sqrt(nom_val_err))
                    denom_err = safe_divide(1., safe_sqrt(denom_val_err))
                else:
                    # if error hists are passed, combine
                    nom_err     = self._nomerr.GetBinContent(i)
                    denom_err   = self._denomerr.GetBinContent(i)

                nom_val     = self._nom.GetBinContent(i)
                denom_val   = self._denom.GetBinContent(i)

                nom_ratio = safe_divide(nom_err, nom_val)
                denom_ratio = safe_divide(denom_err, denom_val)

                eff_val = safe_divide(nom_val, denom_val)
                eff_err = eff_val * safe_sqrt( ma.pow(nom_ratio, 2) + ma.pow(denom_ratio, 2) )
        
            self._errHist.SetBinContent(i, eff_err)

    def makeRelErrorHist(self):
        '''make histo of relative errs'''
        self._relErrHist = self._errHist.Clone()
        self._relErrHist.Divide(self._hist)

    def shiftCentre(self, shift = 0.):
        '''Shift all values by 1, to centre around zero'''
        self._mean  += shift
        self._min   += shift
        self._max   += shift

        for i in range(1, self._nBins+1):
            val = self._hist.GetBinContent(i)
            if val > 0.:
                self._hist.SetBinContent(i, val+shift)
            else:
                # set null points to -666.
                self._hist.SetBinContent(i, -666.)
                # self._hist.SetBinContent(i, 0.)

    def invertHist(self):
        ''' invert an effMap object and all attributes'''

        self._mean  = safe_divide(1., self._mean)
        self._min   = safe_divide(1., self._min)
        self._max   = safe_divide(1., self._max)
                
        for i in range(1, self._hist.GetNbinsX()*
                        self._hist.GetNbinsY()+1000):
            val = self._hist.GetBinContent(i)
            if val > 0.:
                self._hist.SetBinContent(i, 1./val)


###---------------------------------------------------------------------------------------------###

class systMap(object):
    '''Simple systematic plotting'''
    
    ### TO-DO ###
    # 1. implement deltaM plotting
    # 2. top-level pageNum plotting switch?
    # 3. text plotting (add plotString variable) - DONE
    # 4. Add stats print out to each plot - DONE
    # 5. Cap values in each plot (maxVal from run_syst_mapping.py?)
    # 6. implement logging module
    # 7. invert efficiency change for cut_systs
    # 8. cut systs still plotting down (and 3jet syst) - DONE

    def __init__(self, up = None, down = None, central = None, nocuts = None, test = "",
                    model = "", up_noweight = None, down_noweight = None, central_noweight = None,
                    nocuts_noweight = None):

        self._yieldPlots = {
                            "up":       up,
                            "down":     down,
                            "central":  central,
                            "nocuts":   nocuts
                            }

        self._yieldPlots_noweight = {
                                        "up_noweight":       up_noweight,
                                        "down_noweight":     down_noweight,
                                        "central_noweight":  central_noweight,
                                        "nocuts_noweight":   nocuts_noweight
                                        }

        self._effs = {
                        "up":           None,
                        "down":         None,
                        "central":      None,
                        "up_change":    None,
                        "down_change":  None,
                        }

        self._effs_1d = {
                        "up":           None,
                        "down":         None,
                        "central":      None,
                        "up_change":    None,
                        "down_change":  None,
                        }

        self._syst = None
        self._model = model
        self._test = test
        self._cutSyst = True if self._test in ["MHT_MET", "DeadECAL", "LeptonVeto"][:2] else False
        self._ranges = {
                    'x':[],
                    'y':[],
                    'z':[],
        }
        self._logZ = False
        self._plotString = "colz"
        self._nBins = self._yieldPlots['central'].GetNbinsX()* \
                        self._yieldPlots['central'].GetNbinsY() + 1000

        # get model specific plot details
        import plotDetails as pdets
        try:
            self._plotSpec = pdets.modelPlotDetails[self._model]
        except KeyError:
            print ">>> Warning: systMap: Model details not found. Using default values."
            self._plotSpec = {
                    'xRange': [0., 1000.],
                    'yRange': [0., 1000.],
                    'xDMRange': [0., 1000.],
                    'yDMRange': [0., 1000.],
                    'xTitle': "m_{mother} (GeV)",
                    'yTitle': "m_{LSP} (GeV)",
            }

        # get z ranges (varies for each test, not model)
        try:
            self._plotSpec['zRange'] = pdets.systZRanges[self._test]
        except KeyError:
            print ">>> Warning: systMap: Test zRange not found. Using default values."
            self._plotSpec['zRange'] = [0.9,1.1]

        self.makeSystPlots()

    def makeSystPlots(self):
        '''process all yield plots into effs and systs'''

        # zero out the above-diagonal region
        for a in range(1, self._yieldPlots['central'].GetNbinsX()+1):
            xval = self._yieldPlots['central'].GetXaxis().GetBinCenter(a)
            for b in range(1, self._yieldPlots['central'].GetNbinsY()+1):
                yval = self._yieldPlots['central'].GetYaxis().GetBinCenter(b)
                if xval - yval < 0. or a < 0.:
                    bin = self._yieldPlots['central'].FindBin(float(a), float(b))
                    self._yieldPlots['central'].SetBinContent(bin, 0.)
                    self._yieldPlots['up'].SetBinContent(bin, 0.)
                    self._yieldPlots['nocuts'].SetBinContent(bin, 0.)
                    if self._yieldPlots['down']: self._yieldPlots['down'].SetBinContent(bin, 0.)

        # create a load of effMap objects for each variation
        if any(self._yieldPlots_noweight.values()):
            self._effs['central']   = effMap(self._yieldPlots['central'],
                                                self._yieldPlots['nocuts'],
                                                self._yieldPlots_noweight['central_noweight'],
                                                self._yieldPlots_noweight['nocuts_noweight'],
                                                True)
            self._effs['up']        = effMap(self._yieldPlots['up'],
                                                self._yieldPlots['nocuts'],
                                                self._yieldPlots_noweight['up_noweight'],
                                                self._yieldPlots_noweight['nocuts_noweight'],
                                                True)
            self._effs['up_change'] = effMap(self._effs['up']._hist,
                                                self._effs['central']._hist,
                                                self._effs['up']._errHist,
                                                self._effs['central']._errHist)
            if not self._cutSyst:# and self._test != "LeptonVeto":
                # only shift to centre around 0 if it's not a cut systematic
                self._effs['up_change'].shiftCentre(-1)

            for key in ['central', 'up_change']:
                self._effs_1d[key] = self.make_1d_plot(self._effs[key]._hist)

            # if 2-way syst vari exists, do down variation
            if self._yieldPlots['down']:
                self._effs['down']          = effMap(self._yieldPlots['down'],
                                                        self._yieldPlots['nocuts'],
                                                        self._yieldPlots_noweight['down_noweight'],
                                                        self._yieldPlots_noweight['nocuts_noweight'],
                                                        True)
                self._effs['down_change']   = effMap(self._effs['down']._hist,
                                                        self._effs['central']._hist,
                                                        self._effs['down']._errHist,
                                                        self._effs['central']._errHist)
                self._effs['down_change'].shiftCentre(-1)
                self._effs_1d['down_change'] = self.make_1d_plot(self._effs['down_change']._hist)
        else:
            exit("No noweights histograms passed to systMap object.")

        #########################################################
        # here we should have effMap objects for all variations #
        # now go on the create a total systematics map          #
        #########################################################

        if self._cutSyst:
            ### cut hist scenario ###
            # invert up_change to represent a cut efficiency
            # FIXME: THIS MAY NOT WORK BY CHANGING UP HERE
            self._effs['up_change'].invertHist()

        # if self._test == "LeptonVeto":
        #   self._effs['up_change'].shiftCentre(-1)

        # calculate overall syst value
        tmp_hist    = self._effs['central']._hist.Clone() # cheeky copy of some old histo
        tmp_errHist = self._effs['central']._errHist.Clone()

        for i in range(1, self._nBins+1):
            
            # skip null points
            if abs(self._effs['up_change']._hist.GetBinContent(i)) == 666.:
                tmp_hist.SetBinContent(i, -666.)
                # tmp_hist.SetBinContent(i, 0.)
                continue

            # get all systematic values
            pos_syst = abs(self._effs['up_change']._hist.GetBinContent(i))
            pos_err = abs(self._effs['up_change']._errHist.GetBinContent(i))
            if self._effs['down_change']:
                neg_syst = abs(self._effs['down_change']._hist.GetBinContent(i))
                neg_err = abs(self._effs['down_change']._errHist.GetBinContent(i))
            else:
                neg_syst = -1.
                neg_err = 0.000001

            
            # pick the most significant systematic, based on requirements in syst_picker()
            this_syst, this_err = syst_picker( (pos_syst, pos_err), (neg_syst, neg_err) )

            # if pos_syst >= neg_syst:
            #     this_syst = pos_syst
            #     this_err = self._effs['up_change']._errHist.GetBinContent(i)
            #     # print " ", this_syst, this_err, safe_divide(this_err, this_syst)
            #     # print neg_syst, neg_err, safe_divide(neg_err, neg_syst)
            # else:
            #     this_syst = neg_syst
            #     this_err = self._effs['down_change']._errHist.GetBinContent(i)
            #     # print " ", this_syst, this_err, safe_divide(this_err, this_syst)
            #     # print pos_syst, pos_err, safe_divide(pos_err, pos_syst)

            tmp_hist.SetBinContent(i, this_syst)
            tmp_errHist.SetBinContent(i, this_err)
        
        # squash the outliers, if there are any
        tmp_hist = self.squash_outliers(tmp_hist)

        if self._model not in ["T2cc", "T2_4body"]:
            # fill any holes, so every point has a systematic
            tmp_hist = self.fill_holes(tmp_hist)

        # create effMap from new total syst hist
        self._syst = effMap(nom = tmp_hist, nom_err = tmp_errHist)

        self._syst_1d = self.make_1d_plot(self._syst._hist)

        # now need to smooth this final map to account for fluctuations
        # self.syst_smooth(self._syst)

        del tmp_hist, tmp_errHist

    def print_all(self, label = "", plotText = False):
        '''print syst output plots'''

        if plotText:
            if self._model not in ["T2cc", "T2_4body"]:
                print ">>> Warning: systMap: print_all: Text plot will look shit for this model,"
                print "    so it will be turned off. Remember, looks are everything in HEP."
                self._plotString = "colz"
            else:
                self._plotString = "colztext"
                r.gStyle.SetPaintTextFormat("0.7f");

        # create instance of a multi page PDF
        pdf0 = multiPagePDF(outFileName = "out/%s_%s_%s.pdf" % (self._model, self._test, label),
                            title = "Systematics - %s - %s" % (self._test, label),
                            title_page = False)

        # setup fine grain z-axis
        sutils.set_palette()

        if not self._cutSyst:
            # draw total systematic (don't plot for cut systematics)
            self.draw_plot(self._syst,
                            pdf0,
                            "%s Systematic - %s" % (self._test, label), 
                            "Systematic Value",
                            0.)
            self.draw_plot(self._syst,
                            pdf0,
                            "%s Systematic Error - %s" % (self._test, label), 
                            "Systematic Value",
                            0., err=True)

            # also draw relative err value - for debugging purproses
            # self._syst._relErrHist.Draw("colztext")
            # self._syst._relErrHist.SetTitle("Relative Error")
            # pdf0.AddPage()

            pdf0.setLogY(1)
            self._syst_1d.Draw("hist")
            pdf0.AddPage()
            pdf0.setLogY(0)

        # draw each variation
        for key in ["central", "up", "up_change", "down", "down_change"]:
            if "down" not in key or self._effs['down']:
                
                if self._cutSyst and key == "up_change":
                    key_string = "acceptance"
                else:
                    key_string = key

                self.draw_plot(self._effs[key],
                                pdf0,
                                "Efficiency %s %s - %s" % (self._test,
                                    key_string, label),
                                "Acceptance",
                                1. if "change" in key else 0.)
                self.draw_plot(self._effs[key],
                              pdf0,
                              "Efficiency %s %s Error - %s" % (self._test, key_string, label),
                              "Acceptance",
                              1. if "change" in key else 0., err=True)
                if self._effs_1d[key]:

                    pdf0.setLogY(1)
                    r.gStyle.SetOptStat(111111)
                    self._effs_1d[key].Draw("hist")
                    pdf0.AddPage()
                    pdf0.setLogY(0)
                    r.gStyle.SetOptStat(0)
        
        pdf0.close()

    def squash_outliers(self, hist = None):
        vals = []
        points = []
        for i in range(1, self._nBins + 1):
            val = hist.GetBinContent(i)
            centers = sutils.get_bin_centre_vals(hist, i)
            if centers['x'] > centers['y']:
                if val > 0.:
                    vals.append(val)
                    points.append(centers)

        to_squash = sutils.reject_outliers(vals, 3.)

        print "squash length:", len(to_squash)

        for sq in to_squash:
            # new_val = np.median(vals)*2. #new val to replace outliers

            # try replacing with the median, rather than x2
            new_val = np.median(vals) #new val to replace outliers - very non-conservative?

            point = points[sq]
            bin = hist.FindBin(point['x'], point['y'])
            hist.SetBinContent(bin, new_val)

        return hist

    def fill_holes(self, hist = None):

        for xbin in range(1, hist.GetNbinsX()+1):
            for ybin in range(1, hist.GetNbinsY()+1):
                centers = {'x': hist.GetXaxis().GetBinCenter(xbin),
                            'y': hist.GetYaxis().GetBinCenter(ybin)}
                
                # only consider points in the scan
                if not pdets.point_white_list(self._model)(centers['x'], centers['y']): continue

                val = hist.GetBinContent(xbin, ybin)
                
                if val > 0.: continue

                to_avg = []
                # get swiss cross values to avg over
                for xtmp in [-1, 1]:
                    tmpval = hist.GetBinContent(xbin+xtmp, ybin)
                    if tmpval > 0.:
                        to_avg.append(tmpval)

                for ytmp in [-1, 1]:
                    tmpval = hist.GetBinContent(xbin, ybin+ytmp)
                    if tmpval > 0.:
                        to_avg.append(tmpval)

                if len(to_avg) > 0:
                    avgval = np.average(to_avg)
                    hist.SetBinContent(xbin, ybin, avgval)

        return hist

    def make_1d_plot(self, hist = None):
        
        outhist = r.TH1D("1d", "1d", 200, 0, 0.6)

        for i in range(1, self._nBins+1):
            centers = sutils.get_bin_centre_vals(hist, i)
            # only plot values below diagonal
            if centers["y"] > centers["x"]: 
                continue
            val = hist.GetBinContent(i)
            if val>0.6:
                outhist.SetBinContent(200, 1.)
            else:
                outhist.Fill(val)

        return outhist

    def draw_plot(self, effMap = None, pdfFile = None, title = "", zTitle = "", shiftZ = 0., err = False):
        '''draw a single plot on a single page of pdfFile'''

        if not err:
            hist = effMap._hist.Clone()
        else:
            hist = effMap._errHist.Clone()
        hist.Draw(self._plotString)
        hist.SetTitle(title)

        if "change" in title:
            hist.GetZaxis().SetTitle(zTitle + " change")
        else:
            hist.GetZaxis().SetTitle(zTitle)

        # if err:
        #     for i in range(1, effMap._nBins+10):
        #         if effMap._errHist.GetBinContent(i)>0.:
        #             xbin, ybin, zbin = r.Long(0.), r.Long(0.), r.Long(0.)
        #             effMap._errHist.GetBinXYZ(i, xbin, ybin, zbin)
        #             print i, xbin, ybin, effMap._errHist.GetBinContent(i)

        if self._cutSyst:
            # remove the z-axis shift for cut-based systematics
            shiftZ = 0. 
        self.setDetails(hist, shiftZ)

        #hack to fix systematic zRange
        if "Syst" in title and not self._cutSyst:
            hist.GetZaxis().SetRangeUser(0., self._plotSpec['zRange'][1]-1.)

        if "text" in self._plotString:
            hist.SetMarkerSize(0.8)

        # draw all the stats numbers
        num0 = r.TLatex(0.151,0.8,"#scale[0.6]{avg = %.4f}" % effMap._mean)
        num0.SetNDC()
        num0.Draw("same")

        num1 = r.TLatex(0.15,0.77,"#scale[0.6]{RMS= %.4f}" % effMap._rms)
        num1.SetNDC()
        num1.Draw("same")

        num2 = r.TLatex(0.15,0.74,"#scale[0.6]{min = %.4f}" % effMap._min)
        num2.SetNDC()
        num2.Draw("same")

        num3 = r.TLatex(0.15,0.71,"#scale[0.6]{max = %.4f}" % effMap._max)
        num3.SetNDC()
        num3.Draw("same")

        pdfFile.AddPage()

    def setDetails(self, hist=None, zoffset=0.):
        '''set the ranges, titles, and title text specs of a plot'''

        hist.GetXaxis().SetRangeUser(self._plotSpec['xRange'][0], self._plotSpec['xRange'][1])
        hist.GetXaxis().SetTitle(self._plotSpec['xTitle'])
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetXaxis().SetTitleSize(0.04)

        hist.GetYaxis().SetRangeUser(self._plotSpec['yRange'][0], self._plotSpec['yRange'][1])
        hist.GetYaxis().SetTitle(self._plotSpec['yTitle'])
        hist.GetYaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetTitleSize(0.04)

        if zoffset:
            hist.GetZaxis().SetRangeUser(self._plotSpec['zRange'][0]-zoffset,
                                            self._plotSpec['zRange'][1]-zoffset)
        # if self._cutSyst:
            # hist.GetZaxis().SetRangeUser(self., 1.2)
        if self._cutSyst:
            hist.SetMaximum(1.)

        hist.GetZaxis().SetTitleOffset(0.95)
        hist.GetZaxis().SetTitleSize(0.05)

    def syst_smooth(self, emap = None, iterations = 5):
        """smooth the emap"""

        new_hist = emap._hist.Clone()

        n = 0
        # note: ONLY WORKS FOR T2CC AND T24BODY because of binning assumptions
        while n < iterations:
            print "smooth", n
            for xbin in range(1, emap._hist.GetNbinsX()+1):
                for ybin in range(1, emap._hist.GetNbinsY()+1):
                    val = emap._hist.GetBinContent(xbin, ybin)
                    
                    if val <= 0.: continue
                    
                    vals = []
                    errs = []

                    for xtmp in range(-2,3):
                        # print xtmp,xtmp*5
                        tmp_val = emap._hist.GetBinContent( xbin+xtmp, ybin+xtmp*5)
                        if tmp_val > 0.:
                            tmp_err = float(emap._errHist.GetBinContent( xbin+xtmp, ybin+xtmp*5)/tmp_val)
                            vals.append(tmp_val)

                            errs.append(float(1./ma.pow(tmp_err, 2)))
                    # get weighted average considering errs
                    err_sum = np.sum(errs)
                    for j in range(len(errs)):
                        errs[j]/=err_sum

                    # ave_val = np.average(vals, weights=errs)
                    ave_val = np.average(vals)
                    new_hist.SetBinContent(xbin, ybin, ave_val)
            n+=1

        emap._hist = new_hist.Clone()

    def __del__(self):
        '''destructor to deal with effMap objects'''

        del self._effs
        del self._yieldPlots


###---------------------------------------------------------------------------------------------###

class multiPagePDF(object):
    '''class for producing multipage PDF file'''

    def __init__(self, outFileName = "", title = "_NoTitle_", title_page = True):
        self._canv = r.TCanvas() # FIXME: pick good size
        self._fName = outFileName
        self._pageNums = True
        self._doName = True
        self._analyst = "Chris Lucas"
        self._title = title
        self._pageCtr = 1 #initiate page numbers
        self._title_page = title_page
        if self._title_page:
            self.makeTitlePage()

    def setLogY(self, val = None):
        self._canv.SetLogy(val)

    def makeTitlePage(self):
        '''make a title page with title, name, date'''
        from time import strftime

        title_text = r.TLatex(0.07,0.8,"#scale[1.5]{%s}" % self._title)
        title_text.SetNDC()
        title_text.Draw("same")

        if self._doName:
            name_text = r.TLatex(0.07, 0.2, "Analyst: %s" % self._analyst)
            name_text.SetNDC()
            name_text.Draw("same")

        date_text = r.TLatex(0.07, 0.15, "#scale[0.7]{%s}" %
                                strftime("%a, %d %b %Y %H:%M:%S")) # FIXME: make this text smaller
        date_text.SetNDC()
        date_text.Draw("same")

        self._canv.Print(self._fName+"(")
        r.gPad.SetRightMargin(0.15)
        r.gPad.SetLeftMargin(0.10)
        r.gPad.SetTopMargin(0.08)
        r.gPad.SetBottomMargin(0.12)

    def AddPage(self):
        '''Add page with canvas drawn and page number'''
        num = r.TLatex(0.97,0.025,"%d"%(self._pageCtr))
        num.SetNDC()
        if self._pageNums: num.Draw("same")
        if not self._title_page and self._pageCtr == 1:
            self._canv.Print(self._fName + "(")
        else:
            self._canv.Print(self._fName)
        self._pageCtr += 1
        pass

    def close(self):
        '''close the pdf file'''
        self._canv.Print(self._fName+"]")
        pass
