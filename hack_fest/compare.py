import ROOT as r

r.gStyle.SetOptStat(0)

def threeToTwo(h3) :
    name = h3.GetName()
    binsz = h3.GetNbinsZ()

    h2 = r.TH2D(name+"_2D",h3.GetTitle(),
                h3.GetNbinsX(), h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax(),
                h3.GetNbinsY(), h3.GetYaxis().GetXmin(), h3.GetYaxis().GetXmax(),
                )
                
    for iX in range(1, 1+h3.GetNbinsX()) :
        for iY in range(1, 1+h3.GetNbinsY()) :
            content = h3.GetBinContent(iX, iY, 1) + h3.GetBinContent(iX, iY, 2)+ h3.GetBinContent(iX, iY, 0)
            h2.SetBinContent(iX, iY, content)
    h2.GetZaxis().SetTitle(h3.GetZaxis().GetTitle())
    
    return h2

def GetHist(File = None, folder = None, hist = None, Norm = None, rebinX = None, rebinY = None):
    h = None
    for f in folder:
        directory = File.Get(f)
        a = directory.Get(hist)
        if h is None:
            h = a.Clone()
        else: h.Add(a)

    if rebinX:
        h.RebinX(rebinX)
    if rebinY:
        h.RebinY(rebinY)

    return h

def make_comparison_plot(h1 = None, h2 = None, fname = ""):
    canv = r.TCanvas()
    
    lg = r.TLegend(0.6, 0.75, 0.89, 0.89)

    h1.Draw("hist")
    h1.RebinX(2)
    h1.SetLineColor(r.kBlue)
    lg.AddEntry(h1, "T2cc", "L")

    h2.Draw("histsame")
    h2.RebinX(2)
    h2.SetLineColor(r.kRed)
    lg.AddEntry(h2, "T2_4body", "L")

    lg.Draw()
    lg.SetFillColor(0)
    lg.SetLineColor(0)
    canv.SetLogy(1)
    canv.Print(fname)

    div = h1.Clone()
    div.Divide(h2)
    div.Draw("")
    div.SetTitle("T2cc/T2_4body")
    div.SetMaximum(3)
    canv.SetGridy(1)
    canv.SetLogy(0)
    canv.Print(fname.replace(".pdf", "_ratio.pdf"))


def main():

    hists = ["ISRPt0", "ISRPt1", "ISRSumPt", "StopPt0", "StopPt1", "StopSumPt"]
    models = ["T2cc", "T2_4body"]
    threshs = ["73.3", "86.7", "100.0"][-1:]
    rw_proc = ["theirs", "ours", "none"][0]

    in_hists = {}
    # for model in models:
    #     in_hists[model] = {}
    #     file = r.TFile.Open("in/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" % (model))
    #     for hname in hists:
    #         in_hists[model][hname] = GetHist(File = file,folder = ["smsScan_before",],hist = hname, Norm = None, rebinY=None)

    for hname in hists:
        in_hists[hname] = {}
        for model in models:
            if model == "T2cc":
                file = r.TFile.Open("in/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0.root" % (model))
            else:
                file = r.TFile.Open("in/sigScan_%s_had_2012_100.0_bt0.0_MChi-1.0_%s.root" % (model, rw_proc))
            in_hists[hname][model] = GetHist(File = file,folder = ["smsScan_before",],hist = hname, Norm = None, rebinY=None)            
        make_comparison_plot(in_hists[hname]["T2cc"], in_hists[hname]["T2_4body"], "out_%s_%s.pdf" % (rw_proc, hname))
    


if __name__ == "__main__":
    main()