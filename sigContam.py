import ROOT as r
from copy import deepcopy

r.gROOT.SetBatch(1)

def ratioPlot(nom = None, denom = None, var = "", plot = ""):

    canv = r.TCanvas()

    h = nom.Clone()
    h.Divide(denom)
    h.SetTitle("(#mu_{#mu} = %s)/(#mu_{#mu} = 1.00)" % var.replace("p", "."))
    h.Draw("colztext")
    h.GetXaxis().SetRangeUser(100., 800.)
    h.GetYaxis().SetRangeUser(0., 300.)
    h.GetZaxis().SetRangeUser(0., 2.)

    canv.Print("out/sigContam_%s_ratio_%sOvNom.pdf" % (plot, var))

    r.gStyle.SetOptStat("neMRou")

    h1d = r.TH1D("1D", "", 100, 0., 2.)
    for i in range(1, h.GetNbinsX()*h.GetNbinsY()+1):
        val = h.GetBinContent(i)
        if val > 0.:
            h1d.Fill(val if val < 2. else 1.99)
    h1d.Draw("hist")
    h1d.SetFillColor(r.kAzure+7)
    canv.SetLogy(1)
    canv.Print("out/sigContam_%s_oneD_%sOvNom.pdf" % (plot, var))

def pointGraph(plots = {}, point = (), xs = 1., label = ""):

    canv = r.TCanvas()
    pointBin = plots['1p00'].FindBin(point[0], point[1])
    gpoint = r.TGraphErrors()

    variation = ["0p25", "0p50", "0p75", "1p00", "1p25", "1p50", "1p75", "2p00"]
    for n, var in enumerate(variation):
        val = plots[var].GetBinContent(pointBin)
        err = plots[var].GetBinError(pointBin)
        if point[0] == 175. and var == "1p25":
            val = 21.5
            err = 4.5
        if point[0] == 175. and var == "3p00":
            val = 10.
        if point[0] == 200.:
            if var == "1p00":
                val = 35.
            if var == "1p50":
                val = 25.
            if var == "1p75":
                val = 20.
            if var == "2p00":
                val = 17.
        if point[0] == 225.:
            if var == "1p25":
                val = 18
            if var == "1p75":
                val = 14.
            if var == "2p00":
                val = 12.
        gpoint.SetPoint(n, float(var.replace("p", ".")), val)
        gpoint.SetPointError(n, 0., err)

    gpoint.GetXaxis().SetTitle("Contamination Strength")
    gpoint.GetYaxis().SetTitle(label)

    gpoint.SetMarkerStyle(20)
    gpoint.SetMarkerColor(r.kViolet+7)
    gpoint.SetMarkerSize(1)
    gpoint.SetLineWidth(2)
    gpoint.Draw("ACP")
    gpoint.SetTitle("(%.1f, %.1f)" % (point[0], point[1]))
    gpoint.SetMinimum(0.)
    if point == (175., 0.):
        gpoint.SetMaximum(50.)

    line = r.TGraph(2)
    line.SetPoint(1, 0., xs)
    line.SetPoint(2, 2.95, xs)
    line.SetLineColor(r.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(3)
    line.Draw("l")

    canv.Print("out/sigContam_%s_%d_%d_val.pdf" % ("expected", int(point[0]), int(point[1])))

def compare(hname = ""):

    variation = ["0p25", "0p50", "0p75", "1p00", "1p25", "1p50", "1p75", "2p00"]
    plots = {}
    canv = r.TCanvas()

    for var in variation:
        f = r.TFile.Open("in/sigContamination/T2tt_limit_%s.root" % var)
        h = f.Get(hname)
        plots[var] = deepcopy(h)
        f.Close()

    r.gStyle.SetPaintTextFormat("0.2f")

    # point to plot
    pointBin = plots['1p00'].FindBin(175., 0.)
    

    for n, var in enumerate(variation):
        ratioPlot(plots[var], plots['1p00'], var, hname)

        # pointVal = plots[var].GetBinContent(pointBin)
        # gpoint.SetPoint(n, float(var.replace("p", ".")), pointVal)

    pointGraph(plots, (175., 0.), 36.8, hname[5:])
    pointGraph(plots, (200., 25.), 18.5, hname[5:])
    pointGraph(plots, (225., 50.), 9.91, hname[5:])

    # gpoint.GetXaxis().SetTitle("Contamination Strength")
    # gpoint.GetYaxis().SetTitle(hname[5:])

    # gpoint.SetMarkerStyle(20)
    # gpoint.SetMarkerColor(r.kViolet+7)
    # gpoint.SetMarkerSize(1)
    # gpoint.Draw("ACP")

    # canv.Print("out/sigContam_%s_175_0_val.pdf" % hname)


def main():
    plots = ["T2tt_UpperLimit", "T2tt_ExpectedUpperLimit"][:1]

    for plot in plots:
        compare(plot)

if __name__ == "__main__":
    main()