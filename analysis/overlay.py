import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend


file10p = TFile("10p/upper_limits.root")
file50p = TFile("50p/upper_limits.root")
file100p = TFile("100p/upper_limits.root")
file300p = TFile("300p/upper_limits.root")
    
g10p = file10p.Get("Graph")
g50p = file50p.Get("Graph")
g100p = file100p.Get("Graph")
g300p = file300p.Get("Graph")

g50p.SetLineColor(2)
g100p.SetLineColor(3)
g300p.SetLineColor(4)

W = 800
H  = 600
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.04*W
c = TCanvas("c","c",100,100,W,H)
c.SetFillColor(0)
c.SetBorderMode(0)
c.SetFrameFillStyle(0)
c.SetFrameBorderMode(0)
c.SetLeftMargin( L/W )
c.SetRightMargin( R/W )
c.SetTopMargin( T/H )
c.SetBottomMargin( B/H )
c.SetTickx(0)
c.SetTicky(0)
c.SetGrid()
c.cd()
frame = c.DrawFrame(1.4,0.001, 4.1, 10)
frame.GetYaxis().CenterTitle()
frame.GetYaxis().SetTitleSize(0.05)
frame.GetXaxis().SetTitleSize(0.05)
frame.GetXaxis().SetLabelSize(0.04)
frame.GetYaxis().SetLabelSize(0.04)
frame.GetYaxis().SetTitleOffset(0.9)
frame.GetXaxis().SetNdivisions(508)
frame.GetYaxis().CenterTitle(True)
frame.GetYaxis().SetTitle("95% upper limit on BR")
frame.GetXaxis().SetTitle("c#tau [mm]")
frame.SetMinimum(0)
frame.SetMaximum(2.5)
frame.GetXaxis().SetLimits(1,1000)

g10p.Draw('L')
g50p.Draw('Lsame')
g100p.Draw('Lsame')
g300p.Draw('Lsame')

ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLogx()
frame.Draw('sameaxis')


x1 = 0.4
x2 = x1 + 0.24
y2 = 0.76
y1 = 0.5
legend = TLegend(x1,y1,x2,y2)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.041)
legend.SetTextFont(42)
legend.AddEntry(g300p, "dummy systs = 300%",'L')
legend.AddEntry(g100p, "dummy systs = 100%",'L')
legend.AddEntry(g50p, "dummy systs = 50%",'L')
legend.AddEntry(g10p, "dummy systs = 10%",'L')

legend.Draw()

c.SaveAs("overlay.png")
c.SaveAs("overlay.pdf")
