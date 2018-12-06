import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend


file40 = TFile("50p/upper_limits.root")
file55 = TFile("m55_50p/upper_limits.root")
file15 = TFile("m15_50p/upper_limits.root")
    
g40 = file40.Get("Graph")
g55 = file55.Get("Graph")
g15 = file15.Get("Graph")

g40.SetLineColor(2)
g55.SetLineColor(3)
g15.SetLineColor(4)

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
frame.SetMaximum(5)
frame.GetXaxis().SetLimits(1,1000)

g40.Draw('L')
g55.Draw('Lsame')
g15.Draw('Lsame')

ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLogx()
frame.Draw('sameaxis')


x1 = 0.2
x2 = x1 + 0.24
y2 = 0.86
y1 = 0.6
legend = TLegend(x1,y1,x2,y2)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.041)
legend.SetTextFont(42)
legend.AddEntry(g15, "m_{S}=15 GeV, dummy systs = 50%",'L')
legend.AddEntry(g40, "m_{S}=40 GeV, dummy systs = 50%",'L')
legend.AddEntry(g55, "m_{S}=55 GeV, dummy systs = 50%",'L')


legend.Draw()

c.SaveAs("overlay_m.png")
c.SaveAs("overlay_m.pdf")
