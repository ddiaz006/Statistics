import ROOT as rt

############################
# from javier duarte
############################
procNorm = rt.ProcessNormalization('norm','norm',1)
x = rt.RooRealVar('x', 'x', 0, -5, 5)
kappaLo = 1.1
kappaHi = 1.2
procNorm.addAsymmLogNormal(kappaLo, kappaHi, x)
frame = x.frame()
procNorm.plotOn(frame)
c= rt.TCanvas('c','c',500,400)
frame.Draw()
c.Print('procNorm.pdf')

procNorm = rt.ProcessNormalization('norm','norm',1)
x = rt.RooRealVar('x', 'x', 0, -5, 5)
y = rt.RooRealVar('y', 'y', 0, -5, 5)
kappaLoX = 1.1
kappaHiX = 1.2
kappaLoY = 0.9
kappaHiY = 1.1
procNorm.addAsymmLogNormal(kappaLoX, kappaHiX, x)
procNorm.addAsymmLogNormal(kappaLoY, kappaHiY, y)

c= rt.TCanvas('c','c',500,400)
frame = x.frame()
procNorm.plotOn(frame)
frame.Draw()
c.Print('procNorm_x.pdf')
frame = y.frame()
procNorm.plotOn(frame)
frame.Draw()
c.Print('procNorm_y.pdf')

##########################
# ben manual symmetrized
##########################

form = rt.RooFormulaVar("form", "form", "TMath::Power(1+.2,TMath::Abs(@0))", rt.RooArgList(x))



