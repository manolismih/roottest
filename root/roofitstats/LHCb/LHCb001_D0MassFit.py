# RooFit example of D0 mass fit with toy data.

from ROOT import RooWorkspace, RooArgSet, RooDataSet, RooFit, RooAbsData, TCanvas, TColor, TTree, TPad, gStyle, TLine, gSystem, kTRUE

def main():
    ws = RooWorkspace()
    
    import_model(ws)
    model = ws.pdf('model')
    data = generate_data(ws, model)
    do_mass_fit(ws, data)

def import_model(workspace):
    print('Make D0 model.')

    # Fit parameters
    workspace.factory('prod::sigma(frac[0.5,0,1],s[10,8,12])')
    workspace.factory('f[1.7, 0, 3]')
    workspace.factory('prod::sigma1(sigma,f)')
    workspace.factory('a1[0, -0.2, 0.2]')
    workspace.factory('a2[0, -0.2, 0.2]')
    workspace.factory('alpha[3]')

    workspace.factory('Gaussian::sig2(D0_LoKi_DTF_M[1805,1925],mean[1865, 1850, 1870],sigma1)')
    workspace.factory('CBShape::sig1(D0_LoKi_DTF_M,mean, sigma, alpha, n[3])')

    workspace.factory('Chebychev::bkg(D0_LoKi_DTF_M,{a1,a2})')

    workspace.factory('SUM::sig(frac*sig1,sig2)')
    workspace.factory('SUM::model(nsig[10e6, 5e6, 15e6]*sig,nbkg[5e6, 2e6, 15e6]*bkg)')
    workspace.Print()
    print('Making model finished.')

def generate_data(workspace, pdf):
    # Get what we need from workspace
    D0_LoKi_DTF_M = workspace.var('D0_LoKi_DTF_M')
    myset = RooArgSet(D0_LoKi_DTF_M)

    # Generate toy events
    dataset = pdf.generate(myset, 15000000)
    return dataset

def do_mass_fit(workspace, dataset):
    print('Do fit.')

    # Get what we need from workspace
    D0_LoKi_DTF_M = workspace.var('D0_LoKi_DTF_M')

    model = workspace.pdf('model')
    bkg = workspace.pdf('bkg')
    sig = workspace.pdf('sig')
    sig1 = workspace.pdf('sig1')
    sig2 = workspace.pdf('sig2')

    frame = D0_LoKi_DTF_M.frame()
    dataset.plotOn(frame, RooFit.DataError(RooAbsData.SumW2))
    fitresult = model.fitTo(dataset, RooFit.Save(), RooFit.Extended())
    model.plotOn(frame)

    ndof = fitresult.floatParsFinal().getSize()
    print('Fit has chi2/ndof: ' + str(frame.chiSquare(ndof)))

    # Calculate pulls
    frame_p = D0_LoKi_DTF_M.frame(RooFit.Title(''))
    frame_p.addPlotable(frame.pullHist(), 'P')

    model.plotOn(frame, RooFit.Components('bkg'), RooFit.LineColor(419), RooFit.LineStyle(5))
    model.plotOn(frame, RooFit.Components('sig1'), RooFit.LineColor(874), RooFit.LineStyle(2))
    model.plotOn(frame, RooFit.Components('sig2'), RooFit.LineColor(886), RooFit.LineStyle(4))
    model.paramOn(frame, RooFit.Format('NEU', RooFit.AutoPrecision(1)), RooFit.Layout(0.7,0.9,0.9))
    frame.getAttText().SetTextSize(0.03)
    frame.getAttText().SetTextFont(132)
    frame.GetXaxis().SetLabelOffset(999) # Remove X axis labels
    frame.GetYaxis().SetLabelFont(132)
    frame.GetYaxis().SetTitleFont(132)

    # Pull frame
    frame_p.SetMinimum(-5)
    frame_p.SetMaximum(5)
    frame_p.GetXaxis().SetTitle('#font[12]{m}(#font[12]{D^{0}} ) [MeV]')
    frame_p.GetXaxis().SetTitleOffset(1.)
    frame_p.GetXaxis().SetLabelFont(132)
    frame_p.GetXaxis().SetTitleFont(132)
    frame_p.SetLabelSize(0.06, 'XY')
    frame_p.GetYaxis().SetTitle('(data - model)/#font[12]{#sigma}_{data}')
    frame_p.GetYaxis().SetTitleOffset(0.5)
    frame_p.GetYaxis().SetNdivisions(209)
    frame_p.GetYaxis().SetLabelFont(132)
    frame_p.GetYaxis().SetTitleFont(132)
    frame_p.SetTitleSize(0.06, 'XY')

    # Create Canvas and Pads and Draw plot
    c1 = TCanvas('', '', 600, 600)
    upperPad = TPad('upperPad', 'upperPad', .005, .3475, .995, .995)
    lowerPad = TPad('lowerPad', 'lowerPad', .005, .005, .995, .3475)
    upperPad.SetBottomMargin(0.03)
    lowerPad.SetTopMargin(0.02)
    lowerPad.SetBottomMargin(0.3)

    upperPad.Draw()
    lowerPad.Draw()
    upperPad.cd()
    frame.Draw()
    lowerPad.cd()
    gStyle.SetOptTitle(0)
    line = TLine(1805,0,1925,0)
    line.SetLineColor(1)
    frame_p.Draw()
    line.Draw()
    c1.SaveAs('RooFit_example.pdf')


if __name__ == '__main__':
    main()
