import sys
sys.path.append('src')

from Experiment import Experiment
from Distribution import Distribution

import random
import os

import ROOT

random.seed(919)


########################################
# Main
########################################

def main():

    # # ======================================
    # # Set up ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
    # ROOT.gStyle.SetOptFit(1011)
    ROOT.gStyle.SetOptStat(0)

    ROOT.gSystem.Load(
        '/mnt/t3nfs01/data01/shome/tklijnsm/Combination/RooUnfold-1.1.1/libRooUnfold'
        )



    Exp = Experiment()


    Exp.SmearSigma = 0.6
    Exp.SetPlotPrefix( 'SmallSmear_' )
    Exp.GeneratePseudoData( nEvents=10000, setOutput=True )

    # Exp.BinByBinUnfoldTest()

    # Exp.InvertedKUnfoldTest()

    # Exp.BayesUnfoldTest()

    Exp.SVDUnfoldTest()

    Exp.SVDUnfoldTest_TUnfold()




    Exp.SmearSigma = 0.9
    Exp.SetPlotPrefix( 'BigSmear_' )
    Exp.GeneratePseudoData( nEvents=1000, setOutput=True )

    # Exp.BinByBinUnfoldTest()

    # Exp.InvertedKUnfoldTest()

    # Exp.BayesUnfoldTest()

    Exp.SVDUnfoldTest()

    Exp.SVDUnfoldTest_TUnfold()







########################################
# Functions
########################################




########################################
# End of Main
########################################
if __name__ == "__main__":
    main()