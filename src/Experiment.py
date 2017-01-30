#!/usr/bin/env python
"""
Thomas:
"""

########################################
# Imports
########################################

from Distribution import Distribution

import os
import random

from array import array
from math import sqrt, exp, log, pi
from scipy.special import gammaln
from copy import deepcopy
import numpy


import ROOT

# Function that returns a unique root name (if the name of an object is not important anyway)
rootc = 1000
def GetUniqueRootName():
    global rootc
    rootc += 1
    return 'obj{0}'.format(rootc-1)


class Experiment:

    # ======================================
    # Initialization

    def __init__( self ):
        
        self.name = 'test'


        self.plotdir = 'plots/'

        self.c = ROOT.TCanvas( 'c', 'c', 800, 600 )

        self.nBins = 20

        self.SmearSigma = 0.6


        self.Persistence = []

        self.PlotPrefix = False


    def SetPlotPrefix( self, prefix ):
        self.PlotPrefix = True
        self.PlotPrefixText = prefix


    def Smear( self, x ):
        # Smear function

        # Smear = lambda x: x + random.gauss( 0., 0.6 )
        # return x + random.gauss( 0., 0.6 )
        return x + random.gauss( 0., self.SmearSigma )


    def GeneratePseudoData( self, nEvents, setOutput=False ):
        nBins = self.nBins

        # Define xAxis
        xAxis, xAxisCenters = self.MakeAxis( -6., 6., nBins )
        self.xAxis = xAxis
        self.xAxisCenters = xAxisCenters

        # ======================================
        # Make truth spectrum

        # Simple Gaussian function
        Gauss =   lambda x, mu, sigma: 1./sqrt(2*pi*sigma**2) * exp( -0.5* ((x-mu)/sigma)**2 )
        Poisson = lambda k, labda:     exp(  k*log(labda) - labda - gammaln(k+1)  )

        # Central truth distribution function
        truth = lambda x: 0.8 * Gauss( x, 2., 1. ) + 0.2 * Gauss( x, -2., 1. )

        truthAxis = [ truth(x) for x in xAxisCenters ]
        binWidth = float(xAxis[1] - xAxis[0])


        # ======================================
        # Sample truth events

        truthDistribution = Distribution( xAxisCenters, truthAxis )

        sampledTruthAxis = [ 0 for i in xrange(nBins) ]
        smearedTruthAxis = [ 0 for i in xrange(nBins) ]

        events = []

        for iEvent in xrange(nEvents):
            
            iBin, xDrawn, yDrawn = truthDistribution.DrawSample()

            # Get the smeared x value
            xSmeared = self.Smear(xDrawn)

            # Find corresponding bin
            iBinSmeared = self.GetClosestIndex( xSmeared, xAxisCenters )

            # print 'Original bin: {0:4} | Smeared bin: {1:4}'.format( iBin, iBinSmeared )

            sampledTruthAxis[iBin] += 1
            smearedTruthAxis[iBinSmeared] += 1

            events.append(  ( xDrawn, xSmeared )  )


        # Record the pure count lists
        sampledCountAxis = sampledTruthAxis[:]
        smearedCountAxis = smearedTruthAxis[:]

        # Normalize the lists
        sumSampledTruthAxis = float(sum(sampledTruthAxis))
        sampledTruthAxis = [ y/(sumSampledTruthAxis*binWidth) for y in sampledTruthAxis ]

        sumSmearedTruthAxis = float(sum(smearedTruthAxis))
        smearedTruthAxis = [ y/(sumSmearedTruthAxis*binWidth) for y in smearedTruthAxis ]


        if setOutput:
            # Set the output variables as class variables
            self.xAxis = xAxis
            self.truthAxis = truthAxis
            self.sampledTruthAxis = sampledTruthAxis
            self.smearedTruthAxis = smearedTruthAxis
            self.sampledCountAxis = sampledCountAxis
            self.smearedCountAxis = smearedCountAxis
            self.events = events
        else:
            # Return the relevant lists
            return ( xAxis, truthAxis, sampledTruthAxis, smearedTruthAxis )



    def PlotBasicHists( self ):

        self.c.Clear()


        # -------------------
        # Truth plot

        Htruth = self.HistogramFromAxes( self.xAxis, self.truthAxis, normalize=True )
        Htruth.Draw()
        Htruth.SetLineColor(1)
        Htruth.SetLineStyle(2)
        Htruth.SetLineWidth(3)
        Htruth.SetName('truth')
        Htruth.SetMaximum( 1.5*max(self.truthAxis) )

        self.Save( self.c, 'truth' )

        # Set pointer to Htruth in class so the axes can be edited
        self.base = Htruth


        # -------------------
        # Sampled truth plot

        Hsampled = self.HistogramFromAxes( self.xAxis, self.sampledTruthAxis )
        Hsmeared = self.HistogramFromAxes( self.xAxis, self.smearedTruthAxis )

        Hsampled.Draw('SAME')
        Hsampled.SetLineColor(2)

        Hsmeared.Draw('SAME')
        Hsmeared.SetLineColor(3)

        self.Save( self.c, 'truthSampledSmeared' )


        # -------------------
        # Legend

        Htruth.SetName('truth')
        Hsampled.SetName('sampled')
        Hsmeared.SetName('smeared')

        self.leg = ROOT.TLegend( 0.1, 0.5, 0.42, 0.9 )
        self.leg.AddEntry( 'truth', 'truth' )
        self.leg.AddEntry( 'sampled', 'sampled' )
        self.leg.AddEntry( 'smeared', 'smeared' )
        self.leg.Draw('SAME')


        # Add used histograms to python list for persistence
        self.Persistence.extend([ Htruth, Hsampled, Hsmeared ])


    def BayesUnfoldTest( self ):

        # Make the basic plots
        self.PlotBasicHists()

        # Make the histogram with the pure smeared counts
        Hmeas = self.HistogramFromAxes( self.xAxis, self.smearedCountAxis, normalize=False )

        # Plot for several k's
        colors = [ 40, 41, 42, 43, 46 ]
        for k in [ 1000, 100, 10 ]:

            # Response object
            R = ROOT.RooUnfoldResponse( self.nBins, self.xAxis[0], self.xAxis[-1] )

            for ( xTrue, xMeas ) in self.events:
                R.Fill( xMeas, xTrue )


            unfold = ROOT.RooUnfoldBayes( R, Hmeas, k )

            Hunfolded = unfold.Hreco()

            Hunfolded.Scale( 1./Hunfolded.Integral('width') )

            Hunfolded.Draw('SAME')
            Hunfolded.SetLineColor( colors.pop(0) )
            Hunfolded.SetLineWidth(2)

            # Add entry to legend
            Hname = GetUniqueRootName()
            Hunfolded.SetName( Hname )
            self.leg.AddEntry( Hname, 'Bayes (N={0})'.format(k) )

            self.Persistence.append( deepcopy(Hunfolded) )


        self.Save( self.c, 'Bayes' )


    def SVDUnfoldTest( self ):

        # Make the basic plots
        self.PlotBasicHists()

        # Make the histogram with the pure smeared counts
        Hmeas = self.HistogramFromAxes( self.xAxis, self.smearedCountAxis, normalize=False )


        # Plot for several k's
        colors = [ 40, 41, 42, 43, 46 ]
        for k in [ 1, 5, 10 ]:

            # Response object
            R = ROOT.RooUnfoldResponse( self.nBins, self.xAxis[0], self.xAxis[-1] )

            for ( xTrue, xMeas ) in self.events:
                R.Fill( xMeas, xTrue )


            unfold = ROOT.RooUnfoldSvd( R, Hmeas, k )

            Hunfolded = unfold.Hreco()

            Hunfolded.Scale( 1./Hunfolded.Integral('width') )

            Hunfolded.Draw('SAME')
            Hunfolded.SetLineColor( colors.pop(0) )
            Hunfolded.SetLineWidth(2)

            # Add entry to legend
            Hname = GetUniqueRootName()
            Hunfolded.SetName( Hname )
            self.leg.AddEntry( Hname, 'SVD (k={0})'.format(k) )

            self.Persistence.append( deepcopy(Hunfolded) )


        self.Save( self.c, 'SVD' )


    def BinByBinUnfoldTest( self ):

        # Make the basic plots
        self.PlotBasicHists()

        # Make the histogram with the pure smeared counts
        Hmeas = self.HistogramFromAxes( self.xAxis, self.smearedCountAxis, normalize=False )


        # Response object
        R = ROOT.RooUnfoldResponse( self.nBins, self.xAxis[0], self.xAxis[-1] )

        for ( xTrue, xMeas ) in self.events:
            R.Fill( xMeas, xTrue )


        unfold = ROOT.RooUnfoldBinByBin( R, Hmeas )

        Hunfolded = unfold.Hreco()

        Hunfolded.Scale( 1./Hunfolded.Integral('width') )

        Hunfolded.Draw('SAME')
        Hunfolded.SetLineColor( 9 )
        Hunfolded.SetLineWidth(2)

        # Add entry to legend
        Hname = GetUniqueRootName()
        Hunfolded.SetName( Hname )
        self.leg.AddEntry( Hname, 'BinByBin' )

        self.Persistence.append( deepcopy(Hunfolded) )

        self.Save( self.c, 'BinByBin' )





    def InvertedKUnfoldTest( self ):

        # Make the basic plots
        self.PlotBasicHists()

        # ======================================
        # Unfold by inverting K

        self.EstimateK()
        Kinv = numpy.linalg.inv( self.K )

        unfoldedKinvAxis = numpy.dot( Kinv, self.smearedTruthAxis )

        HunfoldedKinv = self.HistogramFromAxes( self.xAxis, unfoldedKinvAxis, normalize=False )

        HunfoldedKinv.Draw('SAME')
        HunfoldedKinv.SetLineColor(6)
        HunfoldedKinv.SetLineStyle(2)


        globalMax = 1.1 * max([
            max(self.truthAxis),
            max(self.sampledTruthAxis),
            max(self.smearedTruthAxis),
            max(unfoldedKinvAxis),
            ])

        globalMin = 1.1 * min([
            min(self.truthAxis),
            min(self.sampledTruthAxis),
            min(self.smearedTruthAxis),
            min(unfoldedKinvAxis),
            ])

        self.base.SetMaximum( globalMax )
        self.base.SetMinimum( globalMin )


        self.Save( self.c, 'Kinversion' )


    def EstimateK( self ):

        nIter = 10000

        # Initialize 2D empty matrix
        self.K = [ [ 0. for j in xrange(self.nBins) ] for i in xrange(self.nBins) ]


        for iBin in xrange(self.nBins):

            # central x value for this bin
            xCentral = self.xAxisCenters[iBin]

            # Smear it nIter times, count to which bins it goes
            for iIter in xrange(nIter):
                jBinSmeared = self.GetClosestIndex( self.Smear(xCentral), self.xAxisCenters )
                self.K[iBin][jBinSmeared] += 1

            self.K[iBin] = [ float(count)/nIter for count in self.K[iBin] ]


    def HistogramFromAxes( self, xAxis, yAxis, normalize=False ):
        
        if len(xAxis) == len(yAxis)+1:
            'This is fine'
        elif len(xAxis) == len(yAxis):
            print '  Warning in HistogramFromAxes: Number of points in xAxis and yAxis equal'
            print '  xAxis should be bin boundaries (length nBins+1) and yAxis bin contents (length nBins)'
            print '  Now ignoring the last entry in yAxis'
            yAxis = yAxis[:-1]
        else:
            print '  Error in HistogramFromAxes: Passed axes have wrong length'
            print '  len(xAxis) = {0}  |  len(yAxis) = {1}'.format( len(xAxis), len(yAxis) )

        nBins = len(xAxis)-1

        Hname = GetUniqueRootName()
        H = ROOT.TH1F( Hname, Hname, nBins, array( 'f', xAxis )   )

        if normalize: sumYAxis = float(sum(yAxis))

        for iBin in xrange(nBins):
            if not normalize:
                H.SetBinContent( iBin+1, yAxis[iBin] )
            else:
                binWidth = float(xAxis[iBin+1]-xAxis[iBin])
                H.SetBinContent( iBin+1, yAxis[iBin] / (sumYAxis*binWidth) )

        H.SetLineWidth(2)
        ROOT.SetOwnership( H, False )
        return H


    def Save( self, c, outname, png=False ):

        if not '.pdf' in outname:
            outname += '.pdf'

        if self.PlotPrefix:
            outname = self.PlotPrefixText + outname

        if not os.path.isdir(self.plotdir): os.makedirs(self.plotdir)

        c.SaveAs( os.path.join( self.plotdir, outname ) )


    def MakeAxis( self, xmin, xmax, nBins ):
        # Convention: Axis is a list of bin BOUNDARIES
        #   nBoundaries = nBins + 1 !
        binWidth = float(xmax-xmin) / (nBins)

        # These are the bounds
        axis = [ xmin + i * binWidth for i in xrange(nBins+1) ]

        # These are the centers
        centers = [ xmin + 0.5*binWidth + i * binWidth for i in xrange(nBins) ]

        return axis, centers



    # Finds the closest x-value and index in a list, given a certain x
    def GetClosestIndex( self, xFit, xAxis ):

        closest_x = 9999999.
        closest_i = None

        for i, x in enumerate( xAxis ):

            if abs(x-xFit) < closest_x:
                closest_x = abs(x-xFit)
                closest_i = i

        return closest_i

