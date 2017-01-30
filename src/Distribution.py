#!/usr/bin/env python
"""
Thomas:
"""

########################################
# Imports
########################################

import os
import ROOT
import random


class Distribution:

    # ======================================
    # Initialization

    def __init__( self, xAxisCenters, yAxis ):
        'Initialization'

        self.xAxisCenters = xAxisCenters
        self.yAxis = yAxis
        self.nBins = len(xAxisCenters)

        if not len(xAxisCenters) == len(yAxis):
            print 'Lengths of axes not the same; this very probably a mistake'


    def DrawSample( self ):

        # # First make the CDF
        # self.CDFAxis = [ 0. for i in xrange(self.nBins) ]

        # Sum = 0.
        # for iBin in xrange(self.nBins):
        #     Sum += self.yAxis[iBin]
        #     self.CDFAxis[iBin] = Sum


        # Use maximum of yAxis as normalization
        yMax = float(max(self.yAxis))


        Accepted = False

        while not Accepted:

            # Select a random bin
            iBin = random.randint( 0., self.nBins-1 )

            acceptanceProbability = self.yAxis[iBin] / yMax

            # Draw random number between 0 and 1
            if random.random() <= acceptanceProbability: Accepted = True


        xDrawn = self.xAxisCenters[iBin]
        yDrawn = self.yAxis[iBin]

        return ( iBin, xDrawn, yDrawn )




