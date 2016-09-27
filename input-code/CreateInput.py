#python2.7

#======================== Welcome to CreateInput ==========================#
# CreateInput is a python based program for creating
# the input files for classifier from the 
# output of LZIFU (Ho et al. 2016)
# This is not a pretty program
#==========================================================================#

#================================ Author ==================================#
# Elise Hampton
# PhD Candiadate Australian National University
# Research School of Astronomy and Astrophysics
# CreateInput written for the S7 collaboration
# (PI: A/Prof Michael Dopita)
# Last update: July 2016
#==========================================================================#

#================================ Updates =================================#
# 17th Dec 2014 - created
# 17th Dec 2014 - runs through several galaxies
# 12th Jan 2015 - Comment definitions and update to not-train
# 13th Jan 2015 - galaxy IDs now read in through galaxies.txt
# 13th Jan 2015 - now loops over dof21 not qualmask so can run
# without testing
# 19th Jan 2015 - Takes in traina nd test as command line arguments
# 5th May 2015  - Rewritten for SAMI data
#==========================================================================#

import pyfits
import numpy as np
import math
import argparse
#import statistics


#============================= Definitions ================================#

def writetostring(A):
    """ writetostring takes A and writes it out to a string
    Input: A
    Output: strings
    """
    
    strings = ''
    for n in range(0,len(A)):

        if math.isnan(A[n]):
            strings = strings+' '+'0.0'
        else:
            strings = strings+' '+str(A[n])

    return strings

def calcradvalue(central_pix, other_pix, i, Re):
    """ calcradvalue calculates the radius value of a particular pixel
    Input: central pix, pixel in questions, inclination of galaxy, and effective radius
    Output: rad
    UPDATE: Now unused Dec 2014
    """
    pixtoarcsec = 0.5
    Re_pix = Re/pixtoarcsec
    #difference between two points
    r = math.sqrt(math.pow(central_pix[0] - other_pix[0],2)+math.pow(central_pix[1] - other_pix[1],2))
    rr = math.fabs(r*math.cos(r))
    rad = rr/Re_pix
    
    return rad

def writeRatstoString(hbeta1,hbeta1_err,hbeta20,hbeta20_err,hbeta21,hbeta21_err,hbeta22,hbeta22_err,hbeta30,hbeta30_err,hbeta31,hbeta31_err,hbeta32,hbeta32_err,
hbeta33,hbeta33_err):
    """ writeRatstoString calculates ratios based on the input and puts them together in a string
    Input: fluxes and flux errors of a particular emission line
    Output: out = SN + totrat
    """
    SN = ''
    totrat = ''
    if math.isnan(hbeta1[n,nn]/hbeta1_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta1[n,nn]/hbeta1_err[n,nn])
                
    if math.isnan(hbeta20[n,nn]/hbeta20_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta20[n,nn]/hbeta20_err[n,nn])
                
    if math.isnan(hbeta21[n,nn]/hbeta21_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta21[n,nn]/hbeta21_err[n,nn])

    if math.isnan(hbeta22[n,nn]/hbeta22_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta22[n,nn]/hbeta22_err[n,nn])

    if math.isnan(hbeta30[n,nn]/hbeta30_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta30[n,nn]/hbeta30_err[n,nn])

    if math.isnan(hbeta31[n,nn]/hbeta31_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta31[n,nn]/hbeta31_err[n,nn])

    if math.isnan(hbeta32[n,nn]/hbeta32_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta32[n,nn]/hbeta32_err[n,nn])

    if math.isnan(hbeta33[n,nn]/hbeta33_err[n,nn]):
        SN = SN + ' 0.0'
    else:
        SN = SN + ' ' + str(hbeta33[n,nn]/hbeta33_err[n,nn])
            

    if math.isnan(hbeta1[n,nn]/hbeta20[n,nn]) or hbeta20[n,nn] == 0.0:
        totrat = totrat + ' 0.0'
    else:
        totrat = totrat + ' ' + str(hbeta1[n,nn]/hbeta20[n,nn])

    if math.isnan(hbeta1[n,nn]/hbeta30[n,nn]) or hbeta30[n,nn] == 0.0:
        totrat = totrat + ' 0.0'
    else:
        totrat = totrat + ' ' + str(hbeta1[n,nn]/hbeta30[n,nn])

    out = SN + ' ' + totrat
        
    return out

#=============== CreateInput ===================#
#IDs = ['IC1816','NGC7590','MARK573','MCG-05-23-004','NGC613','NGC5728','NGC6000']
# make this an input to the program so no adjustments have to be made in the code
# Read in a text file with these IDs listed
galaxylistF = open('galaxies.txt','r')
fie = galaxylistF.read()
IDs = fie.split('\n')
IDs = IDs[0:len(IDs)-1]
print IDs
#stop

#read in command line arguments train and test
parser = argparse.ArgumentParser(description='CreateInput')
parser.add_argument('train',help='train yes or no')
parser.add_argument('test',help='test yes or no')
args = parser.parse_args()
test = int(args.test)
train = int(args.train)


# When in training mode extra files are created for testing the training data
# this includes the use of y files that contain known answers
#train = True
if train == 1:
    outfile2 = open('output.txt','w')
    outfile3 = open('output_pix.txt','w')
    cvfile = open('input_cv.txt','w')
    testfile = open('input_test.txt','w')
    cvfileout = open('output_cv.txt','w')
    testfileout = open('output_test.txt','w')
    cvfilepix = open('pix_cv.txt','w')
    testfilepix = open('pix_test.txt','w')

#=====================================================#

outfile = open('input.txt','w')

zerocomps = []
onecomps = []
twocomps = []
threecomps = []
badcomps = []
zerocompsL = []
onecompsL = []
twocompsL = []
threecompsL = []
badcompsL = []

pix1 = []

pix2 = []

pix3 = []

pix0 = []

pixb = []


for row in IDs:
    ID = row
    #================== lzifu info ======================#
    #open lzifu output files
    #comp1fits = pyfits.open('../data/'+str(ID)+'_1_comp.fits')
    #print str(ID)
    #comp2fits = pyfits.open('../data/'+str(ID)+'_2_comp.fits')
    #comp3fits = pyfits.open('../data/'+str(ID)+'_3_comp.fits')
    comp1fits = pyfits.open('../../SAMI/sami_v0.9/'+str(ID)+'_1_comp.fits')
    print str(ID)
    comp2fits = pyfits.open('../../SAMI/sami_v0.9/'+str(ID)+'_2_comp.fits')
    comp3fits = pyfits.open('../../SAMI/sami_v0.9/'+str(ID)+'_3_comp.fits')
    header2 = comp3fits[0].header
    #print header2
    #stop

    #err1ha = pyfits.open('../datacopy/newerrors/'+str(ID)+'_1_comp_haerr.fits')
    #err2ha = pyfits.open('../datacopy/newerrors/'+str(ID)+'_2_comp_haerr.fits')
    #err3ha = pyfits.open('../datacopy/newerrors/'+str(ID)+'_3_comp_haerr.fits')
    #err1hb = pyfits.open('../datacopy/newerrors/'+str(ID)+'_1_comp_hberr.fits')
    #err2hb = pyfits.open('../datacopy/newerrors/'+str(ID)+'_2_comp_hberr.fits')
    #err3hb = pyfits.open('../datacopy/newerrors/'+str(ID)+'_3_comp_hberr.fits')
    
    #print header2
    #stop
    #pull out values
    halpha1 = comp1fits['HALPHA'].data[0,:,:]
    halpha1_err = comp1fits['HALPHA_ERR'].data[0,:,:]
    halpha20 = comp2fits['HALPHA'].data[0,:,:]
    halpha20_err = comp2fits['HALPHA_ERR'].data[0,:,:]
    halpha21 = comp2fits['HALPHA'].data[1,:,:]
    halpha21_err = comp2fits['HALPHA_ERR'].data[1,:,:]
    halpha22 = comp2fits['HALPHA'].data[2,:,:]
    halpha22_err = comp2fits['HALPHA_ERR'].data[2,:,:]
    halpha30 = comp3fits['HALPHA'].data[0,:,:]
    halpha30_err = comp3fits['HALPHA_ERR'].data[0,:,:]
    halpha31 = comp3fits['HALPHA'].data[1,:,:]
    halpha31_err = comp3fits['HALPHA_ERR'].data[1,:,:]
    halpha32 = comp3fits['HALPHA'].data[2,:,:]
    halpha32_err = comp3fits['HALPHA_ERR'].data[2,:,:]
    halpha33 = comp3fits['HALPHA'].data[3,:,:]
    halpha33_err = comp3fits['HALPHA_ERR'].data[3,:,:]

    hbeta1 = comp1fits['HBETA'].data[0,:,:]
    hbeta1_err = comp1fits['HBETA_ERR'].data[0,:,:]
    hbeta20 = comp2fits['HBETA'].data[0,:,:]
    hbeta20_err = comp2fits['HBETA_ERR'].data[0,:,:]
    hbeta21 = comp2fits['HBETA'].data[1,:,:]
    hbeta21_err = comp2fits['HBETA_ERR'].data[1,:,:]
    hbeta22 = comp2fits['HBETA'].data[2,:,:]
    hbeta22_err = comp2fits['HBETA_ERR'].data[2,:,:]
    hbeta30 = comp3fits['HBETA'].data[0,:,:]
    hbeta30_err = comp3fits['HBETA_ERR'].data[0,:,:]
    hbeta31 = comp3fits['HBETA'].data[1,:,:]
    hbeta31_err = comp3fits['HBETA_ERR'].data[1,:,:]
    hbeta32 = comp3fits['HBETA'].data[2,:,:]
    hbeta32_err = comp3fits['HBETA_ERR'].data[2,:,:]
    hbeta33 = comp3fits['HBETA'].data[3,:,:]
    hbeta33_err = comp3fits['HBETA_ERR'].data[3,:,:]

    nii65831 = comp1fits['NII6583'].data[0,:,:]
    nii65831_err = comp1fits['NII6583_ERR'].data[0,:,:]
    nii658320 = comp2fits['NII6583'].data[0,:,:]
    nii658320_err = comp2fits['NII6583_ERR'].data[0,:,:]
    nii658321 = comp2fits['NII6583'].data[1,:,:]
    nii658321_err = comp2fits['NII6583_ERR'].data[1,:,:]
    nii658322 = comp2fits['NII6583'].data[2,:,:]
    nii658322_err = comp2fits['NII6583_ERR'].data[2,:,:]
    nii658330 = comp3fits['NII6583'].data[0,:,:]
    nii658330_err = comp3fits['NII6583_ERR'].data[0,:,:]
    nii658331 = comp3fits['NII6583'].data[1,:,:]
    nii658331_err = comp3fits['NII6583_ERR'].data[1,:,:]
    nii658332 = comp3fits['NII6583'].data[2,:,:]
    nii658332_err = comp3fits['NII6583_ERR'].data[2,:,:]
    nii658333 = comp3fits['NII6583'].data[3,:,:]
    nii658333_err = comp3fits['NII6583_ERR'].data[3,:,:]

    sii67161 = comp1fits['SII6716'].data[0,:,:]
    sii67161_err = comp1fits['SII6716_ERR'].data[0,:,:]
    sii671620 = comp2fits['SII6716'].data[0,:,:]
    sii671620_err = comp2fits['SII6716_ERR'].data[0,:,:]
    sii671621 = comp2fits['SII6716'].data[1,:,:]
    sii671621_err = comp2fits['SII6716_ERR'].data[1,:,:]
    sii671622 = comp2fits['SII6716'].data[2,:,:]
    sii671622_err = comp2fits['SII6716_ERR'].data[2,:,:]
    sii671630 = comp3fits['SII6716'].data[0,:,:]
    sii671630_err = comp3fits['SII6716_ERR'].data[0,:,:]
    sii671631 = comp3fits['SII6716'].data[1,:,:]
    sii671631_err = comp3fits['SII6716_ERR'].data[1,:,:]
    sii671632 = comp3fits['SII6716'].data[2,:,:]
    sii671632_err = comp3fits['SII6716_ERR'].data[2,:,:]
    sii671633 = comp3fits['SII6716'].data[3,:,:]
    sii671633_err = comp3fits['SII6716_ERR'].data[3,:,:]

    sii67311 = comp1fits['SII6731'].data[0,:,:]
    sii67311_err = comp1fits['SII6731_ERR'].data[0,:,:]
    sii673120 = comp2fits['SII6731'].data[0,:,:]
    sii673120_err = comp2fits['SII6731_ERR'].data[0,:,:]
    sii673121 = comp2fits['SII6731'].data[1,:,:]
    sii673121_err = comp2fits['SII6731_ERR'].data[1,:,:]
    sii673122 = comp2fits['SII6731'].data[2,:,:]
    sii673122_err = comp2fits['SII6731_ERR'].data[2,:,:]
    sii673130 = comp3fits['SII6731'].data[0,:,:]
    sii673130_err = comp3fits['SII6731_ERR'].data[0,:,:]
    sii673131 = comp3fits['SII6731'].data[1,:,:]
    sii673131_err = comp3fits['SII6731_ERR'].data[1,:,:]
    sii673132 = comp3fits['SII6731'].data[2,:,:]
    sii673132_err = comp3fits['SII6731_ERR'].data[2,:,:]
    sii673133 = comp3fits['SII6731'].data[3,:,:]
    sii673133_err = comp3fits['SII6731_ERR'].data[3,:,:]

    #nii65481 = comp1fits['NII6548'].data[0,:,:]
    #nii65481_err = comp1fits['NII6548_ERR'].data[0,:,:]
    #nii654820 = comp2fits['NII6548'].data[0,:,:]
    #nii654820_err = comp2fits['NII6548_ERR'].data[0,:,:]
    #nii654821 = comp2fits['NII6548'].data[1,:,:]
    #nii654821_err = comp2fits['NII6548_ERR'].data[1,:,:]
    #nii654822 = comp2fits['NII6548'].data[2,:,:]
    #nii654822_err = comp2fits['NII6548_ERR'].data[2,:,:]
    #nii654830 = comp3fits['NII6548'].data[0,:,:]
    #nii654830_err = comp3fits['NII6548_ERR'].data[0,:,:]
    #nii654831 = comp3fits['NII6548'].data[1,:,:]
    #nii654831_err = comp3fits['NII6548_ERR'].data[1,:,:]
    #nii654832 = comp3fits['NII6548'].data[2,:,:]
    #nii654832_err = comp3fits['NII6548_ERR'].data[2,:,:]
    #nii654833 = comp3fits['NII6548'].data[3,:,:]
    #nii654833_err = comp3fits['NII6548_ERR'].data[3,:,:]

    oiii50071 = comp1fits['OIII5007'].data[0,:,:]
    oiii50071_err = comp1fits['OIII5007_ERR'].data[0,:,:]
    oiii500720 = comp2fits['OIII5007'].data[0,:,:]
    oiii500720_err = comp2fits['OIII5007_ERR'].data[0,:,:]
    oiii500721 = comp2fits['OIII5007'].data[1,:,:]
    oiii500721_err = comp2fits['OIII5007_ERR'].data[1,:,:]
    oiii500722 = comp2fits['OIII5007'].data[2,:,:]
    oiii500722_err = comp2fits['OIII5007_ERR'].data[2,:,:]
    oiii500730 = comp3fits['OIII5007'].data[0,:,:]
    oiii500730_err = comp3fits['OIII5007_ERR'].data[0,:,:]
    oiii500731 = comp3fits['OIII5007'].data[1,:,:]
    oiii500731_err = comp3fits['OIII5007_ERR'].data[1,:,:]
    oiii500732 = comp3fits['OIII5007'].data[2,:,:]
    oiii500732_err = comp3fits['OIII5007_ERR'].data[2,:,:]
    oiii500733 = comp3fits['OIII5007'].data[3,:,:]
    oiii500733_err = comp3fits['OIII5007_ERR'].data[3,:,:]

    #oiii49591 = comp1fits['OIII4959'].data[0,:,:]
    #oiii49591_err = comp1fits['OIII4959_ERR'].data[0,:,:]
    #oiii495920 = comp2fits['OIII4959'].data[0,:,:]
    #oiii495920_err = comp2fits['OIII4959_ERR'].data[0,:,:]
    #oiii495921 = comp2fits['OIII4959'].data[1,:,:]
    #oiii495921_err = comp2fits['OIII4959_ERR'].data[1,:,:]
    #oiii495922 = comp2fits['OIII4959'].data[2,:,:]
    #oiii495922_err = comp2fits['OIII4959_ERR'].data[2,:,:]
    #oiii495930 = comp3fits['OIII4959'].data[0,:,:]
    #oiii495930_err = comp3fits['OIII4959_ERR'].data[0,:,:]
    #oiii495931 = comp3fits['OIII4959'].data[1,:,:]
    #oiii495931_err = comp3fits['OIII4959_ERR'].data[1,:,:]
    #oiii495932 = comp3fits['OIII4959'].data[2,:,:]
    #oiii495932_err = comp3fits['OIII4959_ERR'].data[2,:,:]
    #oiii495933 = comp3fits['OIII4959'].data[3,:,:]
    #oiii495933_err = comp3fits['OIII4959_ERR'].data[3,:,:]

    v1 = comp1fits['V'].data[0,:,:]
    v1_err = comp1fits['V_ERR'].data[0,:,:]
    v20 = comp2fits['V'].data[0,:,:]
    v20_err = comp2fits['V_ERR'].data[0,:,:]
    v21 = comp2fits['V'].data[1,:,:]
    v21_err = comp2fits['V_ERR'].data[1,:,:]
    v22 = comp2fits['V'].data[2,:,:]
    v22_err = comp2fits['V_ERR'].data[2,:,:]
    v30 = comp3fits['V'].data[0,:,:]
    v30_err = comp3fits['V_ERR'].data[0,:,:]
    v31 = comp3fits['V'].data[1,:,:]
    v31_err = comp3fits['V_ERR'].data[1,:,:]
    v32 = comp3fits['V'].data[2,:,:]
    v32_err = comp3fits['V_ERR'].data[2,:,:]
    v33 = comp3fits['V'].data[3,:,:]
    v33_err = comp3fits['V_ERR'].data[3,:,:]

    vdisp1 = comp1fits['VDISP'].data[0,:,:]
    vdisp1_err = comp1fits['VDISP_ERR'].data[0,:,:]
    vdisp20 = comp2fits['VDISP'].data[0,:,:]
    vdisp20_err = comp2fits['VDISP_ERR'].data[0,:,:]
    vdisp21 = comp2fits['VDISP'].data[1,:,:]
    vdisp21_err = comp2fits['VDISP_ERR'].data[1,:,:]
    vdisp22 = comp2fits['VDISP'].data[2,:,:]
    vdisp22_err = comp2fits['VDISP_ERR'].data[2,:,:]
    vdisp30 = comp3fits['VDISP'].data[0,:,:]
    vdisp30_err = comp3fits['VDISP_ERR'].data[0,:,:]
    vdisp31 = comp3fits['VDISP'].data[1,:,:]
    vdisp31_err = comp3fits['VDISP_ERR'].data[1,:,:]
    vdisp32 = comp3fits['VDISP'].data[2,:,:]
    vdisp32_err = comp3fits['VDISP_ERR'].data[2,:,:]
    vdisp33 = comp3fits['VDISP'].data[3,:,:]
    vdisp33_err = comp3fits['VDISP_ERR'].data[3,:,:]

    chi21 = comp1fits['CHI2'].data
    chi22 = comp2fits['CHI2'].data
    chi23 = comp3fits['CHI2'].data

    contchi21 = comp1fits['CONT_CHI2'].data
    contchi22 = comp2fits['CONT_CHI2'].data
    contchi23 = comp3fits['CONT_CHI2'].data

    dof21 = comp1fits['DOF'].data
    dof22 = comp2fits['DOF'].data
    dof23 = comp3fits['DOF'].data

    #starebv1 = comp1fits['STAR_EBV'].data
    #starebv2 = comp2fits['STAR_EBV'].data
    #starebv3 = comp3fits['STAR_EBV'].data

    hbebv1 = comp1fits['HB_CONT_EW'].data
    hbebv2 = comp2fits['HB_CONT_EW'].data
    hbebv3 = comp3fits['HB_CONT_EW'].data

    haebv1 = comp1fits['HA_CONT_EW'].data
    haebv2 = comp2fits['HA_CONT_EW'].data
    haebv3 = comp3fits['HA_CONT_EW'].data

    #=====================================================#
    comp1fits.close()
    comp2fits.close()
    comp3fits.close()

    #============= mask files for training ===============#
    #In training mode open mask file and outputfile for NN
    #train = True
    if train == 1:
        qualmask1 = pyfits.open('../Compmaps/'+str(ID)+'_compmask-comb.fits')[0].data
        
    #=====================================================#
    #sami = False
    #if train == 1:
    nanmask = pyfits.open('masks/'+str(ID)+'_nanmask.fits')[0].data

    siz = np.shape(dof21)

    m = siz[0]*siz[1]
    #perc50 = int((m/100.0)*50)
    perc50 = int((m/100.0)*75)
    perc25 = int((m/100.0)*25)
    #print m, perc50, perc25
    #print siz[1]
    #stop
    split = 0
    for n in range(0,siz[0]):
        for nn in range(0,siz[1]):
            split = split + 1
            if math.isnan(chi21[n,nn]):
                chis = '0.0 0.0'
            else:
                chis = str(chi21[n,nn]) + ' ' + str(chi21[n,nn]/dof21[n,nn])

            if math.isnan(chi22[n,nn]):
                chis = chis + ' 0.0 0.0'
            else:
                chis = chis + ' ' + str(chi22[n,nn]) + ' ' + str(chi22[n,nn]/dof22[n,nn])

            if math.isnan(chi23[n,nn]):
                chis = chis + ' 0.0 0.0'
            else:
                chis = chis + ' ' + str(chi23[n,nn]) + ' ' + str(chi23[n,nn]/dof23[n,nn])

            if math.isnan(contchi21[n,nn]):
                chis = chis + ' 0.0'
            else:
                chis = chis + ' ' +str(contchi21[n,nn]) 

            if math.isnan(contchi22[n,nn]):
                chis = chis + ' 0.0'
            else:
                chis = chis + ' ' + str(contchi22[n,nn])

            if math.isnan(contchi23[n,nn]):
                chis = chis + ' 0.0'
            else:
                chis = chis + ' ' + str(contchi23[n,nn]) 

            if math.isnan(hbebv1[n,nn]):
                stars = '0.0'
            else:
                stars = str(hbebv1[n,nn])
            if math.isnan(hbebv2[n,nn]):
                stars = stars + ' 0.0'
            else:
                stars = stars + ' ' + str(hbebv2[n,nn])
            if math.isnan(hbebv3[n,nn]):
                stars = stars + ' 0.0'
            else:
                stars = stars + ' ' + str(hbebv3[n,nn])

            if math.isnan(haebv1[n,nn]):
                stars = '0.0'
            else:
                stars = str(haebv1[n,nn])
            if math.isnan(haebv2[n,nn]):
                stars = stars + ' 0.0'
            else:
                stars = stars + ' ' + str(haebv2[n,nn])
            if math.isnan(haebv3[n,nn]):
                stars = stars + ' 0.0'
            else:
                stars = stars + ' ' + str(haebv3[n,nn])

            SN = ''
            totrat = ''
            SN = SN + writeRatstoString(halpha1,halpha1_err,halpha20,halpha20_err,halpha21,halpha21_err,halpha22,halpha22_err,halpha30,halpha30_err,halpha31,halpha31_err,halpha32,halpha32_err,halpha33,halpha33_err)
            SN = SN + writeRatstoString(nii65831,nii65831_err,nii658320,nii658320_err,nii658321,nii658321_err,nii658322,nii658322_err,nii658330,nii658330_err,nii658331,nii658331_err,nii658332,nii658332_err,nii658333,nii658333_err)
            SN = SN + writeRatstoString(sii67161,sii67161_err,sii671620,sii671620_err,sii671621,sii671621_err,sii671622,sii671622_err,sii671630,sii671630_err,sii671631,sii671631_err,sii671632,sii671632_err,sii671633,sii671633_err)
            SN = SN + writeRatstoString(sii67311,sii67311_err,sii673120,sii673120_err,sii673121,sii673121_err,sii673122,sii673122_err,sii673130,sii673130_err,sii673131,sii673131_err,sii673132,sii673132_err,sii673133,sii673133_err)
            #SN = SN + writeRatstoString(nii65481,nii65481_err,nii654820,nii654820_err,nii654821,nii654821_err,nii654822,nii654822_err,nii654830,nii654830_err,nii654831,nii654831_err,nii654832,nii654832_err,nii654833,nii654833_err)
            SN = SN + writeRatstoString(hbeta1,hbeta1_err,hbeta20,hbeta20_err,hbeta21,hbeta21_err,hbeta22,hbeta22_err,hbeta30,hbeta30_err,hbeta31,hbeta31_err,hbeta32,hbeta32_err,hbeta33,hbeta33_err)
            #SN = SN + writeRatstoString(oii37261,oii37261_err,oii372620,oii372620_err,oii372621,oii372621_err,oii372622,oii372622_err,oii372630,oii372630_err,oii372631,oii372631_err,oii372632,oii372632_err,oii372633,oii372633_err)
            #SN = SN + writeRatstoString(oii37291,oii37291_err,oii372920,oii372920_err,oii372921,oii372921_err,oii372922,oii372922_err,oii372930,oii372930_err,oii372931,oii372931_err,oii372932,oii372932_err,oii372933,oii372933_err)
            #SN = SN + writeRatstoString(neiii38691,neiii38691_err,neiii386920,neiii386920_err,neiii386921,neiii386921_err,neiii386922,neiii386922_err,neiii386930,neiii386930_err,neiii386931,neiii386931_err,neiii386932,neiii386932_err,neiii386933,neiii386933_err)
            SN = SN + writeRatstoString(oiii50071,oiii50071_err,oiii500720,oiii500720_err,oiii500721,oiii500721_err,oiii500722,oiii500722_err,oiii500730,oiii500730_err,oiii500731,oiii500731_err,oiii500732,oiii500732_err,oiii500733,oiii500733_err)

            if math.isnan(v1[n,nn]/v1_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v1[n,nn]/v1_err[n,nn])
                
            if math.isnan(v20[n,nn]/v20_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v20[n,nn]/v20_err[n,nn])
                
            if math.isnan(v21[n,nn]/v21_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v21[n,nn]/v21_err[n,nn])

            if math.isnan(v22[n,nn]/v22_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v22[n,nn]/v22_err[n,nn])

            if math.isnan(v30[n,nn]/v30_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v30[n,nn]/v30_err[n,nn])

            if math.isnan(v31[n,nn]/v31_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v31[n,nn]/v31_err[n,nn])

            if math.isnan(v32[n,nn]/v32_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v32[n,nn]/v32_err[n,nn])

            if math.isnan(v33[n,nn]/v33_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(v33[n,nn]/v33_err[n,nn])

            if math.isnan(vdisp1[n,nn]/vdisp1_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' +nstr(vdisp1[n,nn]/vdisp1_err[n,nn])
                
            if math.isnan(vdisp20[n,nn]/vdisp20_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(vdisp20[n,nn]/vdisp20_err[n,nn])
                
            if math.isnan(vdisp21[n,nn]/vdisp21_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(vdisp21[n,nn]/vdisp21_err[n,nn])

            if math.isnan(vdisp22[n,nn]/vdisp22_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(vdisp22[n,nn]/vdisp22_err[n,nn])

            if math.isnan(vdisp30[n,nn]/vdisp30_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(vdisp30[n,nn]/vdisp30_err[n,nn])

            if math.isnan(vdisp31[n,nn]/vdisp31_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(vdisp31[n,nn]/vdisp31_err[n,nn])

            if math.isnan(vdisp32[n,nn]/vdisp32_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(vdisp32[n,nn]/vdisp32_err[n,nn])

            if math.isnan(vdisp33[n,nn]/vdisp33_err[n,nn]):
                SN = SN + ' 0.0'
            else:
                SN = SN + ' ' + str(vdisp33[n,nn]/vdisp33_err[n,nn])
                
            #pos = calcradvalue(centre,[n,nn],i,Re)
            # nolonger used Dec 2014
            
            toPrint = chis + ' ' + stars + ' ' + SN + ' ' + totrat + '\n'


            if train == 1:
                if nanmask[n,nn] == -5:
                    tempor = 1
                else:
                    #outfile.write(toPrint)
                    if int(qualmask1[n,nn]) == 0:
                        #outfile.write(toPrint)
                        #zerocomps.append(toPrint)
                        #outfile2.write(str(4)+'\n')
                        #zerocompsL.append(str(4)+'\n')
                        #pix0.append(str(n)+' '+str(nn)+' '+str(ID)+'\n')
                        onecomps.append(toPrint)
                        onecompsL.append(str(1)+'\n')
                        pix1.append(str(n)+' '+str(nn)+' '+str(ID)+'\n')
                        #zerocomps.append(toPrint)
                        #zerocompsL.append(str(4)+'\n')
                        #pix0.append(str(n)+' '+str(nn)+' '+str(ID)+'\n')
                    elif int(qualmask1[n,nn]) == -1:
                        #outfile.write(toPrint)
                        badcomps.append(toPrint)
                        #outfile2.write(str(5)+'\n')
                        badcompsL.append(str(5)+'\n')
                        pixb.append(str(n)+' '+str(nn)+' '+str(ID)+'\n')
                    elif int(qualmask1[n,nn]) == -2 or int(qualmask1[n,nn]) == -3 or int(qualmask1[n,nn]) == -4 or int(qualmask1[n,nn]) == -5:
                        print 'dont save anything'
                    else:
                        #outfile.write(toPrint)
                        #outfile2.write(str(int(qualmask1[n,nn]))+'\n')
                        if int(qualmask1[n,nn]) == 1:
                            onecomps.append(toPrint)
                            onecompsL.append(str(int(qualmask1[n,nn]))+'\n')
                            pix1.append(str(n)+' '+str(nn)+' '+str(ID)+'\n')
                        elif int(qualmask1[n,nn]) == 2:
                            twocomps.append(toPrint)
                            twocompsL.append(str(int(qualmask1[n,nn]))+'\n')
                            pix2.append(str(n)+' '+str(nn)+' '+str(ID)+'\n')
                        elif int(qualmask1[n,nn]) == 3:
                            threecomps.append(toPrint)
                            threecompsL.append(str(int(qualmask1[n,nn]))+'\n')
                            pix3.append(str(n)+' '+str(nn)+' '+str(ID)+'\n')
                        
                        

            else:
                if nanmask[n,nn] == -5 or nanmask[n,nn] == -4 or nanmask[n,nn] == -3:
                    hello = 0
                else:
                    outfile.write(toPrint)

if train == 0:
    outfile.close()

if train == 1:

    from random import shuffle

    shuf1 = []
    shuf2 = []
    shuf3 = []
    shuf4 = []
    shuf5 = []
    shuf1L = []
    shuf2L = []
    shuf3L = []
    shuf4L = []
    shuf5L = []
    shuf1P = []
    shuf2P = []
    shuf3P = []
    shuf4P = []
    shuf5P = []

    index_shuf = range(len(onecomps))
    shuffle(index_shuf)
    for i in index_shuf:
        shuf1.append(onecomps[i])
        shuf1L.append(onecompsL[i])
        shuf1P.append(pix1[i])

    index_shuf = range(len(twocomps))
    shuffle(index_shuf)
    for i in index_shuf:
        shuf2.append(twocomps[i])
        shuf2L.append(twocompsL[i])
        shuf2P.append(pix2[i])

    index_shuf = range(len(threecomps))
    shuffle(index_shuf)
    for i in index_shuf:
        shuf3.append(threecomps[i])
        shuf3L.append(threecompsL[i])
        shuf3P.append(pix3[i])

    index_shuf = range(len(zerocomps))
    shuffle(index_shuf)
    for i in index_shuf:
        shuf4.append(zerocomps[i])
        shuf4L.append(zerocompsL[i])
        shuf4P.append(pix0[i])

    index_shuf = range(len(badcomps))
    shuffle(index_shuf)
    for i in index_shuf:
        shuf5.append(badcomps[i])
        shuf5L.append(badcompsL[i])
        shuf5P.append(pixb[i])


    print 'ones',np.shape(onecomps)
    print 'twos',np.shape(twocomps)
    print 'threes',np.shape(threecomps)
    print 'bads',np.shape(badcomps)
    print 'zeros',np.shape(zerocomps)
    #stop
    s70 = 2*(np.shape(zerocomps)[0]/4.0)
    s7bad = 0
    s71 = 2*(np.shape(onecomps)[0]/4.0)
    s72 = 2*(np.shape(twocomps)[0]/4.0)
    s73 = 2*(np.shape(threecomps)[0]/4.0)

    #s70 = 0
    #s7bad = 0
    #s71 = 450
    #s72 = 450
    #s73 = 338

    outL = []
    outLL = []
    cvL = []
    cvLL = []
    testL = []
    testLL = []

    cvpix = []
    outpix = []
    testpix = []

    cnts = -1
    for row in shuf1:
        cnts = cnts + 1
        if cnts <= s71:
            outL.append(row)
            outLL.append(shuf1L[cnts])
            outpix.append(shuf1P[cnts])
        elif cnts <= (s71+0.5*s71):
            cvL.append(row)
            cvLL.append(shuf1L[cnts])
            cvpix.append(shuf1P[cnts])
        else:
            testL.append(row)
            testLL.append(shuf1L[cnts])
            testpix.append(shuf1P[cnts])
        

    cnts = -1
    for row in shuf2:
        cnts = cnts + 1
        if cnts <= s72:
            outL.append(row)
            outLL.append(shuf2L[cnts])
            outpix.append(shuf2P[cnts])
        elif cnts <= (s72+0.5*s72):
            cvL.append(row)
            cvLL.append(shuf2L[cnts])
            cvpix.append(shuf2P[cnts])
        else:
            testL.append(row)
            testLL.append(shuf2L[cnts])
            testpix.append(shuf2P[cnts])
        

    cnts = -1
    for row in shuf3:
        cnts = cnts + 1
        if cnts <= s73:
            outL.append(row)
            outLL.append(shuf3L[cnts])
            outpix.append(shuf3P[cnts])
        elif cnts <= (s73+0.5*s73):
            cvL.append(row)
            cvLL.append(shuf3L[cnts])
            cvpix.append(shuf3P[cnts])
        else:
            testL.append(row)
            testLL.append(shuf3L[cnts])
            testpix.append(shuf3P[cnts])
        

    cnts = -1
    for row in shuf4:
        cnts = cnts + 1
        if cnts <= s70:
            outL.append(row)
            outLL.append(shuf4L[cnts])
            outpix.append(shuf4P[cnts])
        elif cnts <= (s70+0.5*s70):
            cvL.append(row)
            cvLL.append(shuf4L[cnts])
            cvpix.append(shuf4P[cnts])
        else:
            testL.append(row)
            testLL.append(shuf4L[cnts])
            testpix.append(shuf4P[cnts])
        

    cnts = -1
    for row in shuf5:
        cnts = cnts + 1
        if cnts <= s7bad:
            outL.append(row)
            outLL.append(shuf5L[cnts])
            outpix.append(shuf5P[cnts])
        elif cnts <= (s7bad+0.5*s7bad):
            cvL.append(row)
            cvLL.append(shuf5L[cnts])
            cvpix.append(shuf5P[cnts])
        else:
            testL.append(row)
            testLL.append(shuf5L[cnts])
            testpix.append(pshuf5P[cnts])
        



    out1_shuf = []
    out2_shuf = []
    out3_shuf = []
    index_shuf = range(len(outL))
    shuffle(index_shuf)
    for i in index_shuf:
        out1_shuf.append(outL[i])
        out2_shuf.append(outLL[i])
        out3_shuf.append(outpix[i])

    cv1_shuf = []
    cv2_shuf = []
    cv3_shuf = []
    index_shuf = range(len(cvL))
    shuffle(index_shuf)
    for i in index_shuf:
        cv1_shuf.append(cvL[i])
        cv2_shuf.append(cvLL[i])
        cv3_shuf.append(cvpix[i])

    test1_shuf = []
    test2_shuf = []
    test3_shuf = []
    index_shuf = range(len(testL))
    shuffle(index_shuf)
    for i in index_shuf:
        test1_shuf.append(testL[i])
        test2_shuf.append(testLL[i])
        test3_shuf.append(testpix[i])



    for row in range(len(out1_shuf)):
        outfile.write(out1_shuf[row])
        outfile2.write(out2_shuf[row])
        outfile3.write(out3_shuf[row])

    for row in range(len(cv1_shuf)):
        cvfile.write(cv1_shuf[row])
        cvfileout.write(cv2_shuf[row])
        cvfilepix.write(cv3_shuf[row])

    for row in range(len(test1_shuf)):
        testfile.write(test1_shuf[row])
        testfileout.write(test2_shuf[row])
        testfilepix.write(test3_shuf[row])
