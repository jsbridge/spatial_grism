import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
from glob import glob
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage.filters import gaussian_filter
import matplotlib.patches as patches
from matplotlib.patches import Circle
from scipy import stats
from sklearn.utils.extmath import weighted_mode
import os
import scipy.ndimage as ndi
from astropy.modeling import models
import photutils
import statmorph
from mpl_toolkits.axes_grid1 import make_axes_locatable
import image_diagnostics
reload(image_diagnostics)
from scipy import stats
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
from astropy.convolution import convolve, Box2DKernel


def read_goods_cats():

    # Read the zout catalogs (contains masses and photometry)
    gn_zcat = 'CLEAR_v2.1/grizli_v2.1_cats/goodsn_3dhst.v4.4.zout.fits'
    gs_zcat = 'CLEAR_v2.1/grizli_v2.1_cats/goodss_3dhst.v4.4.zout.fits'
    
    gn = fits.open(gn_zcat)
    gs = fits.open(gs_zcat)
    
    gn_ids = gn[1].data['id']
    gs_ids = gs[1].data['id']
    gn_mass = gn[1].data['mass']
    gs_mass = gs[1].data['mass']
    gn_sfr = gn[1].data['SFR']
    gs_sfr = gs[1].data['SFR']
    gn_U = gn[1].data['restU']
    gs_U = gs[1].data['restU']
    gn_Uerr = gn[1].data['restU_err']
    gs_Uerr = gs[1].data['restU_err']
    gn_B = gn[1].data['restB']
    gs_B = gs[1].data['restB']
    gn_Berr = gn[1].data['restB_err']
    gs_Berr = gs[1].data['restB_err']
    gn_V = gn[1].data['restV']
    gs_V = gs[1].data['restV']
    gn_Verr = gn[1].data['restV_err']
    gs_Verr = gs[1].data['restV_err']
    gn_J = gn[1].data['restJ']
    gs_J = gs[1].data['restJ']
    gn_Jerr = gn[1].data['restJ_err']
    gs_Jerr = gs[1].data['restJ_err']

    vdwn_ids, gn_r, gn_rerr = np.loadtxt('sizes/goodsn/goodsn_3dhst.v4.1_f125w.galfit', unpack=True, usecols=(0,6,7))
    vdws_ids, gs_r, gs_rerr = np.loadtxt('sizes/goodss/goodss_3dhst.v4.1_f125w.galfit', unpack=True, usecols=(0,6,7))

    f = open('goodsn_mass_color_size.dat', 'w')
    f.write('# Field    IDs       Field_ID     mass    SFR       RF_U       RF_U_err        RF_B       RF_B_err        '
                'RF_V       RF_V_err        RF_J       RF_J_err     R      R_err    \n')

    for i, g in enumerate(gn_ids):
        if g != vdwn_ids[i]:
            print g, vdwn_ids
        f.write('GN         '+str(g)+'     GN'+'_'+str(g)+'   '+str(gn_mass[i])+'   '+str(gn_sfr[i])+'   '
                +'   '+str(gn_U[i])+'   '+str(gn_Uerr[i])+'   '+str(gn_B[i])+'   '+str(gn_Berr[i])+'   '
                +str(gn_V[i])+'   '+str(gn_Verr[i])+'   '+str(gn_J[i])+'   '+str(gn_Jerr[i])+'   '+str(gn_r[i])
                +'   '+str(gn_rerr[i])+'\n')

    f.close()
    
    f = open('goodss_mass_color_size.dat', 'w')
    f.write('# Field    IDs      Field_ID     mass    SFR       RF_U       RF_U_err        RF_B       RF_B_err        '
                'RF_V       RF_V_err        RF_J       RF_J_err      R      R_err  \n')
    
    for i, g in enumerate(gs_ids):
        if g != vdws_ids[i]:
            print g, vdws_ids
        f.write('GS         '+str(g)+'     GS'+'_'+str(g)+'   '+'  '+str(gs_mass[i])+'   '+str(gs_sfr[i])+'   '
                +str(gs_U[i])+'   '+str(gs_Uerr[i])+'   '+str(gs_B[i])+'   '+str(gs_Berr[i])+'   '+str(gs_V[i])
                +'   '+str(gs_Verr[i])+'   '+str(gs_J[i])+'   '+str(gs_Jerr[i])+'   '+str(gs_r[i])+'   '+str(gs_rerr[i])+'\n')

    f.close()

    return


def all_CLEAR_color_mass():

    goodsn = np.genfromtxt('goodsn_mass_color_size.dat', unpack=True, dtype=[('f0', 'S9'), ('f1', '<i8'), ('f12', 'S9'), ('f2', '<f8'), ('f3', '<f8'),
                                                                        ('f4', '<f8'), ('f5', '<f8'), ('f6', '<f8'), ('f7', '<f8'),
                                                                        ('f8', '<f8'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'),
                                                                        ('f13', '<f8'), ('f14', '<f8')])
    field_n, gids_n, field_id_n, mass_n, sfr_n, U_n, U_err_n, B_n, B_err_n, V_n, V_err_n, J_n, J_err_n = [],[],[],[],[],[],[],[],[],[],[],[],[]
    R_n, R_err_n = [],[]
    for i, d in enumerate(goodsn):
        field_n.append(d[0])
        gids_n.append(d[1])
        field_id_n.append(d[2])
        mass_n.append(d[3])
        sfr_n.append(d[4])
        U_n.append(d[5])
        U_err_n.append(d[6])
        B_n.append(d[7])
        B_err_n.append(d[8])
        V_n.append(d[9])
        V_err_n.append(d[10])
        J_n.append(d[11])
        J_err_n.append(d[12])
        R_n.append(d[13])
        R_err_n.append(d[14])

    goodss = np.genfromtxt('goodss_mass_color_size.dat', unpack=True, dtype=[('f0', 'S9'), ('f1', '<i8'), ('f12', 'S9'), ('f2', '<f8'), ('f3', '<f8'),
                                                                        ('f4', '<f8'), ('f5', '<f8'), ('f6', '<f8'), ('f7', '<f8'),
                                                                        ('f8', '<f8'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'),
                                                                        ('f13', '<f8'), ('f14', '<f8')])
    field_s, gids_s, field_id_s, mass_s, sfr_s, U_s, U_err_s, B_s, B_err_s, V_s, V_err_s, J_s, J_err_s = [],[],[],[],[],[],[],[],[],[],[],[],[]
    R_s, R_err_s = [],[]
    for i, d in enumerate(goodss):
        field_s.append(d[0])
        gids_s.append(d[1])
        field_id_s.append(d[2])
        mass_s.append(d[3])
        sfr_s.append(d[4])
        U_s.append(d[5])
        U_err_s.append(d[6])
        B_s.append(d[7])
        B_err_s.append(d[8])
        V_s.append(d[9])
        V_err_s.append(d[10])
        J_s.append(d[11])
        J_err_s.append(d[12])
        R_s.append(d[13])
        R_err_s.append(d[14])
    
    cats = glob('CLEAR_v2.1/grizli_v2.1_cats/G*')

    f = open('all_CLEAR_color_mass.dat', 'w')
    f.write('# Field_IDs    IDs     mass    SFR    RF_U       RF_U_err        RF_B       RF_B_err        '
                'RF_V       RF_V_err        RF_J       RF_J_err     R      R_err\n')

    for cat in cats:

        hdu = fits.open(cat)
        field = (cat.split('/')[-1]).split('_')[0]
        bigfield = field[0:2]
    
        ids = hdu[1].data['ID']
        for i,num in enumerate(ids):
            Cfield_id = field+'_'+str(num)
            bigfield_id = bigfield+'_'+str(num)
 
            if bigfield_id in field_id_s:
                loc = [i for i,x in enumerate(field_id_s) if bigfield_id == x][0]
                f.write(Cfield_id+'   '+str(num)+'   '+str(mass_s[loc])+'   '+str(sfr_s[loc])+'   '+str(U_s[loc])+'   '+str(U_err_s[loc])
                        +'   '+str(B_s[loc])+'   '+str(B_err_s[loc])+'   '+str(V_s[loc])+'   '+str(V_err_s[loc])+'   '+str(J_s[loc])
                        +'   '+str(J_err_s[loc])+'   '+str(R_s[loc])+'   '+str(R_err_s[loc])+'\n')
            if bigfield_id in field_id_n:
                loc = [i for i,x in enumerate(field_id_n) if bigfield_id == x][0]
                f.write(Cfield_id+'   '+str(num)+'   '+str(mass_n[loc])+'   '+str(sfr_n[loc])+'   '+str(U_n[loc])+'   '+str(U_err_n[loc])
                        +'   '+str(B_n[loc])+'   '+str(B_err_n[loc])+'   '+str(V_n[loc])+'   '+str(V_err_n[loc])+'   '+str(J_n[loc])
                        +'   '+str(J_err_n[loc])+'   '+str(R_n[loc])+'   '+str(R_err_n[loc])+'\n')

    f.close()
    return
    

def read_CLEAR_cats(SNcut = 3/np.sqrt(2), zcut=1.3):

    # Read the line catalogs (not the other catalogs with mass and RF colors)
    cats = glob('CLEAR_v2.1/grizli_v2.1_cats/G*')
    f = open('ids_z_SN_3rt2ratio.dat', 'w')
    f.write('# Field   IDs        RA             dec             z              OIII               OIII_err          Hbeta      Hbeta_err \n')
    
    for cat in cats:

        hdu = fits.open(cat)
        field = (cat.split('/')[-1]).split('_')[0]
    
        ids = hdu[1].data['ID']
        z = hdu[1].data['z_50']
        RA = hdu[1].data['RA']
        dec = hdu[1].data['DEC']
        OIII = hdu[1].data['OIII_FLUX']
        OIII_e = hdu[1].data['OIII_FLUX_ERR']
        Hb = hdu[1].data['Hb_FLUX']
        Hb_e = hdu[1].data['Hb_FLUX_ERR']

        signal = OIII/Hb
        noise = signal * np.sqrt((OIII_e/OIII)**2 + (Hb_e/Hb)**2)
        SN = signal/noise
        
        # Make a S/N cut on OIII
        gid = [d for i, d in enumerate(ids) if SN[i] > SNcut]
        gz = [d for i, d in enumerate(z) if SN[i] > SNcut]
        gra = [d for i, d in enumerate(RA) if SN[i] > SNcut]
        gdec = [d for i, d in enumerate(dec) if SN[i] > SNcut]
        goiii = [d for i, d in enumerate(OIII) if SN[i] > SNcut]
        goiii_e = [d for i, d in enumerate(OIII_e) if SN[i] > SNcut]
        ghb = [d for i, d in enumerate(Hb) if SN[i] > SNcut]
        ghb_e = [d for i, d in enumerate(Hb_e) if SN[i] > SNcut]

        # Make a redshift cut so I'm only lookin at G102 data
        gid = [d for i, d in enumerate(gid) if gz[i] < zcut]
        gra = [d for i, d in enumerate(gra) if gz[i] < zcut]
        gdec = [d for i, d in enumerate(gdec) if gz[i] < zcut]
        goiii = [d for i, d in enumerate(goiii) if gz[i] < zcut]
        goiii_e = [d for i, d in enumerate(goiii_e) if gz[i] < zcut]
        ghb = [d for i, d in enumerate(ghb) if gz[i] < zcut]
        ghb_e = [d for i, d in enumerate(ghb_e) if gz[i] < zcut]
        gz = [d for i, d in enumerate(gz) if gz[i] < zcut]

        # Make sure Hbeta is measured
        gid = [d for i, d in enumerate(gid) if ghb[i] != -99.0]
        gra = [d for i, d in enumerate(gra) if ghb[i] != -99.0]
        gdec = [d for i, d in enumerate(gdec) if ghb[i] != -99.0]
        goiii = [d for i, d in enumerate(goiii) if ghb[i] != -99.0]
        goiii_e = [d for i, d in enumerate(goiii_e) if ghb[i] != -99.0]
        ghb_e = [d for i, d in enumerate(ghb_e) if ghb[i] != -99.0]
        gz = [d for i, d in enumerate(gz) if ghb[i] != -99.0]
        ghb = [d for i, d in enumerate(ghb) if ghb[i] != -99.0]
        
        for i, g in enumerate(gid):
            print field, g
            f.write(field+'      '+str(g)+'   '+str(np.round(gra[i],9))+'   '+str(np.round(gdec[i],9))+'   '+
                        str(np.round(gz[i],8))+'   '+str(goiii[i])+'   '+str(goiii_e[i])+'   '+
                        str(ghb[i])+'   '+str(ghb_e[i])+'\n')

    f.close()
        
    return


def collate_CLEAR_goods():

    CLEAR = np.genfromtxt('ids_z_SN_3rt2ratio.dat', unpack=True, dtype=[('f0', 'S9'), ('f1', '<i8'), ('f2', '<f8'), ('f3', '<f8'),
                                                                     ('f4', '<f8'), ('f5', '<f8'),('f6', '<f8'), ('f7', '<f8'),
                                                                     ('f8', '<f8')])
    reg, ids, RA, dec, z, OIII, OIII_err, Hb, Hb_err = [],[],[],[],[],[],[],[],[]
    for i, d in enumerate(CLEAR):
        reg.append(d[0])
        ids.append(d[1])
        RA.append(d[2])
        dec.append(d[3])
        z.append(d[4])
        OIII.append(d[5])
        OIII_err.append(d[6])
        Hb.append(d[7])
        Hb_err.append(d[8])

    goodsn = np.genfromtxt('goodsn_mass_color_size.dat', unpack=True, dtype=[('f0', 'S9'), ('f1', '<i8'), ('f12', 'S9'), ('f2', '<f8'), ('f3', '<f8'),
                                                                        ('f4', '<f8'), ('f5', '<f8'), ('f6', '<f8'), ('f7', '<f8'),
                                                                        ('f8', '<f8'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'),
                                                                        ('f13', '<f8'), ('f14', '<f8')])
    field_n, gids_n, field_id_n, mass_n, sfr_n, U_n, U_err_n, B_n, B_err_n, V_n, V_err_n, J_n, J_err_n = [],[],[],[],[],[],[],[],[],[],[],[],[]
    R_n, R_err_n = [], []
    for i, d in enumerate(goodsn):
        field_n.append(d[0])
        gids_n.append(d[1])
        field_id_n.append(d[2])
        mass_n.append(d[3])
        sfr_n.append(d[4])
        U_n.append(d[5])
        U_err_n.append(d[6])
        B_n.append(d[7])
        B_err_n.append(d[8])
        V_n.append(d[9])
        V_err_n.append(d[10])
        J_n.append(d[11])
        J_err_n.append(d[12])
        R_n.append(d[13])
        R_err_n.append(d[14])

    goodss = np.genfromtxt('goodss_mass_color_size.dat', unpack=True, dtype=[('f0', 'S9'), ('f1', '<i8'), ('f12', 'S9'), ('f2', '<f8'), ('f3', '<f8'),
                                                                        ('f4', '<f8'), ('f5', '<f8'), ('f6', '<f8'), ('f7', '<f8'),
                                                                        ('f8', '<f8'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'),
                                                                        ('f13', '<f8'), ('f14', '<f8')])
    field_s, gids_s, field_id_s, mass_s, sfr_s, U_s, U_err_s, B_s, B_err_s, V_s, V_err_s, J_s, J_err_s = [],[],[],[],[],[],[],[],[],[],[],[],[]
    R_s, R_err_s = [], []
    for i, d in enumerate(goodss):
        field_s.append(d[0])
        gids_s.append(d[1])
        field_id_s.append(d[2])
        mass_s.append(d[3])
        sfr_s.append(d[4])
        U_s.append(d[5])
        U_err_s.append(d[6])
        B_s.append(d[7])
        B_err_s.append(d[8])
        V_s.append(d[9])
        V_err_s.append(d[10])
        J_s.append(d[11])
        J_err_s.append(d[12])
        R_s.append(d[13])
        R_err_s.append(d[14])

    f = open('collated_data_SN3rt2.dat', 'w')
    f.write('# Field    IDs        RA             dec                 z              OIII               OIII_err           '
            'Hbeta                Hbeta_err             mass                 sfr                   U                   U_err'
            '    B                     B_err                V                      V_err                J                 J_err'
            '    R            R_err\n')

    for i, r in enumerate(reg):
        print r, ids[i]
        if 'GN' in r:
            if ids[i] in gids_n:
                ind = [j for j,x in enumerate(gids_n) if ids[i] == x][0]
                f.write(r+'       '+str(ids[i])+'    '+str(RA[i])+'    '+str(dec[i])+'     '+str(z[i])+'     ')
                f.write(str(OIII[i])+'    '+str(OIII_err[i])+'   '+str(Hb[i])+'    '+str(Hb_err[i])+'     ')
                f.write(str(np.log10(mass_n[ind]))+'   '+str(sfr_n[ind])+'    '+str(U_n[ind])+'   '+str(U_err_n[ind])+'   ')
                f.write(str(B_n[ind])+'   '+str(B_err_n[ind])+'   '+str(V_n[ind])+'   '+str(V_err_n[ind])+'    ')
                f.write(str(J_n[ind])+'   '+str(J_err_n[ind])+'    '+str(R_n[ind])+'   '+str(R_err_n[ind])+'\n')

        elif 'GS' in r:
            if ids[i] in gids_s:
                ind = [j for j,x in enumerate(gids_s) if ids[i] == x][0]
                f.write(r+'       '+str(ids[i])+'    '+str(RA[i])+'    '+str(dec[i])+'     '+str(z[i])+'     ')
                f.write(str(OIII[i])+'    '+str(OIII_err[i])+'   '+str(Hb[i])+'    '+str(Hb_err[i])+'     ')
                f.write(str(np.log10(mass_s[ind]))+'   '+str(sfr_s[ind])+'    '+str(U_s[ind])+'   '+str(U_err_s[ind])+'   ')
                f.write(str(B_s[ind])+'   '+str(B_err_s[ind])+'   '+str(V_s[ind])+'   '+str(V_err_s[ind])+'    ')
                f.write(str(J_s[ind])+'   '+str(J_err_s[ind])+'    '+str(R_s[ind])+'   '+str(R_err_s[ind])+'\n')

    f.close()

    return


def plot():

    f = open('ratios_grad_SN3.dat', 'w')
    f.write('#  ID        num        RA             dec            z              OIII               OIII_err            Hbeta    '
            '         Hbeta_err          mass                sfr                  U                   U_err                    '
            'B                     B_err                V                V_err                J                 J_err           '
            'in_ratio      out_ratio    loggradient    R     R_err \n')
    
    data = np.genfromtxt('collated_data_SN3.dat', unpack=True, dtype=[('f0', 'S9'), ('f1', '<i8'), ('f2', '<f8'), ('f3', '<f8'),
                                                                      ('f4', '<f8'), ('f5', '<f8'), ('f6', '<f8'), ('f7', '<f8'),
                                                                      ('f8', '<f8'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'),
                                                                      ('f12', '<f8'), ('f13', '<f8'), ('f14', '<f8'), ('f15', '<f8'),
                                                                      ('f16', '<f8'), ('f17', '<f8'), ('f18', '<f8'), ('f19', '<f8'),
                                                                      ('f20', '<f8')])
    gal, ids, ra, dec, red, OIIIline, OIIIline_err, Hbline, Hbline_err, mass, sfr = [],[],[],[],[],[],[],[],[],[],[]
    U, U_err, B, B_err, V, V_err, J, J_err, R, R_err = [],[],[],[],[],[],[],[],[],[]

    for d in data:
        gal.append(d[0])
        ids.append(d[1])
        ra.append(d[2])
        dec.append(d[3])
        red.append(d[4])
        OIIIline.append(d[5])
        OIIIline_err.append(d[6])
        Hbline.append(d[7])
        Hbline_err.append(d[8])
        mass.append(d[9])
        sfr.append(d[10])
        U.append(d[11])
        U_err.append(d[12])
        B.append(d[13])
        B_err.append(d[14])
        V.append(d[15])
        V_err.append(d[16])
        J.append(d[17])
        J_err.append(d[18])
        R.append(d[19])
        R_err.append(d[20])

    reg, ss = [],[]
    for i, r in enumerate(gal):
        if len(str(ids[i])) == 5:
            reg.append(r+'_'+str(ids[i]))
        else:
            reg.append(r+'_0'+str(ids[i]))

    files = glob('CLEAR_v2.1/fits/*/*')

    for fil in files:
        
        name = fil.split('/')[-1]
        field_id = (name.split('.'))[0]
        name = (name.split('.'))[0].split('_')
        field = name[0]
        num = int(name[1])

        if field_id in reg:
            print field_id
            
            ind = [i for i,x in enumerate(reg) if field_id == x]
            z = red[ind[0]]
            size = R[ind[0]]
            
            hdu = fits.open(fil)
            try:
                img = hdu['DSCI'].data[25:55, 25:55]
                OIII_lim = hdu['LINE', 'OIII'].data
                Hb_lim = hdu['LINE', 'Hb'].data
                OIII = hdu['LINE', 'OIII'].data[25:55, 25:55]
                Hb = hdu['LINE', 'Hb'].data[25:55,25:55]
                OIII_err = hdu['LINEWHT', 'OIII'].data[25:55, 25:55]     # These line weights are the inverse 
                Hb_err = hdu['LINEWHT', 'Hb'].data[25:55, 25:55]         # variance = (1/sigma^2), according to Gabe B.

            except KeyError:
                print fil+' generated a key error'
                continue

            os.system('cp '+fil+' cutouts/')
            
            sigma_OIII = np.sqrt(1/OIII_err)
            sigma_Hb = np.sqrt(1/Hb_err)

            #OIII = convolve(OIII, Box2DKernel(2))
            #Hb = cOIII = convolve(Hb, Box2DKernel(2))
            
            ratio = OIII/Hb
            sigma_ratio = np.abs(ratio) * np.sqrt((sigma_OIII/OIII)**2 + (sigma_Hb/Hb)**2)
            ratio_wht = 1/sigma_ratio**2
            
            scale = 0.1     # arcsec/pix, from the full.fits DSCI header
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            scale = cosmo.kpc_proper_per_arcmin(z).value
            dist = np.round((8 * 0.1/60 * scale), 1)
            arcs = 1/0.1
            y = 10 * 60/(scale * 0.1)   # number of pixels that equals 10 kpc
            size = size/0.1

            img_center = [ratio.shape[0]/2., ratio.shape[1]/2.]
            img_cen_range = img[int(img_center[0])-4:int(img_center[0])+4,int(img_center[1])-4:int(img_center[1])+4]
            OIII_cen_range = OIII[int(img_center[0])-4:int(img_center[0])+4,int(img_center[1])-4:int(img_center[1])+4]
            center_OIII = np.where(OIII == np.nanmax(OIII_cen_range))
            center = np.where(img == np.nanmax(img_cen_range))
            ss.append(size)
            #print center_OIII
            #print center
            
            if size <= 2.5:
                idx = [[center[0]+2, center[1]-2],
                [center[0]+2, center[1]-1],
                [center[0]+2, center[1]],
                [center[0]+2, center[1]+1],
                [center[0]+2, center[1]+2],
                [center[0]+1, center[1]-2],
                [center[0]+1, center[1]-1],
                [center[0]+1, center[1]+1],
                [center[0]+1, center[1]+2],
                [center[0], center[1]-2],
                [center[0], center[1]+2],
                [center[0]-2, center[1]-2],
                [center[0]-2, center[1]-1],
                [center[0]-2, center[1]],
                [center[0]-2, center[1]+1],
                [center[0]-2, center[1]+2],
                [center[0]-1, center[1]-2],
                [center[0]-1, center[1]-1],
                [center[0]-1, center[1]+1],
                [center[0]-1, center[1]+2]]

            if size > 2.5:
                idx = [[center[0]+2, center[1]-2],
                [center[0]+2, center[1]-1],
                [center[0]+2, center[1]],
                [center[0]+2, center[1]+1],
                [center[0]+2, center[1]+2],
                [center[0]+1, center[1]-2],
                #[center[0]+1, center[1]-1],
                #[center[0]+1, center[1]+1],
                [center[0]+1, center[1]+2],
                [center[0], center[1]-2],
                [center[0], center[1]+2],
                [center[0]-2, center[1]-2],
                [center[0]-2, center[1]-1],
                [center[0]-2, center[1]],
                [center[0]-2, center[1]+1],
                [center[0]-2, center[1]+2],
                [center[0]-1, center[1]-2],
                #[center[0]-1, center[1]-1],
                #[center[0]-1, center[1]+1],
                [center[0]-1, center[1]+2],
                [center[0]+3, center[1]-1],
                [center[0]+3, center[1]],
                [center[0]+3, center[1]+1],
                [center[0]-3, center[1]-1],
                [center[0]-3, center[1]],
                [center[0]-3, center[1]+1],
                [center[0]+1, center[1]-3],
                [center[0], center[1]-3],
                [center[0]-1, center[1]-3],
                [center[0]+1, center[1]+3],
                [center[0], center[1]+3],
                [center[0]-1, center[1]+3]]

            fig = plt.figure(figsize=(10,3))
            ax1 = fig.add_subplot(142)
            ax2 = fig.add_subplot(143)
            ax3 = fig.add_subplot(144)
            ax4 = fig.add_subplot(141)

            cax4 = ax4.imshow(img, origin='lower', vmax = 4*np.nanstd(img))
            cb4 = fig.colorbar(cax4, ax=ax4, orientation='horizontal', pad=0.03)
            ax4.set_xticks([])
            ax4.set_yticks([])
            ax4.errorbar([22], [2], xerr=[y/2.], color='white', elinewidth=2, mew=2.5)
            ax4.errorbar([center[1][0]+size/2], [center[0][0]], xerr=[size/2], color='m')
            ax4.text(18, 3, '10 kpc', color='white', fontweight='bold')
            ax4.text(2, 25, 'F105W', color='white', fontweight='bold', fontsize=12)
            ax4.add_artist(Circle((img_center[0],img_center[1]), y/2., clip_on=False, zorder=10,linewidth=0.8,
                                edgecolor='white', facecolor='none', linestyle='--'))
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='r',facecolor='none')
                ax4.add_patch(rect)
            ax4.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))
            #ax4.add_patch(patches.Rectangle((0,10),1,1, linewidth = 0.5, edgecolor='r',facecolor='none'))

            ratio_temp = np.zeros(ratio.shape)
            #print OIII/OIII_err
            for i, rr in enumerate(OIII/sigma_OIII):
                for j, r in enumerate(rr):
                    if r < 2 or ~np.isfinite(r) == True:
                        continue
                    else:
                        ratio_temp[i,j] = r

            cax1 = ax1.imshow(OIII, origin='lower', vmin=0, vmax=2*np.std(OIII_lim[60:80, 0:20]))
            ax1.add_artist(Circle((img_center[0],img_center[1]), y/2., clip_on=False, zorder=10,linewidth=0.8,
                                edgecolor='white', facecolor='none', linestyle='--'))
            cb1 = fig.colorbar(cax1, ax=ax1, orientation='horizontal', pad=0.03)
            cb1.set_label(r'10$^{-17}$ erg/s/cm$^2$')
            ax1.text(2, 25, '[OIII]', color='white', fontweight='bold', fontsize=12)
            ax1.set_xticks([])
            ax1.set_yticks([])
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='r',facecolor='none')
                ax1.add_patch(rect)
            ax1.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))
            
            cax2 = ax2.imshow(Hb, origin='lower', vmin=0, vmax=2*np.std(Hb_lim[60:80, 0:20]))
            cb2 = fig.colorbar(cax2, ax=ax2, orientation='horizontal', pad=0.03)
            cb2.set_label(r'10$^{-17}$ erg/s/cm$^2$')
            ax2.text(2, 25, r'H$\beta$', color='white', fontweight='bold', fontsize=12)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.add_artist(Circle((img_center[0],img_center[1]), y/2., clip_on=False, zorder=10,linewidth=0.8,
                                edgecolor='white', facecolor='none', linestyle='--'))
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='r',facecolor='none')
                ax2.add_patch(rect)
            ax2.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))

            ratio_temp = np.zeros(ratio.shape)
            for i, rr in enumerate(ratio/ratio_wht):
                for j, r in enumerate(rr):
                    if r < 2 or ~np.isfinite(r) == True:
                        continue
                    else:
                        ratio_temp[i,j] = r

            cax3 = ax3.imshow(ratio, origin='lower', vmin=-2, vmax=8, cmap='magma')
            cb3 = fig.colorbar(cax3, ax=ax3, orientation='horizontal', pad=0.03)
            cb3.set_label(r'[OIII]/H$\beta$')
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.add_artist(Circle((img_center[0],img_center[1]), y/2., clip_on=False, zorder=10,linewidth=0.8,
                                edgecolor='white', facecolor='none', linestyle='--'))
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='r',facecolor='none')
                ax3.add_patch(rect)
            ax3.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))
                        
            fig.tight_layout()
            plt.savefig('ratio_plots/'+field_id+'_lines_ratio.pdf', dpi=300)
            plt.close(fig)
                    
            in_ratio = (ratio[center[0],center[1]])[0]
            print 'inner ratio:', in_ratio

            ot, otw = [],[]
            for i in idx:
                #if Hb[i] <= 0:
                #    continue
                if np.isnan(ratio_wht[i][0]) == True:
                    continue
                if np.isnan(ratio[i][0]) == True:
                    continue
                ot.append(ratio[i][0])
                otw.append(ratio_wht[i][0])

            #ot = np.array(ot)
            #otw = np.array(otw)
            #out_ratio = np.nansum(ot * otw)/np.nansum(otw)

            out_ratio = weighted_quantile(ot, quantiles=[0.5], sample_weight=otw)[0]
            #out_ratio = weighted_mode(ot, otw)[0][0]
            
            print 'outer ratio:', out_ratio
            grad = np.log10(in_ratio) - np.log10(out_ratio)
            
            f.write(field_id+'   '+str(num)+'   '+str(ra[ind[0]])+'   '+str(dec[ind[0]])+'   '+str(z)+'   '+str(OIIIline[ind[0]])+'   '
                    +str(OIIIline_err[ind[0]])+'   '+str(Hbline[ind[0]])+'   '+str(Hbline_err[ind[0]])+'   '+str(mass[ind[0]])+'   '
                    +str(sfr[ind[0]])+'   '+str(U[ind[0]])+'   '+str(U_err[ind[0]]) +'   '+str(B[ind[0]])+'   '+str(B_err[ind[0]])+'   '
                    +str(V[ind[0]])+'   '+str(V_err[ind[0]])+'   '+str(J[ind[0]])+'   '+str(J_err[ind[0]]) +'   '+str(in_ratio)+'    '
                    +str(out_ratio)+'    '+str(grad)+'    '+str(R[ind[0]])+'    '+str(R_err[ind[0]])+'\n')
            
            #in_err = np.sqrt(np.nansum(np.array(in_wt) * np.array(in_sig)**2)/np.nansum(np.array(in_wt)))
            #out_err = np.sqrt(np.nansum(np.array(out_wt) * np.array(out_sig)**2)/np.nansum(np.array(out_wt)))

            #pix_rad = 3
            #r, ratios, whts, sigs = [], [],[],[]
            #for i, n in enumerate(ratio):
            #    for j, el in enumerate(n):
            #        rval=(np.sqrt((i-center[0])**2+(j-center[1])**2))
            #        if rval <3:
            #            r.append(rval)
            #            if Hb[i,j] > 0:
            #                ratios.append(OIII[i,j]/1e-4)
            #                continue
            #            ratios.append(el)
            #            #whts.append(ratio_wht[i,j])
            #            #sigs.append(sigma_ratio[i,j])
            
            #kpc = np.array(r) * 0.1/60 * scale    # pixel to kpc

            #plt.plot(kpc, ratios, 'k.', ms = 5)
            ##plt.errorbar([np.mean(rad_i)*0.13/60*scale], [in_ratio], yerr=[in_err], color='darkmagenta',
            #                 #marker='o', ms = 10, capsize=3, elinewidth=2, mew=2)
            ##plt.errorbar([np.mean(rad_o)*0.13/60*scale], [out_ratio], yerr=[out_err], color='seagreen',
            #                 #marker='o', ms = 10, capsize=3, elinewidth=2, mew=2)
            #plt.ylim(-5, 20)
            #plt.xlim(0,3)
            #plt.xlabel('Radius (kpc)')
            #plt.ylabel(r'[OIII]/H$\beta$')
            #plt.text(0, 17, r'Inner [OIII]/H$\beta$ = '+str(in_ratio))
            #plt.text(0, 15, r'Outer [OIII]/H$\beta$ = '+str(out_ratio))
            #plt.text(0, 13, r'$\Delta$log([OIII]/H$\beta$) = '+str(grad))
            #plt.savefig('ratio_plots/'+field_id+'_radial.pdf', dpi=300)
            #plt.close()

    f.close()

    return 


def probe_chandra():

    cdfn = np.genfromtxt('CDF-N_2Ms.dat', unpack=True, usecols=(1,2,3,4,5,6,-2,-3), dtype=None, skip_header=217)
    cdfn_sc, class_n, pos, lum_n = [],[],[], []
    for line in cdfn:
        h = line[0]
        m = line[1]
        s = line[2]
        d = line[3]
        am = line[4]
        ase = line[5]
        class_n.append(line[6])
        lum_n.append(line[7])
        pos.append(str(h)+'h'+str(m)+'m'+str(s)+'+'+str(d)+'d'+str(am)+'m'+str(ase))
        
    cdfn_sc = SkyCoord(pos)

    sra_deg, sdec_deg = np.loadtxt('CDF-S_7Ms.dat', unpack=True, usecols=(1,2), dtype=None, skiprows=174)
    cdfs_sc = SkyCoord(sra_deg, sdec_deg, unit=(u.deg, u.deg))

    class_s, lum_s = [],[]
    with open('CDF-S_7Ms.dat') as f:
        lines_after = f.readlines()[174:]
       
    for line in lines_after:
        if 'AGN' in line:
            class_s.append('AGN')
            a = line.split()
            idx = [i for i, s in enumerate(a) if s == 'AGN']
            lum_s.append(float(a[idx[0]-1]))
        if 'Galaxy' in line:
            class_s.append('Galaxy')
            a = line.split()
            idx = [i for i, s in enumerate(a) if s == 'Galaxy']
            lum_s.append(float(a[idx[0]-1]))
        if 'Star' in line:
            class_s.append('Star')
            a = line.split()
            idx = [i for i, s in enumerate(a) if s == 'Star']
            lum_s.append(float(a[idx[0]-1]))

    ids, ra, dec = np.loadtxt('ids_z_SN_3ratio.dat', unpack=True, usecols=(1,2,3))

    clear_sc = SkyCoord(ra, dec, unit=(u.deg, u.deg))

    f = open('AGN_class.dat', 'w')
    f.write('# ID    class    Lum \n')
    xray_ids, clas = [], []
    for i, sc in enumerate(clear_sc):
        sep_s = (cdfs_sc.separation(sc)).arcsecond
        sep_n = (cdfn_sc.separation(sc)).arcsecond
        if min(sep_s) < 0.5:
            xray_ids.append(ids[i])
            clas.append(class_s[np.argmin(sep_s)])
            f.write(str(int(ids[i]))+'   '+class_s[np.argmin(sep_s)]+'    '+str(lum_s[np.argmin(sep_s)])+'\n')
        if min(sep_n) < 0.5:
            xray_ids.append(ids[i])
            clas.append(class_n[np.argmin(sep_n)])
            f.write(str(int(ids[i]))+'   '+class_n[np.argmin(sep_n)]+'    '+str(lum_n[np.argmin(sep_n)])+'\n')

    f.close()

    return 


def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def curate_list():

    bad_ids = np.genfromtxt('bad_obj.dat', dtype='S9')

    data = np.genfromtxt('ratios_grad_SN3.dat', unpack=True, dtype=('S9', '<i8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8',
                        '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8'))
    
    reg, ids, ra, dec, red, OIII, OIII_err, Hb, Hb_err, mass, sfr = [],[],[],[],[],[],[],[],[],[],[]
    U, U_err, B, B_err, V, V_err, J, J_err,in_ratio, out_ratio, grad = [],[],[],[],[],[],[],[],[],[],[]
    R, R_err = [],[]
    for d in data:
        reg.append(d[0])
        ids.append(d[1])
        ra.append(d[2])
        dec.append(d[3])
        red.append(d[4])
        OIII.append(d[5])
        OIII_err.append(d[6])
        Hb.append(d[7])
        Hb_err.append(d[8])
        mass.append(d[9])
        sfr.append(d[10])
        U.append(d[11])
        U_err.append(d[12])
        B.append(d[13])
        B_err.append(d[14])
        V.append(d[15])
        V_err.append(d[16])
        J.append(d[17])
        J_err.append(d[18])
        in_ratio.append(d[19])
        out_ratio.append(d[20])
        grad.append(d[21])
        R.append(d[22])
        R_err.append(d[23])
    
    f = open('ratios_grad_curated_SN3.dat', 'w')
    f.write('#  ID        num        RA             dec            z              OIII               OIII_err            Hbeta    '
            '         Hbeta_err          mass                sfr                  U                   U_err                    '
            'B                     B_err                V                V_err                J                 J_err           '
            'in_ratio      out_ratio    loggradient     R      R_err\n')

    for i, field_id in enumerate(reg):
        if field_id in bad_ids:
            os.system('rm cutouts/'+field_id+'.full.fits')
            continue
        else:
            f.write(field_id+'   '+str(ids[i])+'   '+str(ra[i])+'   '+str(dec[i])+'   '+str(red[i])+'   '+str(OIII[i])+'   '
                    +str(OIII_err[i])+'   '+str(Hb[i])+'   '+str(Hb_err[i])+'   '+str(mass[i])+'   '
                    +str(sfr[i])+'   '+str(U[i])+'   '+str(U_err[i]) +'   '+str(B[i])+'   '+str(B_err[i])+'   '
                    +str(V[i])+'   '+str(V_err[i])+'   '+str(J[i])+'   '+str(J_err[i]) +'   '+str(in_ratio[i])+'    '
                    +str(out_ratio[i])+'    '+str(grad[i])+'   '+str(R[i])+'    '+str(R_err[i])+'\n')

    f.close()
    return


def global_MEx():

    data = np.genfromtxt('ratios_grad_curated_SN3rt2.dat', unpack=True, usecols=(0,1,4,5,6,7,8,9,19,20,21),
                             dtype=[('f0', 'S9'), ('f1', '<i8'), ('f2', '<f8'), ('f3', '<f8'),
                                    ('f4', '<f8'), ('f5', '<f8'), ('f6', '<f8'), ('f7', '<f8'),
                                    ('f8', '<f8'),('f9', '<f8'),('f10', '<f8') ])
    
    reg, ids, red, OIII, OIII_err, Hb, Hb_err, mass,inner, outer, grad = [],[],[],[],[],[],[], [],[],[],[]
    for i, d in enumerate(data):
        reg.append(d[0])
        ids.append(d[1])
        red.append(d[2])
        OIII.append(d[3])
        OIII_err.append(d[4])
        Hb.append(d[5])
        Hb_err.append(d[6])
        mass.append(d[7])
        inner.append(d[8])
        outer.append(d[9])
        grad.append(d[10])

    #ids, red, OIII, OIII_err, Hb, Hb_err, mass, inner, outer, grad = np.loadtxt('ratios_grad_curated.dat', unpack=True,
    #                                                                                usecols=(1,4,5,6,7,8,9,19,20,21))
    
    reg = np.array(reg)[np.isfinite(grad)]
    ids = np.array(ids)[np.isfinite(grad)]
    OIII = np.array(OIII)[np.isfinite(grad)]
    OIII_err = np.array(OIII_err)[np.isfinite(grad)]
    Hb = np.array(Hb)[np.isfinite(grad)]
    Hb_err = np.array(Hb_err)[np.isfinite(grad)]
    inner = np.array(inner)[np.isfinite(grad)]
    outer = np.array(outer)[np.isfinite(grad)]
    mass = np.array(mass)[np.isfinite(grad)]
    red = np.array(red)[np.isfinite(grad)]
    grad = np.array(grad)[np.isfinite(grad)]

    reg = np.array(reg)[np.isfinite(mass)]
    ids = np.array(ids)[np.isfinite(mass)]
    OIII = np.array(OIII)[np.isfinite(mass)]
    OIII_err = np.array(OIII_err)[np.isfinite(mass)]
    Hb = np.array(Hb)[np.isfinite(mass)]
    Hb_err = np.array(Hb_err)[np.isfinite(mass)]
    inner = np.array(inner)[np.isfinite(mass)]
    outer = np.array(outer)[np.isfinite(mass)]
    red = np.array(red)[np.isfinite(mass)]
    grad = np.array(grad)[np.isfinite(mass)]
    mass = np.array(mass)[np.isfinite(mass)]
    
    sigma_OIII = np.sqrt(1/np.array(OIII_err))  # errors in catalog are inverse weights
    sigma_Hb = np.sqrt(1/np.array(Hb_err))
    
    tot_ratio = np.array(OIII)/np.array(Hb)
    tot_ratio_err = np.abs(tot_ratio) * np.sqrt((sigma_OIII/np.array(OIII))**2 + (sigma_Hb/np.array(Hb))**2)

    data = np.genfromtxt('AGN_class.dat', unpack=True, dtype=None)
    xray_ids, xray_class = [],[]
    for d in data:
        xray_ids.append(d[0])
        xray_class.append(d[1])

    xray_rat_agn, xray_mass_agn, xray_rat_gal, xray_mass_gal = [],[],[],[]
    for i, x in enumerate(ids):
        if x in xray_ids:
            if xray_class[(np.where(xray_ids == x)[0])[0]] == 'AGN':
                xray_rat_agn.append(tot_ratio[i])
                xray_mass_agn.append(mass[i])
            if  xray_class[(np.where(xray_ids == x)[0])[0]] == 'Galaxy':
                xray_rat_gal.append(tot_ratio[i])
                xray_mass_gal.append(mass[i])

    x1 = np.arange(8,10, 0.01)
    x2 = np.arange(8,9.6, 0.01)
    x3 = np.arange(10.01,12, 0.01)
    x4 = np.arange(9.61,12, 0.01)
    
    # Upper curve J11
    y1 = 0.375/(x1 - 10.5) + 1.14  # x <= 10
    y3 = 410.24 + -109.333*x3 + 9.71731*x3**2 + -0.288244*x3**3 # otherwise

    # Lower curve J11
    y2 = 0.375/(x2 - 10.5) + 1.14  # x <= 9.6
    y4 = 352.066 + -93.8249*x4 + 8.32651*x4**2 + -0.246416*x4**3 # otherwise

    x11 = np.arange(8 ,9.9, 0.01)
    x22 = np.arange(9.91, 12, 0.01)
    x33 = np.arange(9.9, 11.2, 0.01)
    
    # Upper curve J14
    y11 = 0.37/(x11 - 10.5) + 1.1  # x <= 9.9
    y22 = 594.753 + -167.074*x22 +15.6748*x22**2 + -0.491215*x22**3 # otherwise

    # Lower curve J14
    y33 = 800.492 + -217.328*x33 + 19.6431*x33**2 + -0.591349*x33**3 # 9.9 < x < 11.2

    fig = plt.figure(figsize=(10,4))
    ax = fig.add_subplot(121)
    ax1 = fig.add_subplot(122)
    
    ax.plot(x11, y11, 'k--', lw = 2)
    ax.plot(x22, y22, 'k--', lw = 2)
    ax.plot(x33, y33, 'k--', lw = 2)

    f = open('inner_global_moveup.dat', 'w')
    f.write('# ID          num          mass             log(inner)           log(global)\n')
    new_rats, new_error = [],[]
    for k, m in enumerate(mass):
        if m <= 9.9:
            if k == 21:
                kdx = (np.abs(x11 - m)).argmin()
                if np.log10(tot_ratio[k]) < y11[kdx] and np.log10(inner[k]) > y11[kdx]:
                    f.write(reg[k]+'    '+str(ids[k])+'   '+str(m)+'  '+str(np.log10(inner[k]))+'    '+str(np.log10(tot_ratio[k]))+'\n')
                    ax.plot([m], [np.log10(tot_ratio[k])], marker ='o', color='darkslateblue', label='Global ratio',linestyle = 'None')
                    ax.plot([m], [np.log10(inner[k])],  marker ='o', color='mediumseagreen', label='Inner ratio',linestyle = 'None')
                    ax.vlines(m, np.log10(inner)[k], np.log10(tot_ratio)[k], linewidth=0.5)
                    new_error.append(tot_ratio_err[k])
                    new_rats.append(tot_ratio[k])
            else:
                kdx = (np.abs(x11 - m)).argmin()
                if np.log10(tot_ratio[k]) < y11[kdx] and np.log10(inner[k]) > y11[kdx]:
                    f.write(reg[k]+'    '+str(ids[k])+'   '+str(m)+'  '+str(np.log10(inner[k]))+'    '+str(np.log10(tot_ratio[k]))+'\n')
                    ax.plot([m], [np.log10(tot_ratio[k])], marker ='o', color='darkslateblue')
                    ax.plot([m], [np.log10(inner[k])],  marker ='o', color='mediumseagreen')
                    ax.vlines(m, np.log10(inner)[k], np.log10(tot_ratio)[k], linewidth=0.5)
                    new_error.append(tot_ratio_err[k])
                    new_rats.append(tot_ratio[k])
        if m > 9.9:
            kdx = (np.abs(x22 - m)).argmin()
            if np.log10(tot_ratio[k]) < y22[kdx] and np.log10(inner[k]) > y22[kdx]:
                f.write(reg[k]+'    '+str(ids[k])+'   '+str(m)+'  '+str(np.log10(inner[k]))+'    '+str(np.log10(tot_ratio[k]))+'\n')
                ax.plot(m, np.log10(tot_ratio[k]),  marker ='o', color='darkslateblue')
                ax.plot(m, np.log10(inner[k]),  marker ='o', color='mediumseagreen')
                ax.vlines(m, np.log10(inner)[k], np.log10(tot_ratio)[k], linewidth=0.5)
                new_error.append(tot_ratio_err[k])
                new_rats.append(tot_ratio[k])
                
    ax.set_xlim(8,12)
    ax.set_ylim(-0.6,2)
    ax.set_xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    ax.set_ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    #log_ratio_err = np.abs(tot_ratio_err/(tot_ratio*np.log(10)))
    log_ratio_err = np.abs(np.array(new_error)/(np.array(new_rats)*np.log(10)))
    print np.mean(log_ratio_err)
    ax.errorbar([8.2],[-0.29], yerr=[np.mean(log_ratio_err)], color='black', capsize=2, marker='o', ms=3)
    ax.plot(xray_mass_agn, np.log10(xray_rat_agn), 'ks', mfc='none', ms=10, label='X-Ray AGN')
    ax.plot(xray_mass_gal, np.log10(xray_rat_gal), 'kD', mfc='none', ms=10, label='X-Ray Galaxies')
    legend = ax.legend()
    for legend_handle in legend.legendHandles[2:4]:
        legend_handle._legmarker.set_markersize(8)
    #plt.savefig('MEx_inner_above.pdf', dpi=300)
    #plt.close(fig)
    #f.close()

    #fig, ax = plt.subplots(figsize=(5.5,4))
    ax1.plot(x11, y11, 'k--', lw = 2)
    ax1.plot(x22, y22, 'k--', lw = 2)
    ax1.plot(x33, y33, 'k--', lw = 2)

    f = open('inner_global_movedown.dat', 'w')
    f.write('# ID          num          mass             log(inner)           log(global)\n')
    new_rats1, new_error1 = [],[]
    for k, m in enumerate(mass):
        if m <= 9.9:
            if k == 10:
                kdx = (np.abs(x11 - m)).argmin()
                if np.log10(tot_ratio[k]) > y11[kdx] and np.log10(inner[k]) < y11[kdx]:
                    f.write(reg[k]+'    '+str(ids[k])+'   '+str(m)+'  '+str(np.log10(inner[k]))+'    '+str(np.log10(tot_ratio[k]))+'\n')
                    ax1.plot([m], [np.log10(tot_ratio[k])], marker ='o', color='darkslateblue', label='Global ratio',linestyle = 'None')
                    ax1.plot([m], [np.log10(inner[k])],  marker ='o', color='mediumseagreen', label='Inner ratio',linestyle = 'None')
                    ax1.vlines(m, np.log10(inner)[k], np.log10(tot_ratio)[k], linewidth=0.5)
                    new_error1.append(tot_ratio_err[k])
                    new_rats1.append(tot_ratio[k])
            else:
                kdx = (np.abs(x11 - m)).argmin()
                if np.log10(tot_ratio[k]) > y11[kdx] and np.log10(inner[k]) < y11[kdx]:
                    f.write(reg[k]+'    '+str(ids[k])+'   '+str(m)+'  '+str(np.log10(inner[k]))+'    '+str(np.log10(tot_ratio[k]))+'\n')
                    ax1.plot([m], [np.log10(tot_ratio[k])], marker ='o', color='darkslateblue')
                    ax1.plot([m], [np.log10(inner[k])],  marker ='o', color='mediumseagreen')
                    ax1.vlines(m, np.log10(inner)[k], np.log10(tot_ratio)[k], linewidth=0.5)
                    new_error1.append(tot_ratio_err[k])
                    new_rats1.append(tot_ratio[k])
        if m > 9.9:
            kdx = (np.abs(x22 - m)).argmin()
            if np.log10(tot_ratio[k]) > y22[kdx] and np.log10(inner[k]) < y22[kdx]:
                f.write(reg[k]+'    '+str(ids[k])+'   '+str(m)+'  '+str(np.log10(inner[k]))+'    '+str(np.log10(tot_ratio[k]))+'\n')
                ax1.plot(m, np.log10(tot_ratio[k]),  marker ='o', color='darkslateblue')
                ax1.plot(m, np.log10(inner[k]),  marker ='o', color='mediumseagreen')
                ax1.vlines(m, np.log10(inner)[k], np.log10(tot_ratio)[k], linewidth=0.5)
                new_error1.append(tot_ratio_err[k])
                new_rats1.append(tot_ratio[k])
                
    ax1.set_xlim(8,12)
    ax1.set_ylim(-0.6,2)
    ax1.set_xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    #ax1.set_ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    log_ratio_err = np.abs(np.array(new_error1)/(np.array(new_rats1)*np.log(10)))
    print np.mean(log_ratio_err)
    ax1.errorbar([8.2],[-0.29], yerr=[np.mean(log_ratio_err)], color='black', capsize=2, marker='o', ms=3)
    ax1.plot(xray_mass_agn, np.log10(xray_rat_agn), 'ks', mfc='none', ms=10, label='X-Ray AGN')
    ax1.plot(xray_mass_gal, np.log10(xray_rat_gal), 'kD', mfc='none', ms=10, label='X-Ray Galaxies')
    legend1 = ax1.legend()
    for legend_handle in legend1.legendHandles[2:4]:
        legend_handle._legmarker.set_markersize(8)
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    ax1.set_yticklabels([])
    ax.set_xticklabels([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11, 11.5]) 
    plt.savefig('MEx_innervglobal_movers.pdf', dpi=300)
    plt.close(fig)
    f.close()
    
    fig, ax = plt.subplots(figsize=(6.5,4))
    #ax.plot(x1, y1, 'k--', lw = 2)
    #ax.plot(x2, y2, 'k--', lw = 2)
    #ax.plot(x3, y3, 'k--', lw = 2)
    #ax.plot(x4, y4, 'k--', lw = 2)
    ax.plot(x11, y11, 'k--', lw = 2)
    ax.plot(x22, y22, 'k--', lw = 2)
    ax.plot(x33, y33, 'k--', lw = 2)

    log_ratio_err = np.abs(tot_ratio_err/(tot_ratio*np.log(10)))
    ax.plot(xray_mass_agn, np.log10(xray_rat_agn), 'ks', mfc='none', ms=10, label='X-Ray AGN')
    ax.plot(xray_mass_gal, np.log10(xray_rat_gal), 'kD', mfc='none', ms=10, label='X-Ray Galaxies')
    #ax.errorbar(mass, np.log10(tot_ratio), yerr=log_ratio_err, color='black', capsize=2, ls='none', elinewidth=0.5)
    pt = ax.scatter(mass, np.log10(tot_ratio), c = grad, vmin=0, vmax=2, marker='o', s = 15, cmap='viridis')
    ax.errorbar([8.2],[-0.19], yerr=[np.mean(log_ratio_err)], color='black', capsize=2, marker='o', ms=3)
    cb = fig.colorbar(pt)
    cb.set_label(r'$\Delta$log$_{10}$([OIII]/H$\beta$)')
    plt.xlim(8,12)
    plt.ylim(-0.5,1.4)
    plt.xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    plt.ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    legend = plt.legend()
    for legend_handle in legend.legendHandles[0:2]:
        legend_handle._legmarker.set_markersize(8)
    plt.savefig('MEx_cb.pdf', dpi=300)
    plt.close(fig)
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(mass, np.log10(tot_ratio), grad)
    #ax.set_xlabel('mass')
    #ax.set_ylabel('total ratio')
    #ax.set_zlabel('gradient')
    #plt.show()

    #print np.median(grad), np.mean(grad)
    #print np.min(grad)
    #scgrad = np.array(grad)+0.7
    #cutt = np.arange(0.8, 2.7, 0.1)
    #g_up, g_low = [],[]
    #for c in cutt:
    #    sorted_upper = np.sort(scgrad[scgrad>=c])
    #    sorted_lower = np.sort(scgrad[scgrad<c])
    #    num_up = np.arange(1, len(sorted_upper)+1, 1)
    #    num_low = np.arange(1, len(sorted_lower)+1, 1)
    #    n_up = len(sorted_upper)
    #    n_low = len(sorted_lower)
    #    g_up.append((1/(np.mean(sorted_upper)*n_up*(n_up - 1))) * np.sum((2*num_up - n_up - 1) * sorted_upper))
    #    g_low.append((1/(np.mean(sorted_lower)*n_low*(n_low - 1))) * np.sum((2*num_low - n_low - 1) * sorted_lower))


    #kmeans = KMeans(n_clusters=2, random_state=0).fit(zip(mass, np.log10(tot_ratio), grad))
    #kmeans_labels = kmeans.labels_

    #fig = plt.figure(figsize=(10,3))
    
    #ax1 = fig.add_subplot(131)
    #ax2 = fig.add_subplot(132)
    #ax3 = fig.add_subplot(133)
    #ax1.scatter(mass, np.log10(tot_ratio), c = kmeans_labels)
    #ax1.set_xlabel('mass')
    #ax1.set_ylabel('tot ratio')
    #ax2.scatter(mass, grad, c = kmeans_labels)
    #ax2.set_xlabel('mass')
    #ax2.set_ylabel('gradient')
    #ax3.scatter(grad, np.log10(tot_ratio), c = kmeans_labels)
    #ax3.set_xlabel('gradient')
    #ax3.set_ylabel('tot ratio')
    #plt.tight_layout()
    #plt.savefig('kmeans.pdf', dpi=300)
    #plt.close(fig)
    #a = [g for i, g in enumerate(grad) if kmeans_labels[i]==1]
    #print max(a), min(a)
    #b  = [g for i, g in enumerate(grad) if kmeans_labels[i]==0]
    #print max(b), min(b)

    cut = 0.75
    
    fig, ax = plt.subplots(figsize=(6,9))
    #ax.plot(x1, y1, 'k--', lw = 2)
    #ax.plot(x2, y2, 'k--', lw = 2)
    #ax.plot(x3, y3, 'k--', lw = 2)
    #ax.plot(x4, y4, 'k--', lw = 2)
    ax.plot(x11, y11, 'k--', lw = 2)
    ax.plot(x22, y22, 'k--', lw = 2)
    ax.plot(x33, y33, 'k--', lw = 2)

    log_ratio_err = np.abs(tot_ratio_err/(tot_ratio*np.log(10)))
    ax.plot(xray_mass_agn, np.log10(xray_rat_agn), 'ks', mfc='none', ms=10, label='X-Ray AGN')
    ax.plot(xray_mass_gal, np.log10(xray_rat_gal), 'kD', mfc='none', ms=10, label='X-Ray Galaxies')
    ax.scatter(mass[grad>=cut], np.log10(tot_ratio)[grad>=cut], marker='o',
                   s = 15, color='cornflowerblue', label=r'$\Delta$log$_{10}$([OIII]/H$\beta$) > '+str(cut))
    ax.scatter(mass[grad<cut], np.log10(tot_ratio)[grad<cut], marker='o',
                   s = 15, color='darkblue', label=r'$\Delta$log$_{10}$([OIII]/H$\beta$) < '+str(cut))
    ax.errorbar([8.2],[-0.19], yerr=[np.mean(log_ratio_err)], color='black', capsize=2, marker='o', ms=3)

    #ax.set_xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    ax.set_ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    ax.legend(ncol=2, fontsize=8)
    ax.xaxis.set_tick_params(labelbottom=False)

    legend = ax.legend(frameon=True, ncol=2, fontsize=8)
    for legend_handle in legend.legendHandles[0:2]:
        legend_handle._legmarker.set_markersize(7)

    divider = make_axes_locatable(ax)
    
    ax2 = divider.append_axes("bottom", 3, pad=0.1, sharex=ax)
    ax2.set_xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    ax2.set_ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    ax2.set_xlim(8,12)
    ax2.set_ylim(-0.4, 1.25)
    #ax2.plot(x1, y1, 'k--', lw = 2)
    #ax2.plot(x2, y2, 'k--', lw = 2)
    #ax2.plot(x3, y3, 'k--', lw = 2)
    #ax2.plot(x4, y4, 'k--', lw = 2)
    ax2.plot(x11, y11, 'k--', lw = 2)
    ax2.plot(x22, y22, 'k--', lw = 2)
    ax2.plot(x33, y33, 'k--', lw = 2)
    ax2.set_xlim(8,12)
    ax2.set_ylim(-0.5, 1.45)

    xx = mass[grad<cut]
    yy = np.log10(tot_ratio)[grad<cut]
    xmin = np.min(xx)
    xmax = np.max(xx)
    ymin = np.min(yy)
    ymax = np.max(yy)
    
    X, Y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([xx, yy])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

    print np.min(Z), np.max(Z)
    lev = [0.27, 0.32, 0.38, 0.45, 0.51, 0.57, 0.63, 0.69, 0.75, 0.82, 0.88,0.905]
    ax2.contour(X, Y, Z, levels =lev, colors='darkblue')
    
    xx = mass[grad>=cut]
    yy = np.log10(tot_ratio)[grad>=cut]
    xmin = np.min(xx)
    xmax = np.max(xx)
    ymin = np.min(yy)
    ymax = np.max(yy)
    
    X, Y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([xx, yy])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

    print np.min(Z), np.max(Z)
    lev =  np.arange(np.min(Z), np.max(Z), 0.1)[2:] + 0.09
    lev = [0.295,0.37, 0.45, 0.53, 0.6, 0.69, 0.79, 0.89, 0.95, 0.98]
    ax2.contour(X, Y, Z, levels=lev, colors='cornflowerblue')
    
    axHistx = divider.append_axes("top", 0.6, pad=0.1, sharex=ax)
    axHisty = divider.append_axes("right", 0.6, pad=0.1, sharey=ax)
    axHistx.xaxis.set_tick_params(labelbottom=False)
    axHisty.yaxis.set_tick_params(labelleft=False)

    axHistx.hist(mass[grad<cut], bins=15, histtype="stepfilled", color='darkblue', edgecolor='cornflowerblue', alpha=0.6)
    axHistx.hist(mass[grad>=cut], bins =12, histtype="stepfilled",color = 'cornflowerblue', edgecolor='darkblue', alpha=0.5)
    axHisty.hist(np.log10(tot_ratio)[grad<cut], bins =10, histtype="stepfilled", color='darkblue',
                     edgecolor='cornflowerblue', orientation='horizontal', alpha=0.7)
    axHisty.hist(np.log10(tot_ratio)[grad>=cut], bins =10, histtype="stepfilled", color='cornflowerblue',
                     edgecolor='darkblue', orientation='horizontal', alpha=0.4)

    axHistx.set_xlim(8,12)
    axHisty.set_ylim(-0.5, 1.45)
    
    plt.savefig('MEx_cut.pdf', dpi=300)
    plt.close(fig)

    
    fig, ax = plt.subplots()
    #plt.plot(x1, y1, 'k--', lw = 2)
    #plt.plot(x2, y2, 'k--', lw = 2)
    #plt.plot(x3, y3, 'k--', lw = 2)
    #plt.plot(x4, y4, 'k--', lw = 2)
    ax.plot(x11, y11, 'k--', lw = 2)
    ax.plot(x22, y22, 'k--', lw = 2)
    ax.plot(x33, y33, 'k--', lw = 2)

    for i, m in enumerate(mass):
        if i ==0:
            plt.plot([m], [np.log10(inner)[i]], 'o', color='navy', ms = 3, label='Inner')
            plt.plot([m], [np.log10(tot_ratio)[i]], 'o', color = 'orchid', ms = 3, label='Global')
        else:
            plt.plot([m], [np.log10(inner)[i]], 'o', color='navy', ms = 3)
            plt.plot([m], [np.log10(tot_ratio)[i]], 'o', color = 'orchid', ms = 3)
        ax.vlines(m, np.log10(inner)[i], np.log10(tot_ratio)[i], linewidth=0.5)
        
    #plt.plot(mass, np.log10(outer), 'D', color = 'orchid', ms = 3, label=r'Outer')

    plt.xlim(8, 12)
    plt.xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    plt.ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    lgnd = plt.legend()
    for legend_handle in lgnd.legendHandles:
        legend_handle._legmarker.set_markersize(6)
    plt.savefig('MEx_inner_global.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    #plt.plot(x1, y1, 'k--', lw = 2)
    #plt.plot(x2, y2, 'k--', lw = 2)
    #plt.plot(x3, y3, 'k--', lw = 2)
    #plt.plot(x4, y4, 'k--', lw = 2)
    ax.plot(x11, y11, 'k--', lw = 2)
    ax.plot(x22, y22, 'k--', lw = 2)
    ax.plot(x33, y33, 'k--', lw = 2)
    
    for i, m in enumerate(mass):
        if i ==0:
            plt.plot([m], [np.log10(inner)[i]], 'o', color='navy', ms = 3, label='Inner')
            plt.plot([m], [np.log10(outer)[i]], 'o', color = 'cornflowerblue', ms = 3, label='Outer')
        else:
            plt.plot([m], [np.log10(inner)[i]], 'o', color='navy', ms = 3)
            plt.plot([m], [np.log10(outer)[i]], 'o', color = 'cornflowerblue', ms = 3)
        ax.vlines(m, np.log10(inner)[i], np.log10(outer)[i], linewidth=0.5)
        
    #plt.plot(mass, np.log10(outer), 'D', color = 'orchid', ms = 3, label=r'Outer')

    plt.xlim(8, 12)
    plt.xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    plt.ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    lgnd = plt.legend()
    for legend_handle in lgnd.legendHandles:
        legend_handle._legmarker.set_markersize(6)
    plt.savefig('MEx_inner_outer.pdf', dpi=300)
    plt.close(fig)

    return 


def basic_plots():

    ids, red, OIII, OIII_err, Hb, Hb_err, mass, sfr = np.loadtxt('ratios_grad_curated_SN3rt2.dat', unpack=True, usecols=(1,4,5,6,7,8,9,10))
    U, U_err, B, B_err, V, V_err, J, J_err = np.loadtxt('ratios_grad_curated_SN3rt2.dat', unpack=True, usecols=(11,12,13,14,15,16,17,18))
    inner, outer, grad, R, R_err = np.loadtxt('ratios_grad_curated_SN3rt2.dat', unpack=True, usecols=(19,20,21, 22,23))

    sigma_OIII = np.sqrt(1/np.array(OIII_err))  # errors in catalog are inverse weights
    sigma_Hb = np.sqrt(1/np.array(Hb_err))
    
    ratio = np.array(OIII)/np.array(Hb)
    ratio_err = np.abs(ratio) * np.sqrt((sigma_OIII/np.array(OIII))**2 + (sigma_Hb/np.array(Hb))**2)

    movers = np.loadtxt('inner_global_moveup.dat', unpack=True, usecols=(1))

    idx_mov = []
    for i, n in enumerate(movers):
        idx_mov.append(np.where(ids == n)[0][0])
    
    cut = 0.75

    #grad_cut = np.array(grad)[grad > cut]
    #ids_cut = np.array(ids)[grad > cut]
    #inner_cut = np.array(inner)[grad > cut]
    #outer_cut = np.array(outer)[grad > cut]
    #red_cut = np.array(red)[grad > cut]
    #mass_cut = np.array(mass)[grad > cut]
    #sfr_cut = np.array(sfr)[grad > cut]
    #U_cut = np.array(U)[grad > cut]
    #U_err_cut = np.array(U_err)[grad > cut]
    #B_cut = np.array(B)[grad > cut]
    #B_err_cut = np.array(B_err)[grad > cut]
    #V_cut = np.array(V)[grad > cut]
    #V_err_cut = np.array(V_err)[grad > cut]
    #J_cut = np.array(J)[grad > cut]
    #J_err_cut = np.array(J_err)[grad > cut]
    #OIII_cut = np.array(OIII)[grad > cut]
    #OIII_err_cut = np.array(OIII_err)[grad > cut]
    #R_cut = np.array(R)[grad > cut]
    #R_err_cut = np.array(R_err)[grad > cut]

    grad_cut = np.array(grad)[idx_mov]
    ids_cut = np.array(ids)[idx_mov]
    inner_cut = np.array(inner)[idx_mov]
    outer_cut = np.array(outer)[idx_mov]
    red_cut = np.array(red)[idx_mov]
    mass_cut = np.array(mass)[idx_mov]
    sfr_cut = np.array(sfr)[idx_mov]
    U_cut = np.array(U)[idx_mov]
    U_err_cut = np.array(U_err)[idx_mov]
    B_cut = np.array(B)[idx_mov]
    B_err_cut = np.array(B_err)[idx_mov]
    V_cut = np.array(V)[idx_mov]
    V_err_cut = np.array(V_err)[idx_mov]
    J_cut = np.array(J)[idx_mov]
    J_err_cut = np.array(J_err)[idx_mov]
    OIII_cut = np.array(OIII)[idx_mov]
    OIII_err_cut = np.array(OIII_err)[idx_mov]
    R_cut = np.array(R)[idx_mov]
    R_err_cut = np.array(R_err)[idx_mov]

    data = np.genfromtxt('AGN_class.dat', unpack=True, dtype=None)
    xray_ids, xray_class, xray_lum = [],[],[]
    for d in data:
        xray_ids.append(d[0])
        xray_class.append(d[1])
        xray_lum.append(d[2])

    xray_rat_agn, xray_mass_agn, xray_rat_gal, xray_mass_gal, = [],[],[],[],
    xray_OIII_agn, xray_OIII_gal, xray_z_gal, xray_z_agn = [],[],[],[]
    xray_grad_agn, xray_grad_gal, xray_U_agn, xray_U_gal = [],[],[],[]
    xray_lum_gal, xray_lum_agn, xray_agn_ids, xray_gal_ids = [],[],[],[]
    xray_B_agn, xray_B_gal = [],[]
    for i, x in enumerate(ids):
        if x in xray_ids:
            if xray_class[(np.where(xray_ids == x)[0])[0]] == 'AGN':
                xray_agn_ids.append(x)
                xray_rat_agn.append(ratio[i])
                xray_mass_agn.append(mass[i])
                xray_OIII_agn.append(OIII[i])
                xray_z_agn.append(red[i])
                xray_grad_agn.append(grad[i])
                xray_B_agn.append(B[i])
                xray_U_agn.append(U[i])
                xray_lum_agn.append(xray_lum[(np.where(xray_ids == x)[0])[0]])
            if  xray_class[(np.where(xray_ids == x)[0])[0]] == 'Galaxy':
                xray_gal_ids.append(x)
                xray_rat_gal.append(ratio[i])
                xray_mass_gal.append(mass[i])
                xray_OIII_gal.append(OIII[i])
                xray_z_gal.append(red[i])
                xray_grad_gal.append(grad[i])
                xray_B_gal.append(B[i])
                xray_U_gal.append(U[i])
                xray_lum_gal.append(xray_lum[(np.where(xray_ids == x)[0])[0]])
                
    xray_rat_agn_grad, xray_mass_agn_grad, xray_rat_gal_grad, xray_mass_gal_grad,xray_OIII_agn_grad,xray_OIII_gal_grad, xray_z_gal_grad, xray_z_agn_grad = [],[],[],[],[],[],[],[]
    xray_lum_gal_grad, xray_lum_agn_grad = [],[]
    for i, x in enumerate(ids):
        if x in xray_ids and x in ids_cut:
            if xray_class[(np.where(xray_ids == x)[0])[0]] == 'AGN':
                xray_rat_agn_grad.append(ratio[i])
                xray_mass_agn_grad.append(mass[i])
                xray_OIII_agn_grad.append(OIII[i])
                xray_z_agn_grad.append(red[i])
                xray_lum_agn_grad.append(xray_lum[(np.where(xray_ids == x)[0])[0]])
            if  xray_class[(np.where(xray_ids == x)[0])[0]] == 'Galaxy':
                xray_rat_gal_grad.append(ratio[i])
                xray_mass_gal_grad.append(mass[i])
                xray_OIII_gal_grad.append(OIII[i])
                xray_z_gal_grad.append(red[i])
                xray_lum_gal_grad.append(xray_lum[(np.where(xray_ids == x)[0])[0]])

    all_mass, gU, gV, gJ = np.loadtxt('all_CLEAR_color_mass.dat', unpack=True, usecols=(2,4,8,10))
    all_mass = np.log10(all_mass)

    fig, ax = plt.subplots()
    ax.plot(R[np.logical_and(R>0, R<2)], grad[np.logical_and(R>0, R<2)], 'bo')
    plt.xlabel(r'Effective Radius (arcseconds)')
    plt.ylabel(r'$\Delta$log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    plt.savefig('Rvgrad.pdf', dpi=300)
    plt.close(fig)
    
    fig, ax = plt.subplots()
    ax.plot(xray_mass_agn, xray_grad_agn, 'ks', mfc='none', ms=10, label='X-Ray AGN')
    ax.plot(xray_mass_gal, xray_grad_gal,'kD', mfc='none', ms=10, label='X-Ray Galaxies')
    ax.plot(mass, grad, color='darkblue',marker='o', ms=4, ls='None')
    ax.hlines(0.5, 7.5, 12, linestyles='dashed')
    plt.xlim(8.1,12.1)
    plt.xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    plt.ylabel(r'$\Delta$log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    plt.legend()
    plt.savefig('gradvmass.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    #ax.hist(all_mass[np.where(all_mass>0)], histtype='step', color='black', fill=False, ls='dashed', bins=50, lw=2, label='All CLEAR')
    ax.hist(mass[np.isfinite(mass)], histtype='step',color='darkblue', bins=20, lw=2, ls='dashed', label=r'All S/N$_{\mathrm{[OIII]/H}\beta} > 3/\sqrt{2}$')
    ax.hist(mass_cut[np.isfinite(mass_cut)], histtype='stepfilled', color='cornflowerblue', bins=15, lw=2)
    ax.set_xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    ax.set_ylabel('Number')
    plt.legend()
    plt.savefig('mass_hist.pdf', dpi=300)
    plt.close(fig)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    dist = cosmo.luminosity_distance(red).value * 3.086e24    # Mpc to cm
    dist_cut = cosmo.luminosity_distance(red_cut).value * 3.086e24
    OIII_lum = np.array(OIII) * 1e-17 * 4 * np.pi * dist**2# /3.848e33               # Divide by L_sol
    OIII_lum_cut = np.array(OIII_cut) * 1e-17 * 4* np.pi * dist_cut**2# /3.848e33
    Hb_lum = np.array(Hb) * 1e-17 * 4 * np.pi * dist**2

    #print r_ids[np.where(np.log10(OIII_lum) ==np.min(np.log10(OIII_lum)))]
    
    fig, ax = plt.subplots()
    ax.hist(np.log10(OIII_lum), histtype='step',color='darkblue', bins=25, fill=False,
                lw=2, ls='dashed', label=r'All S/N$_{\mathrm{[OIII]}}$ > 3')
    ax.hist(np.log10(OIII_lum_cut), histtype='stepfilled', color='cornflowerblue', bins=15, lw=2)#, label=r'$\Delta$log([OIII]/H$\beta$) > '+str(cut))
    ax.set_xlabel(r'$L_{\mathrm{[OIII]}}$ [log$_{10}$(erg/s)]')
    ax.set_ylabel('Number')
    plt.legend(loc='upper left')
    plt.savefig('OIII_hist.pdf', dpi=300)
    plt.close(fig)

    n1, bins1, p = ax.hist(np.sort(np.log10(OIII_lum)), bins = 60)#, histtype='step',color='black', bins=50, fill=False,
                #lw=2, ls='dashed', label=r'All S/N$_{\mathrm{[OIII]}}$ > 10')
    n2, bins2, p = ax.hist(np.sort(np.log10(OIII_lum_cut)), bins = 60, range=(bins1[0], bins1[-1]))
    #print n1, n2
    #print bins1, bins2
    fig, ax = plt.subplots()
    rat = n2/n1
    rat[np.isnan(rat)]=0
    ax.plot(bins1[:-1],rat)
    ax.set_xlabel(r'$L_{\mathrm{[OIII]}}$ [log$_{10}$(erg/s)]')
    ax.set_ylabel('Fractions of AGN Identified')
    
    plt.savefig('OIII_hist2.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(np.log10(OIII), grad, 'bo')
    ax.set_xlabel(r'log OIII Flux')
    plt.ylabel(r'$\Delta$log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    #plt.xlim(40,43)
    plt.savefig('OIIIvgrad.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(np.log10(Hb), grad, 'bo')
    ax.set_xlabel(r'log Hbeta Flux')
    plt.ylabel(r'$\Delta$log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    plt.savefig('Hbvgrad.pdf', dpi=300)
    plt.close(fig)

    dist_agn = cosmo.luminosity_distance(xray_z_agn).value * 3.086e24    # Mpc to cm
    dist_gal = cosmo.luminosity_distance(xray_z_gal).value * 3.086e24
    OIII_lum_agn = np.array(xray_OIII_agn) * 1e-17 * 4 * np.pi * dist_agn**2# /3.848e33               # Divide by L_sol
    OIII_lum_gal = np.array(xray_OIII_gal) * 1e-17 * 4* np.pi * dist_gal**2
    #dist_agn_cut = cosmo.luminosity_distance(xray_z_agn_grad).value * 3.086e24    # Mpc to cm
    #dist_gal_cut = cosmo.luminosity_distance(xray_z_gal_grad).value * 3.086e24
    #OIII_lum_agn_cut = np.array(xray_OIII_agn_grad) * 1e-17 * 4 * np.pi * dist_agn_cut**2# /3.848e33   # Divide by L_sol
    #OIII_lum_gal_cut = np.array(xray_OIII_gal_grad) * 1e-17 * 4* np.pi * dist_gal_cut**2

    fig, ax = plt.subplots()
    ax.plot(np.log10(OIII_lum_agn),np.log10(xray_lum_agn),  'ko', mfc='none',
                label=r'AGN $\Delta$log([OIII]/H$\beta$) < '+str(cut), ms = 11)
    ax.plot( np.log10(OIII_lum_gal),np.log10(xray_lum_gal), 'ks', mfc='none',
                 label=r'Galaxies $\Delta$log([OIII]/H$\beta$) < '+str(cut), ms = 11)
    #ax.plot( np.log10(OIII_lum_agn_cut),np.log10(xray_lum_agn_grad), color='darkmagenta', marker='o',
    #             ls='None', label=r'AGN $\Delta$log([OIII]/H$\beta$) > '+str(cut), ms = 11)
    #ax.plot( np.log10(OIII_lum_gal_cut),np.log10(xray_lum_gal_grad), color='darkmagenta', marker='s',
    #             ls='None', label=r'Galaxies $\Delta$log([OIII]/H$\beta$) > '+str(cut), ms = 11)
    ax.set_xlim(40.9, 43.05)
    ax.set_ylim(40.5, 45)
    ax.plot([40.5, 43.06], [41.2, 43.06], 'k--')
    ax.plot([41.45, 43.06], [40.5, 41.55], 'k:')
    plt.text(41.1, 41.85, 'AGN', rotation=15, fontsize=12)
    plt.text(41.5, 40.67, 'SF', rotation=17, fontsize=12)
    ax.set_ylabel(r'$L_{\mathrm{(0.5-7 keV)}}$ [log$_{10}$(erg/s)]')
    ax.set_xlabel(r'$L_{\mathrm{[OIII]}}$ [log$_{10}$(erg/s)]')
    lgnd = plt.legend(loc='upper left', fontsize=8)
    for legend_handle in lgnd.legendHandles:
        legend_handle._legmarker.set_markersize(8)
    plt.savefig('xraylumvOIIIlum.pdf', dpi=300)
    plt.close(fig)
    
    gU = gU[np.nonzero(gV)]
    gJ = gJ[np.nonzero(gV)]
    gV = gV[np.nonzero(gV)]

    gV = 25. - 2.5*np.log10(gV)
    gJ = 25. - 2.5*np.log10(gJ)
    gU = 25. - 2.5*np.log10(gU)

    xmin = np.min((gV-gJ))
    xmax = np.max((gV-gJ))
    ymin = np.min((gU-gV))
    ymax = np.max((gU-gV))
    
    X, Y = np.mgrid[xmin:xmax:300j, ymin:ymax:300j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([gV-gJ, gU-gV])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

    print np.max(Z)
    print np.min(Z)
    fig, ax = plt.subplots(figsize=(5,5))
    cfset = ax.contourf(X, Y, Z, levels=[0.13,0.2, 0.3,0.4,0.45, 0.5,0.7,0.9,1.2, 1.5, 1.7, 1.9, 2.192178724589], cmap='Blues')
    #cset = ax.contour(X, Y, Z, , colors='k', linewidths = 0.5, label='All CLEAR')
    
    cV = (25. - 2.5*np.log10(V))[np.isfinite(mass)]
    cU = (25. - 2.5*np.log10(U))[np.isfinite(mass)]
    cJ = (25. - 2.5*np.log10(J))[np.isfinite(mass)]
    cB = (25. - 2.5*np.log10(B))[np.isfinite(mass)]

    xmin = np.min((cV-cJ))
    xmax = np.max((cV-cJ))+0.4
    ymin = np.min((cU-cV))
    ymax = np.max((cU-cV))
    
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([cV-cJ, cU-cV])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

    #cfset = ax.contourf(X, Y, Z, cmap='Blues', label=r'S/N$_{\mathrm{[OIII]}}$ > 3')
    cset = ax.contour(X, Y, Z, colors='k', linewidths = 0.5)

    #ax.scatter(gV-gJ, gU-gV, marker='s', color='gray', s=1, label='All CLEAR')

    V = (25. - 2.5*np.log10(V_cut))
    U = (25. - 2.5*np.log10(U_cut))
    J = (25. - 2.5*np.log10(J_cut))
    ax.plot(V-J, U-V, marker='o', color='orchid', ms = 4, ls='None', mec='k', mew=0.5)#, label=r'$\Delta$log([OIII]/H$\beta$) < '+str(cut))
    ax.set_xlim(-0.5,2.3)
    ax.set_ylim(-0.3,2.2)
    ax.set_xlabel('(V - J)')
    ax.set_ylabel('(U - V)')
    plt.legend(loc='lower right')
    plt.savefig('colors.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5,5))
    plt.plot([0.1, 1.25], [1.3, -0.1], 'k-')
    plt.hlines(-0.1, 1.25, 1.7)
    plt.plot(cU-cB, np.log10(inner[np.isfinite(mass)]), 'o', color='navy', mfc = 'white', ms = 5, label=r'Inner')
    plt.plot(cU-cB, np.log10(outer[np.isfinite(mass)]), 'D', color = 'orchid', mfc = 'none', ms = 3, label=r'Outer')
    plt.plot(cU-cB, np.log10(ratio[np.isfinite(mass)]), 's', color='black', ms=2, label='Global')
    plt.legend()
    ax.set_xlabel('(U - B)')
    ax.set_ylabel(r'log([OIII]/H$\beta$)')
    plt.savefig('U-Bvsratio.pdf', dpi=300)
    plt.close(fig)

    xU_gal = (25. - 2.5*np.log10(xray_U_gal))
    xU_agn = (25. - 2.5*np.log10(xray_U_agn))
    xB_gal = (25. - 2.5*np.log10(xray_B_gal))
    xB_agn = (25. - 2.5*np.log10(xray_B_agn))
    
    fig, ax = plt.subplots(figsize=(6,5))
    plt.plot([0, 1.25], [1.35, -0.1], 'k-')
    plt.hlines(-0.1, 1.25, 1.7)
    pt = ax.scatter(cU-cB, np.log10(ratio[np.isfinite(mass)]), c = grad[np.isfinite(mass)], vmin=0, vmax=1.6)
    ax.plot(xU_agn-xB_agn, np.log10(xray_rat_agn), 'ks', mfc='none', ms=10, label='X-Ray AGN')
    ax.plot(xU_gal-xB_gal, np.log10(xray_rat_gal), 'kD', mfc='none', ms=10, label='X-Ray Galaxies')
    plt.legend()
    cb = fig.colorbar(pt)
    cb.set_label(r'$\Delta$log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    ax.set_xlabel('(U - B)')
    ax.set_ylabel(r'log([OIII]/H$\beta$)')
    ax.set_xlim(0, 1.7)
    plt.savefig('U-Bvsratio_grad.pdf', dpi=300)
    plt.close(fig)
    
    return


def morph():

    files = glob('cutouts/*')
    f = open('morphology.dat', 'w')
    f.write('#Field_ID    ID         Gini                    M20                  A           Flag\n')

    for fil in files:

        name = fil.split('/')[-1]
        field_id = (name.split('.'))[0]
        name = (name.split('.'))[0].split('_')
        field = name[0]
        num = int(name[1])

        hdu = fits.open(fil)
        seg = hdu['SEG'].data#[25:55, 25:55]
        img = hdu['DSCI'].data#[25:55, 25:55]
        wht = hdu['DWHT'].data#[25:55, 25:55]

        # Isolate only the segmentation map of the object in question
        for i, ss in enumerate(seg):
            for j, s in enumerate(ss):
                if s == num:
                    continue
                else:
                    seg[i,j]=0
                    
        #plt.imshow(seg, origin='lower', vmin=np.min(seg), vmax=np.max(seg))
        #plt.show()
        
        #threshold = photutils.detect_threshold(img, snr=1.5)
        #npixels = 5  # minimum number of connected pixels
        #segm = photutils.detect_sources(img, threshold, npixels)

        #label = np.argmax(segm.areas) + 1
        #segmap = segm.data == label

        try:
            source_morphs = statmorph.source_morphology(img, seg, weightmap=wht)
        except AssertionError:
            print field_id
            continue
        m = source_morphs[0]

        fig = image_diagnostics.make_figure(m)
        fig.savefig('morphology/'+str(field_id), dpi=300)

        f.write(field_id+'   '+str(num)+'   '+str(m.gini)+'   '+str(m.m20)+'   '+str(m.asymmetry)+'   '+str(m.flag)+'\n')

    f.close()

    return 


def morph_plots():

    data = np.genfromtxt('ratios_grad_curated.dat', unpack=True, usecols=(0,1,4,5,6,7,8,9,21),
                             dtype=[('f0', 'S9'), ('f1', '<i8'), ('f2', '<f8'), ('f3', '<f8'),
                                    ('f4', '<f8'), ('f5', '<f8'), ('f6', '<f8'), ('f7', '<f8'),('f8', '<f8') ])
    
    reg, ids, red, OIII, OIII_err, Hb, Hb_err, mass, grad = [],[],[],[],[],[],[], [],[]
    for i, d in enumerate(data):
        reg.append(d[0])
        ids.append(d[1])
        red.append(d[2])
        OIII.append(d[3])
        OIII_err.append(d[4])
        Hb.append(d[5])
        Hb_err.append(d[6])
        mass.append(d[7])
        grad.append(d[8])

    sigma_OIII = np.sqrt(1/np.array(OIII_err))  # errors in catalog are inverse weights
    sigma_Hb = np.sqrt(1/np.array(Hb_err))
    
    tot_ratio = np.array(OIII)/np.array(Hb)
    tot_ratio_err = np.abs(tot_ratio) * np.sqrt((sigma_OIII/np.array(OIII))**2 + (sigma_Hb/np.array(Hb))**2)

    morph = np.genfromtxt('morphology.dat', unpack=True, usecols=(0,2,3,4),
                          dtype=[('f0', 'S9'), ('f2', '<f8'), ('f3', '<f8'), ('f4', '<f8')])
    regm, gini, m20, asym = [],[],[],[]
    for m in morph:
        regm.append(m[0])
        gini.append(m[1])
        m20.append(m[2])
        asym.append(m[3])

    new_gini, new_m20, new_asym = [],[],[]
    new_mass, new_tot_ratio, new_grad = [],[],[]
    for ii, n in enumerate(regm):
        for j, num in enumerate(reg):
            if n == num:
                new_gini.append(gini[ii])
                new_m20.append(m20[ii])
                new_asym.append(asym[ii])
                new_mass.append(mass[j])
                new_tot_ratio.append(tot_ratio[j])
                new_grad.append(grad[j])

    print len(new_mass), len(new_gini)

    data = np.genfromtxt('AGN_class.dat', unpack=True, dtype=None)
    xray_ids, xray_class = [],[]
    for d in data:
        xray_ids.append(d[0])
        xray_class.append(d[1])

    xray_rat_agn, xray_mass_agn, xray_rat_gal, xray_mass_gal = [],[],[],[]
    for i, x in enumerate(ids):
        if x in xray_ids:
            if xray_class[(np.where(xray_ids == x)[0])[0]] == 'AGN':
                xray_rat_agn.append(tot_ratio[i])
                xray_mass_agn.append(mass[i])
            if  xray_class[(np.where(xray_ids == x)[0])[0]] == 'Galaxy':
                xray_rat_gal.append(tot_ratio[i])
                xray_mass_gal.append(mass[i])

    x1 = np.arange(7,10, 0.01)
    x2 = np.arange(7,9.6, 0.01)
    x3 = np.arange(10.01,13, 0.01)
    x4 = np.arange(9.61,13, 0.01)
    
    # Upper curve
    y1 = 0.375/(x1 - 10.5) + 1.14  # x <= 10
    y3 = 410.24 + -109.333*x3 + 9.71731*x3**2 + -0.288244*x3**3 # otherwise

    # Lower curve
    y2 = 0.375/(x2 - 10.5) + 1.14  # x <= 9.6
    y4 = 352.066 + -93.8249*x4 + 8.32651*x4**2 + -0.246416*x4**3 # otherwise

    fig, ax = plt.subplots(figsize=(6.5,4))
    ax.plot(x1, y1, 'k--', lw = 2)
    ax.plot(x2, y2, 'k--', lw = 2)
    ax.plot(x3, y3, 'k--', lw = 2)
    ax.plot(x4, y4, 'k--', lw = 2)

    log_ratio_err = np.abs(tot_ratio_err/(tot_ratio*np.log(10)))
    #ax.plot(xray_mass_agn, np.log10(xray_rat_agn), 'ks', mfc='none', ms=10, label='X-Ray AGN')
    #ax.plot(xray_mass_gal, np.log10(xray_rat_gal), 'kD', mfc='none', ms=10, label='X-Ray Galaxies')
    pt = ax.scatter(new_mass, np.log10(new_tot_ratio), c = np.array(new_gini), vmin=np.min(new_gini),
                        vmax=0.65, marker='o', s = 15, cmap='viridis')
    ax.errorbar([8.4],[-0.32], yerr=[np.mean(log_ratio_err)], color='black', capsize=2)
    cb = fig.colorbar(pt)
    cb.set_label(r'Gini')
    plt.xlim(8.1,12.1)
    plt.ylim(-0.6,1.5)
    plt.xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    plt.ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    plt.legend()
    plt.savefig('MEx_gini.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5,4))
    ax.plot(x1, y1, 'k--', lw = 2)
    ax.plot(x2, y2, 'k--', lw = 2)
    ax.plot(x3, y3, 'k--', lw = 2)
    ax.plot(x4, y4, 'k--', lw = 2)
    
    pt = ax.scatter(new_mass, np.log10(new_tot_ratio), c = new_m20, vmin=np.min(new_m20),
                        vmax=-0.8, marker='o', s = 15, cmap='viridis')
    ax.errorbar([8.4],[-0.32], yerr=[np.mean(log_ratio_err)], color='black', capsize=2)
    cb = fig.colorbar(pt)
    cb.set_label(r'M20')
    plt.xlim(8.1,12.1)
    plt.ylim(-0.6,1.5)
    plt.xlabel(r'Stellar Mass [log$_{10}$(M$_{\odot}$)]')
    plt.ylabel(r'Log$_{10}$([OIII]$\lambda5007$/H$\beta$)')
    plt.legend()
    plt.savefig('MEx_m20.pdf', dpi=300)
    plt.close(fig)
    
    fig, ax = plt.subplots(figsize=(6.5,4))
    pt = plt.scatter(new_m20, new_gini, c = new_grad, marker='o', s = 15, cmap='viridis')
    cb = fig.colorbar(pt)
    cb.set_label(r'Gradient')
    plt.ylabel('Gini')
    plt.xlabel('M20')
    plt.savefig('gini_m20.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    pt = plt.scatter(new_gini, new_grad, marker='o', s = 15, cmap='viridis')
    plt.xlabel('Gini')
    plt.ylabel('Gradient')
    plt.savefig('gini_grad.pdf', dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots()
    pt = plt.scatter(new_m20, new_grad, marker='o', s = 15, cmap='viridis')
    plt.xlabel('M20')
    plt.ylabel('Gradient')
    plt.savefig('m20_grad.pdf', dpi=300)
    plt.close(fig)
    
    return
