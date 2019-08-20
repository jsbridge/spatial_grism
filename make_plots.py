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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.convolution import convolve, Box2DKernel, Gaussian2DKernel


def make_plots(way):

    if way == 'up':
        want = np.genfromtxt('inner_global_moveup.dat', unpack=True, usecols=(0), dtype='S9')
        folder = 'ratio_plots_moveup'
    if way == 'down':
        want = np.genfromtxt('inner_global_movedown.dat', unpack=True, usecols=(0), dtype='S9')
        folder = 'ratio_plots_movedown'

    files = glob('CLEAR_v2.1/fits/*/*')

    for fil in files:
        
        name = fil.split('/')[-1]
        field_id = (name.split('.'))[0]
        name = (name.split('.'))[0].split('_')
        field = name[0]
        num = int(name[1])

        if field_id in want:
            print field_id
            
            hdu = fits.open(fil)
            try:
                img = hdu['DSCI'].data[32:48, 32:48]
                OIII_lim = hdu['LINE', 'OIII'].data
                Hb_lim = hdu['LINE', 'Hb'].data
                OIII = hdu['LINE', 'OIII'].data[32:48, 32:48]
                Hb = hdu['LINE', 'Hb'].data[32:48, 32:48]
                OIII_err = hdu['LINEWHT', 'OIII'].data[32:48, 32:48]     # These line weights are the inverse 
                Hb_err = hdu['LINEWHT', 'Hb'].data[32:48, 32:48]         # variance = (1/sigma^2), according to Gabe B.

            except KeyError:
                print fil+' generated a key error'
                continue
            
            sigma_OIII = np.sqrt(1/OIII_err)
            sigma_Hb = np.sqrt(1/Hb_err)

            reg_OIII = np.copy(OIII)
            reg_Hb = np.copy(Hb)

            OIII = convolve(OIII_lim, Box2DKernel(2))
            Hb = convolve(Hb_lim, Box2DKernel(2))
            OIII = OIII[32:48, 32:48]
            Hb = Hb[32:48, 32:48]
            
            ratio = OIII/Hb
            reg_ratio = reg_OIII/reg_Hb

            sigma_ratio = np.abs(reg_OIII/reg_Hb) * np.sqrt((sigma_OIII/reg_OIII)**2 + (sigma_Hb/reg_Hb)**2)
            ratio_wht = 1/sigma_ratio**2

            img_center = [ratio.shape[0]/2., ratio.shape[1]/2.]
            img_cen_range = img[int(img_center[0])-2:int(img_center[0])+2,int(img_center[1])-4:int(img_center[1])+2]
            OIII_cen_range = OIII[int(img_center[0])-2:int(img_center[0])+2,int(img_center[1])-4:int(img_center[1])+2]
            center_OIII = np.where(OIII == np.nanmax(OIII_cen_range))
            center = np.where(img == np.nanmax(img_cen_range))

            idx = [[center[0]+2, center[1]-2],
            [center[0]+2, center[1]-1],
            [center[0]+2, center[1]],
            [center[0]+2, center[1]+1],
            [center[0]+2, center[1]+2],
            [center[0]+1, center[1]-2],
            [center[0]+1, center[1]+2],
            [center[0], center[1]-2],
            [center[0], center[1]+2],
            [center[0]-2, center[1]-2],
            [center[0]-2, center[1]-1],
            [center[0]-2, center[1]],
            [center[0]-2, center[1]+1],
            [center[0]-2, center[1]+2],
            [center[0]-1, center[1]-2],
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

            y = 15
            
            fig = plt.figure()
            ax1 = fig.add_subplot(132)
            ax2 = fig.add_subplot(133)
            ax4 = fig.add_subplot(131)

            im4 = ax4.imshow(img, origin='lower', vmin = 0, vmax =4*np.nanstd(img))
            ax4_divider = make_axes_locatable(ax4)
            cax4 = ax4_divider.append_axes("bottom", size = '6%', pad='3%')
            cb4 = fig.colorbar(im4, cax=cax4, orientation="horizontal")
            cb4.set_label(r'Counts')
            ax4.set_xticks([])
            ax4.set_yticks([])
            ax4.text(0.5, 13.5, 'F105W', color='white', fontweight='bold', fontsize=12)
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='black',facecolor='none')
                ax4.add_patch(rect)
            ax4.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))

            im1 = ax1.imshow(OIII, origin='lower', vmin=0, vmax=2*np.std(OIII_lim[60:80, 0:20]))
            ax1_divider = make_axes_locatable(ax1)
            cax1 = ax1_divider.append_axes("bottom", size = '6%', pad='3%')
            cb1 = fig.colorbar(im1, cax=cax1, orientation="horizontal")
            cb1.set_label(r'10$^{-17}$ erg/s/cm$^2$')
            ax1.text(0.5, 13.5, '[OIII]', color='white', fontweight='bold', fontsize=12)
            ax1.set_xticks([])
            ax1.set_yticks([])
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='black',facecolor='none')
                ax1.add_patch(rect)
            ax1.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))
            
            im2 = ax2.imshow(Hb, origin='lower', vmin=0, vmax=2*np.std(Hb_lim[60:80, 0:20]))
            ax2_divider = make_axes_locatable(ax2)
            cax2 = ax2_divider.append_axes("bottom", size = '6%', pad='3%')
            cb2 = fig.colorbar(im2, cax=cax2, orientation="horizontal")
            cb2.set_label(r'10$^{-17}$ erg/s/cm$^2$')
            ax2.text(0.5, 13.5, r'H$\beta$', color='white', fontweight='bold', fontsize=12)
            ax2.set_xticks([])
            ax2.set_yticks([])
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='white',facecolor='none')
                ax2.add_patch(rect)
            ax2.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))
            
            fig.tight_layout()
            fig.subplots_adjust(wspace=0.1)
            plt.savefig(folder+'/'+field_id+'_lines.pdf', dpi=300)
            plt.close(fig)
            

            fig = plt.figure(figsize=(6,3))
            ax3 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)

            im3 = ax3.imshow(reg_ratio, origin='lower', vmin=-1, vmax=8, cmap='magma')
            ax3.set_xticks([])
            ax3.set_yticks([])
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='white',facecolor='none')
                ax3.add_patch(rect)
            ax3.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))

            im2 = ax2.imshow(ratio, origin='lower', vmin=-1, vmax=8, cmap='magma')
            ax2.set_xticks([])
            ax2.set_yticks([])
            for r in idx:
                rect = patches.Rectangle((r[1]-0.5, r[0]-0.5), 1, 1, linewidth = 0.5, edgecolor='white',facecolor='none')
                ax2.add_patch(rect)
            ax2.add_patch(patches.Rectangle((center[1]-0.5, center[0]-0.5),1,1, linewidth = 0.5, edgecolor='k',facecolor='none'))

            # Hacky code to get plots the same size - make colorbar for first subplot then delete it
            ax3_divider = make_axes_locatable(ax3)
            cax3 = ax3_divider.append_axes("right", size = '6%', pad='5%')
            cb3 = fig.colorbar(im3, cax = cax3)
            fig.delaxes(fig.axes[2])
            
            ax2_divider = make_axes_locatable(ax2)
            cax2 = ax2_divider.append_axes("right", size = '6%', pad='5%')
            cb2 = fig.colorbar(im2, cax=cax2)
            cb2.set_label(r'[OIII]/H$\beta$')
            
            fig.tight_layout()
            fig.subplots_adjust(wspace=-0.05)
            plt.savefig(folder+'/'+field_id+'_ratio.pdf', dpi=300)
            plt.close(fig)

            
            in_ratio = (ratio[center[0],center[1]])[0]
            print 'inner ratio:', in_ratio

            ot, otw = [],[]
            for i in idx:
                if np.isnan(ratio_wht[i][0]) == True:
                    continue
                if np.isnan(ratio[i][0]) == True:
                    continue
                ot.append(ratio[i][0])
                otw.append(ratio_wht[i][0])

            print ot
            
