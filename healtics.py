from setup_matplotlib import *
import healpy as hp
import scipy
import scipy.stats
import numpy as np

from matplotlib.projections.geo import GeoAxes

class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi

    Shifts labelling from -180,180 to 0-360"""
    def __call__(self, x, pos=None):
        if x != 0:
            x *= -1
        if x < 0:
            x += 2*np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)

#def plot_with_ticks(m,vmin=0,vmax=None,cmap="YlOrBr",cmap2="bone",unit="$erg cm^{-2} s^{-1}$",title="",overplot=None):
def plot_with_ticks(m,vmin=0,vmax=None,cmap="YlOrBr",unit="$10^{-7} \mathrm{erg^{}cm^{-2} s^{-1}}$",title="",overplot=[],ilevels=None,ovmin=0.1,ovmax=1):
    if vmax is None:
        vmax=m.max()    

    xsize = 2000
    ysize = xsize/2.
    nside=hp.npix2nside(m.shape[0])
    
    theta = np.linspace(np.pi, 0, ysize)
    phi   = np.linspace(-np.pi, np.pi, xsize)
    longitude = np.radians(np.linspace(-180, 180, xsize))
    latitude = np.radians(np.linspace(-90, 90, ysize))

    # project the map to a rectangular matrix xsize x ysize
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(nside, THETA, PHI)
    grid_map = m[grid_pix]


    width=18
# for width in [ 8.8]:
# for width in [18., 12., 8.8]:
    fig = plt.figure(figsize=(cm2inch(width), cm2inch(width)/(3./2.)))
    fig.suptitle(title)
    # matplotlib is doing the mollveide projection
    ax = fig.add_subplot(111,projection='mollweide')
    

    # rasterized makes the map bitmap while the labels remain vectorial
    # flip longitude to the astro convention
    image = plt.pcolormesh(longitude[::-1], latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=True, cmap=cmap)

    for op,cmapop,onelevel in overplot:
        grid_map_op=hp.get_interp_val(op,THETA,PHI)
        print "overplotting",grid_map_op



        if onelevel is None:
            #levels=np.logspace(-3,0,5)*grid_map_op.max()
            #levels=[0.5*0.1*grid_map_op.max(),0.5*0.5*grid_map_op.max()]
            levels=[(1-scipy.stats.norm(0,1.).cdf(sigma))*2*grid_map_op.max() for sigma in [1,2,3]]

            level_containment=[]
            for n in np.linspace(1,0,500):
                #l=grid_map_op.max()*10**-n
                l=grid_map_op.max()*n
                containment=grid_map_op[grid_map_op>=l].sum()/grid_map_op.sum()
                #print grid_map_op[grid_map_op>=l].shape[0],grid_map_op.flatten().shape[0]
                area=grid_map_op[grid_map_op>=l].shape[0]*4*np.pi*(180/np.pi)**2/grid_map_op.flatten().shape[0]
                #print "level",l/grid_map_op.max(),"containment",containment,"area",area
                level_containment.append([l,l/grid_map_op.max(),containment])

            get_c=lambda l:sorted(level_containment,key=lambda x:abs(x[2]-l))[0][0]
            levels=[get_c(0.75),get_c(0.95)]

        else:
            levels=[op[op>0].min()*onelevel]

        print "levels:",levels

        if ilevels is not None:
            levels=ilevels
            ovmin=levels.min()
            ovmax=levels.max()

        levels=sorted(levels)

        plt.contour(longitude[::-1], latitude, grid_map_op, vmin=ovmin, vmax=ovmax, rasterized=True, cmap=cmapop,levels=levels)
        #image2 = plt.pcolormesh(longitude[::-1], latitude, grid_map, vmin=0.1, vmax=1., rasterized=True, cmap=cmap,alpha=0.3)


        # project the map to a rectangular matrix xsize x ysize
        

    # graticule
    ax.set_longitude_grid(60)
    ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))
    if width < 10:
        ax.set_latitude_grid(45)
        ax.set_longitude_grid_ends(90)

    # colorbar
    cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05, ticks=[vmin, vmax])
    cb.ax.xaxis.set_label_text(unit)
    cb.ax.xaxis.labelpad = -8
    # workaround for issue with viewers, see colorbar docstring
    cb.solids.set_edgecolor("face")

    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)

    # remove tick labels
#    ax.xaxis.set_ticklabels([])
 #   ax.yaxis.set_ticklabels([])
    # remove grid
  #  ax.xaxis.set_ticks([])
    #ax.yaxis.set_ticks([])

     #remove white space around figure
    spacing = 0.05
    plt.subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing)

    plt.grid(True)
    return fig

