import matplotlib
matplotlib.use("agg")

def test_arange():
    import healpy
    import healtics
    import numpy as np

    healtics.plot_with_ticks(
                np.arange(healpy.nside2npix(8)),
                overplot=[[]]
            )

