#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 21:03:33 2022

@author: Max
"""


import numpy as np
import pandas as pd
import matplotlib as mpl
import re
mpl.rcParams.update({'mathtext.default': 'regular'})
mpl.rcParams.update({'lines.markeredgecolor': 'black'})
mpl.rcParams.update({'lines.markersize': 10})
#mpl.rcParams.update({'font.sans-serif': ""})
#mpl.rcParams.update({'font.family': 'sans-serif'})
#mpl.rcParams.update({'mathtext.fontset': 'Arial'})

from matplotlib.legend_handler import HandlerTuple
import matplotlib.pyplot as plt

import os
plt.close('all')
cwd = os.getcwd()

# generate rayleight data
n = int(2e2)
f = np.linspace(1-1e-8, 0, num=n, endpoint=False)
R13_VPDB = 0.011228
p = {'d13C_init': 0,
           'epsilon': -25}
p['alpha'] = p['epsilon']/1000 + 1
p['R13_init'] = (p['d13C_init']/1000 + 1)*R13_VPDB
d = pd.DataFrame(data={'F': f})
d['1-F'] = 1-d['F']
d['R13_rf'] = p['R13_init']*(d['F']**(p['alpha']-1))
d['R13_pf'] = (p['R13_init'] - d['F']*d['R13_rf'])/(1-d['F'])
d['d13C_rf'] = (d['R13_rf']/R13_VPDB-1)*1000
d['d13C_pf'] = (d['R13_pf']/R13_VPDB-1)*1000

cmap = mpl.cm.get_cmap('cool', n)
upper_sigma = 2
normalizer = mpl.colors.Normalize(vmin=p['d13C_init']+p['epsilon'],
                                  vmax=p['d13C_init']-upper_sigma*p['epsilon'])
d['rf_normed'] = normalizer(d['d13C_rf'])
d['pf_normed'] = normalizer(d['d13C_pf'])
d['color_rf'] = cmap(d['rf_normed']).tolist()
d['color_pf'] = cmap(d['pf_normed']).tolist()

# create plot
os.chdir(os.path.join(cwd, 'slides'))
fig, [ax0, ax1] = plt.subplots(nrows=2, figsize=(5,6))
ax0.axis('off')
# draw the Reactant box
# get bounds of bin
bin_floor = 0.1
bin_top = 0.85
rx, ry = ([0.05, 0.05, 0.4, 0.4], [bin_top, bin_floor, bin_floor, bin_top])
rbin = mpl.lines.Line2D(rx, ry, lw=4., color='grey')
ax0.add_line(rbin)

# draw Product box
px, py = ([0.6, 0.6, 0.95, 0.95], [bin_top, bin_floor, bin_floor, bin_top])
pbin = mpl.lines.Line2D(px, py, lw=4., color='grey')
ax0.add_line(pbin)

bin_height = bin_top-bin_floor
# draw an arrow between them
arr = mpl.patches.Arrow(0.4, 0.5, 0.2, 0, width=0.2, color='grey')
ax0.add_patch(arr)

ax1.set_xlim(0, 1)
ax1.set_ylim(d.loc[:,['d13C_rf', 'd13C_pf']].values.min(),p['d13C_init'] - upper_sigma*p['epsilon'])
ax1.set_xlabel('Fraction reacted')
ax1.set_ylabel(r'$\delta^{13}C$ (‰)')
ax0.text(rx[0], bin_top*1.05, "Reactant (r)")
ax0.text(px[0], bin_top*1.05, "Product (p)")
# set the first one to get the party going
i=0
r_amnt = ax0.fill_between(rx[1:3],
                          np.ones(2)*bin_floor,
                          np.ones(2)*(d.loc[i, 'F']*bin_height + bin_floor),
                          color=d.loc[i, 'color_rf'], zorder=0)

p_amnt = ax0.fill_between(px[1:3],
                          np.ones(2)*bin_floor,
                          np.ones(2)*((1-d.loc[i, 'F'])*bin_height + bin_floor),
                          color=d.loc[i, 'color_pf'], zorder=0)

ax1.scatter('1-F', 'd13C_rf', s=20, color=d.loc[i, 'color_rf'], data=d.loc[i,:])
ax1.scatter('1-F', 'd13C_pf', s=20, color=d.loc[i, 'color_pf'], data=d.loc[i,:])
ax1.text(d.loc[int(n/100),'1-F'], d.loc[int(n/100),'d13C_rf']+2, 'r')
ax1.text(d.loc[int(n/100),'1-F'], d.loc[int(n/100),'d13C_pf']+2, 'p')
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=normalizer, cmap=cmap), ax=ax1)
cbar.set_label(r'$\delta^{13}C$ (‰)')
plt.tight_layout()
fig.savefig('rayleigh_series_{0}.png'.format(i))

for i in range(1,n):
    # close the last fig
    plt.close()
    # remove the fills
    r_amnt.remove()
    p_amnt.remove()
    # draw the new ones
    r_amnt = ax0.fill_between(rx[1:3],
                              np.ones(2)*bin_floor,
                              np.ones(2)*(d.loc[i, 'F']*bin_height + bin_floor),
                              color=d.loc[i, 'color_rf'], zorder=0)
    
    p_amnt = ax0.fill_between(px[1:3],
                              np.ones(2)*bin_floor,
                              np.ones(2)*((1-d.loc[i, 'F'])*bin_height + bin_floor),
                              color=d.loc[i, 'color_pf'], zorder=0)
    # add the rayleigh
    ax1.scatter('1-F', 'd13C_rf', s=20, color=d.loc[i, 'color_rf'], data=d.loc[i,:])
    ax1.scatter('1-F', 'd13C_pf', s=20, color=d.loc[i, 'color_pf'], data=d.loc[i,:])
    fig.savefig('rayleigh_series_{0}.png'.format(i))

