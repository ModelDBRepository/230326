# -*- coding: utf-8 -*-
"""
(C) Asaph Zylbertal 01.07.17, HUJI, Jerusalem, Israel
Based on: Hay, E., Hill, S., Schürmann, F., Markram, H., and Segev, I. (2011). Models of Neocortical Layer 5b Pyramidal
Cells Capturing a Wide Range of Dendritic and Perisomatic Active Properties. PLoS Comput. Biol. 7, e1002107. doi:10.1371/journal.pcbi.1002107.

Run a modified pyramidal cell model with a train of synaptic inputs in the dendritic hot zone, concurent with evoked calcium spikes
(article Fig 5E)

****************

"""


import neuron
import pyramidal
import numpy as np
import pickle
import matplotlib.pyplot as plt


def get_site(dist):
    sl = neuron.h.L5PC.locateSites('apic', dist)
    maxdiam = 0
    for i in range(len(sl)):
        dd1 = sl[i][1]
        sec = int(sl[i][0])
        dd = neuron.h.L5PC.apic[sec](dd1).diam
        if dd > maxdiam:
            j = i
            maxdiam = dd

    sec = int(sl[j][0])
    seg = sl[j][1]
    return neuron.h.L5PC.apic[sec](seg)


##########
syn_tau_onset = 0.5
syn_tau_offset = 10.
syn_gmax = 0.005
site_hot = 700.
DNa_coeff_dend = 0.1
soma_stim_amp = 3.0
#########

f = open('params', 'r')
params = pickle.load(f)
f.close()

py = pyramidal.pyramidal(
    params,
    space_factor=3,
    rest_file='rest_state',
    DNa_coeff_dend=DNa_coeff_dend)
py.convert_mechs()
syn_type = neuron.h.naSyn
stim_seg = get_site(site_hot).sec(0.5)

syns = []
soma_stims = []

for i in range(20):
    syns.append(syn_type(stim_seg))
    syns[-1].onset = 2000. + i * 200.
    syns[-1].tau_onset = syn_tau_onset
    syns[-1].tau_offset = syn_tau_offset
    syns[-1].gmax = syn_gmax

    if i < 8:
        soma_stims.append(neuron.h.IClamp(neuron.h.L5PC.soma[0](0.5)))
        soma_stims[-1].dur = 35
        soma_stims[-1].delay = 2020. + i * 200.
        soma_stims[-1].amp = soma_stim_amp

syns.append(syn_type(stim_seg))
syns[-1].onset = 7000.
syns[-1].tau_onset = syn_tau_onset
syns[-1].tau_offset = syn_tau_offset
syns[-1].gmax = syn_gmax

syns.append(syn_type(stim_seg))
syns[-1].onset = 15000.
syns[-1].tau_onset = syn_tau_onset
syns[-1].tau_offset = syn_tau_offset
syns[-1].gmax = syn_gmax

res = {}
res['t'] = neuron.h.Vector()
res['t'].record(neuron.h._ref_t)

res['v_soma'] = neuron.h.Vector()
res['v_soma'].record(neuron.h.L5PC.soma[0](0.5)._ref_v)
res['cai_dend'] = neuron.h.Vector()
res['cai_dend'].record(stim_seg._ref_cai)
res['nai_dend'] = neuron.h.Vector()
res['nai_dend'].record(stim_seg._ref_nai)
res['v_dend'] = neuron.h.Vector()
res['v_dend'].record(stim_seg._ref_v)
res['sec_len'] = stim_seg.sec.L

py.run_model(30000.)

plt.figure()
plt.subplot(311)
plt.plot(res['t'], res['v_dend'], label='Dendrite')
plt.plot(res['t'], res['v_soma'], label='Soma')
plt.legend()
plt.ylabel('Vm (mV)')
plt.subplot(312)
plt.plot(res['t'], res['nai_dend'])
plt.ylabel('Dend [Na+]i (mM)')
plt.subplot(313)
plt.plot(res['t'], 10e3 * np.array(res['cai_dend']))
plt.ylabel('Dend [Ca2+]i (uM)')
