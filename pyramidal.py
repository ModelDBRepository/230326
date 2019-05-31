# -*- coding: utf-8 -*-
"""
(C) Asaph Zylbertal 01.07.17, HUJI, Jerusalem, Israel
Based on: Hay, E., Hill, S., Schürmann, F., Markram, H., and Segev, I. (2011). Models of Neocortical Layer 5b Pyramidal
Cells Capturing a Wide Range of Dendritic and Perisomatic Active Properties. PLoS Comput. Biol. 7, e1002107. doi:10.1371/journal.pcbi.1002107.

Introduction of [Na+]i related mechanisms to the mitral cell model

****************

"""
import neuron
import sys
from neuron import gui
import numpy as np
import matplotlib as mpl
from numpy import mean
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import pickle
import copy


def over_sample_space(factor):
    for sec in neuron.h.allsec():
        sec.nseg *= factor


class pyramidal(object):
    def __init__(
            self,
            params,
            space_factor=1,
            rest_file=None,
            cv=1,
            atol=0.003,
            DNa_coeff_dend=1.0):

        neuron.load_mechanisms('./mod')
        neuron.h.xopen('loadModel4.hoc')
        neuron.h.define_shape()
        over_sample_space(space_factor)
        self.cv = neuron.h.CVode()
        self.cv.active(cv)
        self.fixna = 0
        self.DNa_coeff_dend = DNa_coeff_dend
        self.event_fun = None
        if cv:
            self.cv.atol(atol)
        self.params = params
        self.converted = False
        if rest_file is not None:
            f = open(rest_file, 'r')
            self.rest_vals = pickle.load(f)
            f.close()
            self.fih = neuron.h.FInitializeHandler(1, self.restore_states)

    def events(self, ev_fun):
        self.event_fun = ev_fun

    def run_model(self, run_time, parts=1000):

        neuron.h.finitialize()
        neuron.h.fcurrent()

        if self.cv.active() == 1:
            self.cv.re_init()

        part_len = run_time / parts
        for part in range(parts):

            neuron.run((part + 1) * part_len)
            sys.stdout.write("\r%d / %d" % (part + 1, parts))
            sys.stdout.flush()

    def convert_mechs(self, nadp=True):
        neuron.h.celsius = 35.0
        if nadp:
            self.insert_na_mech()
        self.nadp = nadp
        self.remove_ca_dynamics()
        self.insert_ca_mech()
        neuron.h.use_ghk_Ca_HVA = 1
        neuron.h.use_ghk_Ca_LVAst = 1
        self.converted = True

    def insert_na_mech(self):
        for sec in neuron.h.allsec():
            sec.insert('nadp')
            sec.TotalPump_nadp = self.params['TotalPump_nadp_dend']

        for sec in neuron.h.L5PC.axon:
            sec.TotalPump_nadp = self.params['TotalPump_nadp_axon']

        for sec in neuron.h.L5PC.soma:
            sec.TotalPump_nadp = self.params['TotalPump_nadp_soma']

        neuron.h.k1_nadp = self.params['k1_nadp']
        neuron.h.k2_nadp = self.params['k2_nadp']
        neuron.h.k3_nadp = self.params['k3_nadp']
        neuron.h.DNa_nadp = self.params['DNa']

        neuron.h.nao0_na_ion = self.params['nao']
        den = (8.314e3 * (273.15 + neuron.h.celsius)) / 9.6485e4
        neuron.h.nai0_na_ion = neuron.h.nao0_na_ion * \
            np.exp(-self.params['ena'] / den)

    def remove_ca_dynamics(self):
        mechs = neuron.h.MechanismType(0)
        mechs.select("CaDynamics_E2")

        for sec in neuron.h.allsec():
            sec.push()
            mechs.remove()

    def insert_ca_mech(self):
        neuron.h.cai0_ca_ion = self.params['pump_ca_eq']
        neuron.h.cao0_ca_ion = self.params['cao']

        for sec in neuron.h.allsec():
            sec.insert('cadp')
            sec.insert('ncx')
            sec.TotalPump_cadp = self.params['TotalPump_cadp_dend']
            sec.imax_ncx = self.params['imax_ncx_dend']
        for sec in neuron.h.L5PC.axon:
            sec.TotalPump_cadp = self.params['TotalPump_cadp_axon']
            sec.imax_ncx = self.params['imax_ncx_axon']
        for sec in neuron.h.L5PC.soma:
            sec.TotalPump_cadp = self.params['TotalPump_cadp_soma']
            sec.imax_ncx = self.params['imax_ncx_soma']

        neuron.h.TotalEndBuffer_cadp = self.params['TotalEndBuffer']
        neuron.h.k2bufend_cadp = neuron.h.k1bufend_cadp * \
            self.params['EndBufferKd']
        neuron.h.kna_ncx = self.params['kna_ncx']
        neuron.h.kca_ncx = self.params['kca_ncx']
        neuron.h.gamma_ncx = self.params['gamma_ncx']
        neuron.h.ksat_ncx = self.params['ksat_ncx']

    def save_states(self, filename, new_mechs=True):
        vals = []
        for sec in neuron.h.allsec():
            for seg in sec:
                vals = vals + [seg.v, seg.Ih.m]
                if new_mechs:
                    vals = vals + [seg.nadp.pump,
                                   seg.nadp.pumpna,
                                   seg.nadp.na,
                                   seg.cadp.pump,
                                   seg.cadp.pumpca,
                                   seg.cadp.ca,
                                   seg.cadp.CaEndBuffer,
                                   seg.cadp.EndBuffer]

                if 'dend' not in sec.name():
                    vals = vals + [seg.Ca_HVA.h, seg.Ca_HVA.m,
                                   seg.Ca_LVAst.h, seg.Ca_LVAst.m]
                    vals = vals + [seg.Im.m, seg.SK_E2.z, seg.SKv3_1.m]
                    if 'axon' not in sec.name():
                        vals = vals + [seg.NaTs2_t.h, seg.NaTs2_t.m]

                    if not new_mechs:
                        vals = vals + [seg.cai]

                if 'axon' in sec.name():
                    vals = vals + [seg.K_Pst.h,
                                   seg.K_Pst.m,
                                   seg.K_Tst.h,
                                   seg.K_Tst.m,
                                   seg.Nap_Et2.h,
                                   seg.Nap_Et2.m,
                                   seg.NaTa_t.h,
                                   seg.NaTa_t.m]

        if new_mechs:
            vals = vals + [neuron.h.k4_nadp,
                           neuron.h.k2_cadp, neuron.h.k4_cadp]

        if filename is not None:
            f = open(filename, 'w')
            pickle.dump(vals, f)
            f.close()
        return vals

    def restore_states(self):

        vals = copy.deepcopy(self.rest_vals)
        for sec in neuron.h.allsec():

            for seg in sec:
                seg.v = vals.pop(0)
                seg.Ih.m = vals.pop(0)

                if self.converted:
                    if self.nadp:
                        seg.nadp.pump = vals.pop(0)
                        seg.nadp.pumpna = vals.pop(0)
                        seg.nadp.na = vals.pop(0)
                        seg.nai = seg.nadp.na
                        if ('dend' in sec.name() or 'apic' in sec.name()):
                         #                           seg.k4_coeff_nadp = self.k4_coeff_dend
                            seg.DNa_coeff_nadp = self.DNa_coeff_dend

                    else:
                        vals.pop(0)
                        vals.pop(0)
                        seg.nai = vals.pop(0)

                    seg.cadp.pump = vals.pop(0)
                    seg.cadp.pumpca = vals.pop(0)
                    seg.cadp.ca = vals.pop(0)
                    seg.cai = seg.cadp.ca
                    seg.cadp.CaEndBuffer = vals.pop(0)
                    seg.cadp.EndBuffer = vals.pop(0)

                if 'dend' not in sec.name():
                    seg.Ca_HVA.h = vals.pop(0)
                    seg.Ca_HVA.m = vals.pop(0)
                    seg.Ca_LVAst.h = vals.pop(0)
                    seg.Ca_LVAst.m = vals.pop(0)
                    seg.Im.m = vals.pop(0)
                    seg.SK_E2.z = vals.pop(0)
                    seg.SKv3_1.m = vals.pop(0)
                    if 'axon' not in sec.name():
                        seg.NaTs2_t.h = vals.pop(0)
                        seg.NaTs2_t.m = vals.pop(0)

                    if not self.converted:
                        seg.cai = vals.pop(0)
#                        if 'soma' in sec.name():
#                            print seg.cai

                if 'axon' in sec.name():
                    seg.K_Pst.h = vals.pop(0)
                    seg.K_Pst.m = vals.pop(0)

                    seg.K_Tst.h = vals.pop(0)
                    seg.K_Tst.m = vals.pop(0)

                    seg.Nap_Et2.h = vals.pop(0)
                    seg.Nap_Et2.m = vals.pop(0)

                    seg.NaTa_t.h = vals.pop(0)
                    seg.NaTa_t.m = vals.pop(0)

        if self.converted:
            if self.nadp:
                neuron.h.k4_nadp = vals.pop(0)
            else:
                vals.pop(0)
            neuron.h.k2_cadp = vals.pop(0)
            neuron.h.k4_cadp = vals.pop(0)

            neuron.h.fix_na_nadp = self.fixna
        if self.event_fun is not None:
            self.event_fun()

    def save_rest_state(self, init_run_duration, filename):
        self.run_model(init_run_duration, record=False)
        self.save_states(filename, self.converted)
