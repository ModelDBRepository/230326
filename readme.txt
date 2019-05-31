Model files from the paper:

Zylbertal et al., "The Slow Dynamics of Intracellular Sodium Concentration Increase the Time Window of Neuronal Integration: A Simulation Study", XXXXXXX (2017)
The file synaptic_train.py reproduces the protocol used in Fig. 5E of the article
by calling the module pyramidal.py.
The model is adapted from Hay, E., Hill, S., Schürmann, F., Markram, H., and Segev, I. (2011). Models of Neocortical Layer 5b Pyramidal
Cells Capturing a Wide Range of Dendritic and Perisomatic Active Properties. PLoS Comput. Biol. 7, e1002107. doi:10.1371/journal.pcbi.1002107.

Questions on how to use this model should be directed to
asaph.zylbertal@mail.huji.ac.il

Synopsis:

Changes in intracellular Na+ concentration ([Na+]i) are rarely taken into account when neuronal activity is examined. As opposed to Ca2+, [Na+]i dynamics are strongly affected by longitudinal diffusion, and therefore they are governed by the morphological structure of the neurons, in addition to the localization of influx and efflux mechanisms. Here we examined [Na+]i dynamics and their effects on neuronal computation in three multi-compartmental neuronal models, representing three distinct cell types: accessory olfactory bulb (AOB) mitral cells, cortical layer V pyramidal cells, and cerebellar Purkinje cells. We added [Na+]i as a state variable to these models, and allowed it to modulate the Na+ Nernst potential, the Na+-K+ pump current, and the Na+-Ca2+ exchanger rate. Our results indicate that in most cases [Na+]i dynamics are significantly slower than [Ca2+]i dynamics, and thus may exert a prolonged influence on neuronal computation in a neuronal type specific manner. We show that [Na+]i dynamics affect neuronal activity via three main processes: reduction of EPSP amplitude in repeatedly active synapses due to reduction of the Na+ Nernst potential; activity-dependent hyperpolarization due to increased activity of the Na+-K+ pump; specific tagging of active synapses by extended Ca2+ elevation, intensified by concurrent back-propagating action potentials or complex spikes. Thus, we conclude that [Na+]i dynamics should be considered whenever synaptic plasticity, extensive synaptic input, or bursting activity are examined.

The example protocol simulates the modified pyramidal cell with a train of synaptic inputs in the dendritic hot zone, concurent with evoked calcium spikes (see article figure 5E).

Example use:

Extract the archive, run nrnivmodl in the channels directory
(linux/unix) or mknrndll (mswin or mac os x) (see
http://senselab.med.yale.edu/ModelDB/NEURON_DwnldGuide.html for more
help) to compile the channels, and run the file synaptic_train.py After a
while, it will plot the membrane potential (top), the dendritic sodium concentration (middle) and the dendritic calcium concentration (bottom).
