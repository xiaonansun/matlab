
t = readNPY('Z:\simon_ephys\Optotest1\_spikeglx_ephysData_probe00_g1\spike_times.npy');
c = readNPY('Z:\simon_ephys\Optotest1\_spikeglx_ephysData_probe00_g1\spike_clusters.npy');
sample_rate = 30000;

spikes = t(c==1);
