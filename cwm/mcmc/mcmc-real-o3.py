import os
import numpy as np
import pyfstat
from pyfstat.utils import get_predict_fstat_parameters_from_dict
label = "G347Real398Hzv6-03"
outdir = os.path.join("mcmc", label)
logger = pyfstat.set_up_logger(label=label, outdir=outdir)

# Properties of the GW data
data_parameters = {
"tstart": 1238166483,
"duration": 360 * 86400,
"detectors": "H1,L1",
}

alpha = 4.5093705
delta = -0.6951891
tend = data_parameters["tstart"] + data_parameters["duration"]
mid_time = 0.5 * (data_parameters["tstart"] + tend)
signal_parameters = {
"Alpha": alpha,
"Delta": delta,
"tref": mid_time
}


f0=398.02177518
df0=8e-09
f1=-7.40622425e-09
df1=8.1140004e-14
f2=5.64491762e-19
df2=2.93091845e-20


"""
f0=284.03648756
df0=7.26464104e-09
f1=-3.71331841e-10
df1=1.1140004e-14
f2=5.76806436e-20
df2=1.26172298e-20
f3=4.73914781e-27
df3=7.51342517e-33
f4=-7.60072233e-34
df4=2.57921688e-43
"""
"""
f0=46.09999612
df0=7.23953934e-09
f1=-5.17427737e-10
df1=1.04663628e-14
f2=4.79583446e-21
df2=1.2487763e-20

"""

df0 = 1e-6
df1 = 1e-13
df2 = 1e-19
ds = 1e-4
theta_prior = {
"Alpha": {
"type": "unif",
"lower": signal_parameters["Alpha"] - ds/2.0,
"upper": signal_parameters["Alpha"] + ds/2.0,
},
"Delta": {
"type": "unif",
"lower": signal_parameters["Delta"] - ds,
"upper": signal_parameters["Delta"] + ds,
},
"F0": {
"type": "unif",
"lower": f0 - df0,
"upper": f0 + df0,
},
"F1": {
"type": "unif",
"lower": f1 - df1 * 1.0,
"upper": f1 + df1 * 1.0,
},
"F2": {
"type": "unif",
"lower": f2 - df2 * 1.0,
"upper": f2 + df2 * 1.0,
},
#"F3": {
#"type": "unif",
#"lower": f3 - df3 * 2.0,
#"upper": f3 + df3 * 2.0,
#},
#"F4": {
#"type": "unif",
#"lower": f4 - df4 * 2.0,
#"upper": f4 + df4 * 2.0,
#}
}
"""
theta_prior = {
"Alpha": {
"type": "unif",
"lower": signal_parameters["Alpha"] - 1e-5,
"upper": signal_parameters["Alpha"] + 1e-5,
},
"Delta": {
"type": "unif",
"lower": signal_parameters["Delta"] - 1e-5,
"upper": signal_parameters["Delta"] + 1e-5,
},
"F0": {
"type": "unif",
"lower": f0 - df0 * 1.0,
"upper": f0 + df0 * 1.0,
},
"F1": {
"type": "unif",
"lower": f1 - df1 * 1.0,
"upper": f1 + df1 * 1.0,
},
"F2": {
"type": "unif",
"lower": f2 - df2 * 1.0,
"upper": f2 + df2 * 1.0,
},
#"F3": {
#"type": "unif",
#"lower": f3 - df3 * 2.0,
#"upper": f3 + df3 * 2.0,
#},
#"F4": {
#"type": "unif",
#"lower": f4 - df4 * 2.0,
#"upper": f4 + df4 * 2.0,
#}
}
"""
ntemps = 2
log10beta_min = -0.5
nwalkers = 200
nsteps = [300, 300]
mcmc = pyfstat.MCMCSearch(
label=label,
outdir=outdir,
sftfilepattern='/home/hoitim.cheung/snrsearch/test_osg/SFTs/NarrowBand/360days/H1/398/H-6772_H1_1800SFT_NBF0396Hz0W0005Hz0_O3_Gated_Sub60Hz-1238166483-20539351.sft;/home/hoitim.cheung/snrsearch/test_osg/SFTs/NarrowBand/360days/L1/398/L-10129_L1_1800SFT_NBF0396Hz0W0005Hz0_O3_Gated_Sub60Hz-1238184901-31172922.sft',
#sftfilepattern='/home/hoitim.cheung/snrsearch/o4/SFTs/narrowBand/240days/H1/284/H-7338_H1_1800SFT_O4RUN+R1+CGDSCALIBSTRAINCLEANGATEDG01+WTKEY5_NBF0283Hz0W0003Hz0-1368980712-20472783.sft;/home/hoitim.cheung/snrsearch/o4/SFTs/narrowBand/240days/L1/284/L-7622_L1_1800SFT_O4RUN+R1+CGDSCALIBSTRAINCLEANGATEDG01+WTKEY5_NBF0283Hz0W0003Hz0-1368976350-20451259.sft',
#sftfilepattern='/home/hoitim.cheung/snrsearch/o4/SFTs/narrowBand/240days/H1/46/H-7338_H1_1800SFT_O4RUN+R1+CGDSCALIBSTRAINCLEANGATEDG01+WTKEY5_NBF0045Hz0W0003Hz0-1368980712-20472783.sft;/home/hoitim.cheung/snrsearch/o4/SFTs/narrowBand/240days/L1/46/L-7622_L1_1800SFT_O4RUN+R1+CGDSCALIBSTRAINCLEANGATEDG01+WTKEY5_NBF0045Hz0W0003Hz0-1368976350-20451259.sft',  
theta_prior=theta_prior,
tref=mid_time,
minStartTime=data_parameters["tstart"],
maxStartTime=tend,
nsteps=nsteps,
nwalkers=nwalkers,
ntemps=ntemps,
log10beta_min=log10beta_min,
)
#mcmc.transform_dictionary = dict(
#F0=dict(subtractor=f0, symbol="$f-f^\\mathrm{s}$"),
#F1=dict(subtractor=f1, symbol="$f1-f1^\\mathrm{s}$"),
#F2=dict(subtractor=f2, symbol="$f2-f2^\\mathrm{s}$"),
#F3=dict(subtractor=f3, symbol="$f3-f3^\\mathrm{s}$"),
#F4=dict(subtractor=f4, symbol="$f4-f4^\\mathrm{s}$"),
#)
mcmc.transform_dictionary = dict(
F0=dict(subtractor=0, symbol="$f$"),
F1=dict(subtractor=0, symbol="$f1$"),
F2=dict(subtractor=0, symbol="$f2$"),
Alpha=dict(subtractor=0, symbol="$alpha$"),
Delta=dict(subtractor=0, symbol="$delta$"),

#F3=dict(subtractor=f3, symbol="$f3$"),
#F4=dict(subtractor=f4, symbol="$f4-f4^\\mathrm{s}$"),
)
mcmc.run(
walker_plot_args={"plot_det_stat": True}
)
mcmc.print_summary()
mcmc.plot_corner(add_prior=True)
mcmc.plot_prior_posterior()

mcmc.generate_loudest()

# plot cumulative 2F, first building a dict as required for PredictFStat
#d, maxtwoF = mcmc.get_max_twoF()
#for key, val in mcmc.theta_prior.items():
#    if key not in d:
#        d[key] = val
#d["h0"] = data.h0
#d["cosi"] = data.cosi
#d["psi"] = data.psi
#PFS_input = get_predict_fstat_parameters_from_dict(d)
#mcmc.plot_cumulative_max(PFS_input=PFS_input)
