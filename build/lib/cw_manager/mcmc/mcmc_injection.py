#!/home/hoitim.cheung/.conda/envs/cw/bin/python
import os
import numpy as np
import pyfstat
from pyfstat.utils import get_predict_fstat_parameters_from_dict
parser = argparse.ArgumentParser(description='MCMC for outlier follow-up.')
parser.add_argument('--target', type=str, help='Target to be analyzed', default='CassA')
parser.add_argument('--obsDay',type=float,help='Total observation time',default=240)  
parser.add_argument('--freq',type=float,help='Nominal frequency',default=20)
parser.add_argument('--f1dot',type=float,help='First order of frequency evolution',default=-1e-9)
parser.add_argument('--f2dot',type=float,help='Second order of frequency evolution',default=1e-18)
parser.add_argument('--f1dot',type=float,help='First order of frequency evolution',default=1e-9)
parser.add_argument('--f2dot',type=float,help='Second order of frequency evolution',default=1e-18)
#parser.add_argument('--f3dot',type=float,help='Third order of frequency evolution',default=1e-27)
#parser.add_argument('--f4dot',type=float,help='Forth order of frequency evolution',default=-1e-34)

parser.add_argument('--df',type=float,help='Search range for  nominal frequency',default=1e-6)
parser.add_argument('--df1dot',type=float,help='Search range for first order of frequency evolution',default=1e-13)
parser.add_argument('--df2dot',type=float,help='Search range second order of frequency evolution',default=1e-19)
parser.add_argument('--startTime',type=int,help='GPS time for the start of the observation run',default=1368970000)
parser.add_argument('--ra',type=float,help='Right ascension of the target.',default=6.1237704239609)
parser.add_argument('--dec',type=float,help='Declination of the target.',default=1.0264578036951)
parser.add_argument('--dra',type=float,help='Search range for right ascension of the target.',default=1e-5)
parser.add_argument('--ddec',type=float,help='Search range for declination of the target.',default=1e-5)
parser.add_argument('--hnoise',type=float,help='',default=1e-23)
parser.add_argument('--h0',type=float,help='',default=1e-25)

args = parser.parse_args()
target = args.target
obsDay = args.obsDay
freq = args.freq
f1dot = args.f1dot
f2dot = args.f2dot
df = args.df
df1dot = args.df1dot
df2dot = args.df2dot
startTime = args.startTime
ra = args.ra
dec = args.dec
dra = args.dra
ddec = args.ddec
hnoise = args.hnoise
h0 = args.h0


label = "{}Mock{}Hz".format(target, freq)
outdir = os.path.join("mcmc", label)
logger = pyfstat.set_up_logger(label=label, outdir=outdir)

# Properties of the GW data
data_parameters = {
"tstart": startTime,
"sqrtSX": hnoise,
"duration": obsDay * 86400,
"detectors": "H1,L1",
}
tend = data_parameters["tstart"] + data_parameters["duration"]
mid_time = 0.5 * (data_parameters["tstart"] + tend)



signal_parameters = {
"F0": freq,
"F1": f1dot,
"F2": f2dot,
"Alpha": ra,
"Delta":  dec,
"tref": mid_time,
"h0": h0,
"cosi": 1.0,
}
data = pyfstat.Writer(
label=label, outdir=outdir, **data_parameters, **signal_parameters
)

data.make_data()
theta_prior = {
"Alpha": {
"type": "unif",
"lower": ra - dra,
"upper": ra + dra,
},
"Delta": {
"type": "unif",
"lower": dec - ddec,
"upper": dec + ddec,
},
"F0": {
"type": "unif",
"lower": freq - df,
"upper": freq + df,
},
"F1": {
"type": "unif",
"lower": f1dot - df1dot,
"upper": f1dot + df1dot,
},
"F2": {
"type": "unif",
"lower": f2dot - df2dot,
"upper": f2dot + df2dot,
},
#"F3": {
#"type": "unif",
#"lower": f3dot - df3dot,
#"upper": f3dot + df3dot,
#},
#"F4": {
#"type": "unif",
#"lower": f4dot - df4dot,
#"upper": f4dot + df4dot,
#}
}

ntemps = 2
log10beta_min = -0.5
nwalkers = 200
nsteps = [300, 300]
mcmc = pyfstat.MCMCSearch(
label=label,
outdir=outdir,
sftfilepattern=data.sftfilepath,
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
