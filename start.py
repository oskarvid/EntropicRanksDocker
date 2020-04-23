#!/usr/bin/env python2.7

import yaml
import argparse
import docker

# This argparse thing is not in use, I will keep it for now
parser = argparse.ArgumentParser(description='Define job parameters for query submission')
parser.add_argument('--data', type = str, action = 'store', help = 'GSE data table')
parser.add_argument('--vector', type = str, action = 'store', help = 'Population vector')
parser.add_argument('--granularity', type = str, action = 'store', help = 'Data granularity (?)')
parser.add_argument('--origin', type = str, action = 'store', help = 'Data origin')
parser.add_argument('--supervised', type = str, action = 'store', help = 'Supervised mode, valid options are TRUE or FALSE')
parser.add_argument('--log', type = str, action = 'store', help = 'Set to TRUE to generate log file')
parser.add_argument('--plots', type = str, action = 'store', help = 'Set to TRUE to generate plots')
parser.add_argument('--outputs', type = str, action = 'store', help = 'Set to TRUE to save output files')
parser.add_argument('--islogged', type = str, action = 'store', help = 'Set to TRUE if values are logarithmic')
parser.add_argument('--logbase', type = str, action = 'store', help = 'Integer that sets log base')
parser.add_argument('--hugefeaturelist', type = str, action = 'store', help = 'Set to true to generate an extensive feature list')
args = parser.parse_args()

# Open the config file
with open("config.yml") as f:
        configYml = yaml.load(f, Loader=yaml.FullLoader)

# Assign variables from config file
entropicranksdir = configYml[0]["entropicranksdir"]
datatable = configYml[0]["datatable"]
vector = configYml[0]["vector"]
origin = configYml[0]["origin"]
granularity = configYml[0]["granularity"]
supervised = configYml[0]["supervised"]
logfile = configYml[0]["logfile"]
plots = configYml[0]["plots"]
outputs = configYml[0]["outputs"]
islogged = configYml[0]["islogged"]
logbase = configYml[0]["logbase"]
hugefeaturelist = configYml[0]["hugefeaturelist"]

# entropic ranks arguments
arguments =  "/data/" + datatable + " /data/" + vector + " " + origin +  " " + granularity + " " + supervised + " " + logfile + " " + plots + " " + outputs + " " + islogged + " " + logbase + " " + hugefeaturelist

client = docker.from_env()
# Run docker and r-script
container = client.containers.run(
	'oskarv/entropicrank',
	'Rscript /data/EntropicRanks.R ' + arguments,
	volumes = {
		entropicranksdir: {
			'bind':'/data',
			'mode':'rw',
		},
	}, detach = True
)
for line in container.logs(stream=True):
	print line.strip()
