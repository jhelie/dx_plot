#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'dx_plot', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/dx_plot
**********************************************

[ DESCRIPTION ]
 
This loads a data file containing electrostatic potential values stored in the OpenDX
format and plots it.


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - numpy
 - scipy


[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: dx file
-a		[z]	: axis along which to produce the 1D graph (x,y or z)
-s		[xz]	: axes of slice for 2D graphs (xz,yz or xy)
--vmax			: upper limit of scale
--vmin			: lower limit of scale

Volume to process
-----------------------------------------------------
--xmin		[0]	: position of lower delimiter on the x axis (as a %)
--ymin		[0]	: position of lower delimiter on the y axis (as a %)
--zmin		[0]	: position of lower delimiter on the z axis (as a %)
--xmax		[100]	: position of upper delimiter on the x axis (as a %)
--ymax		[100]	: position of upper delimiter on the y axis (as a %)
--zmax		[100]	: position of upper delimiter on the z axis (as a %)

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs=1, dest='dxfilename', help=argparse.SUPPRESS, required=True)
parser.add_argument('-a', dest='axis1D', choices=['x','y','z'], default='z', help=argparse.SUPPRESS)
parser.add_argument('-s', dest='axis2D', choices=['xz','yz','xy'], default='xz', help=argparse.SUPPRESS)
parser.add_argument('--vmax', nargs=1, dest='vmax', default=['auto'], help=argparse.SUPPRESS)
parser.add_argument('--vmin', nargs=1, dest='vmin', default=['auto'], help=argparse.SUPPRESS)

#volume to process
parser.add_argument('--xmin', dest='xmin', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--ymin', dest='ymin', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--zmin', dest='zmin', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--xmax', dest='xmax', default=[100], type=float, help=argparse.SUPPRESS)
parser.add_argument('--ymax', dest='ymax', default=[100], type=float, help=argparse.SUPPRESS)
parser.add_argument('--zmax', dest='zmax', default=[100], type=float, help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.dxfilename = args.dxfilename[0]
args.vmax = args.vmax[0]
args.vmin = args.vmin[0]
args.xmin = args.xmin[0]
args.ymin = args.ymin[0]
args.zmin = args.zmin[0]
args.xmax = args.xmax[0]
args.ymax = args.ymax[0]
args.zmax = args.zmax[0]

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.density
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if not os.path.isfile(args.dxfilename):
	print "Error: file " + str(args.dxfilename) + " not found."
	sys.exit(1)

if args.vmin > args.vmax:
	print "Error: --vmin must be smaller than --vmax."
	sys.exit(1)

if args.xmin < 0:
	print "Error: --xmin must be > 0."
	sys.exit(1)
if args.ymin < 0:
	print "Error: --ymin must be > 0."
	sys.exit(1)
if args.zmin < 0:
	print "Error: --zmin must be > 0."
	sys.exit(1)

if args.xmax > 100:
	print "Error: --xmax must be < 100."
	sys.exit(1)
if args.ymax > 100:
	print "Error: --ymax must be < 100."
	sys.exit(1)
if args.zmax > 100:
	print "Error: --zmax must be < 100."
	sys.exit(1)

if args.xmin > args.xmax:
	print "Error: --xmin must be smaller than --xmax."
	sys.exit(1)
if args.ymin > args.ymax:
	print "Error: --ymin must be smaller than --ymax."
	sys.exit(1)
if args.zmin > args.zmax:
	print "Error: --zmin must be smaller than --zmax."
	sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

global dims
global data
global data_1D
global data_2D
global coords_x
global coords_y
global coords_z

#=========================================================================================
# data loading
#=========================================================================================

def load_dx():

	global dims
	global data
	global coords_x
	global coords_y
	global coords_z

	#potential value in each voxel
	g = MDAnalysis.analysis.density.Grid()
	g.load(str(args.dxfilename))	
	data = g.grid
	
	#coords of bins along each dimensions
	dims = np.shape(g.grid)
	coords_x = np.zeros(dims[0])
	coords_y = np.zeros(dims[1])
	coords_z = np.zeros(dims[2])
	tmp_x = g.edges[0]
	tmp_y = g.edges[1]
	tmp_z = g.edges[2]
	for nx in range(0,dims[0]):
		coords_x[nx] = (tmp_x[nx+1] + tmp_x[nx])/float(2)
	for ny in range(0,dims[1]):
		coords_y[ny] = (tmp_y[ny+1] + tmp_y[ny])/float(2)
	for nz in range(0,dims[2]):
		coords_z[nz] = (tmp_z[nz+1] + tmp_z[nz])/float(2)

	
	#center z coordinate around box center
	coords_z -= np.average(coords_z)
	
	return
	
#=========================================================================================
# averages
#=========================================================================================

def calc_profiles():
	
	global data_1D
	global data_2D
	
	# 1D average: along z
	#--------------------
	data_1D = np.zeros(dims[2])
	for nz in range(0,dims[2]):
		data_1D[nz] = np.average(data[:,:,nz])

	# 2D average: z versus x
	#-----------------------
	data_2D = np.zeros((dims[0],dims[2]))
	for nz in range(0,dims[2]):
		for nx in range(0,dims[0]):
			data_2D[nx,nz] = np.average(data[nx,:,nz])
	
	#convert units to V
	#------------------
	#convert to kT to kJ (in PMEpot C++ code a temperature of 300 K is used to obtain kT)
	factor = 8.3144621 * 300 / float(1000) * 0.010364272
	data_1D *= factor
	data_2D *= factor

	#sets potential to 0 V in solvent (using left / extracellular side of the membrane)
	#--------------------------------
	offset = np.average(data_1D[10:15])
	data_1D -= offset
	data_2D -= offset

	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.dxfilename[:-4]) + '_1D.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [1D average content of " + str(args.dxfilename) + " - written by dx_plot v" + str(version_nb) + "]\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label \"z distance from bilayer center (A)\"\n")
	output_xvg.write("@ yaxis label \"potential (V)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 0\n")
	output_xvg.write("@ s0 legend \"potential\"\n")
	
	#data
	for r in range(0, len(data_1D)):
		results = str(round(coords_z[r],2)) + "	" + "{:.6e}".format(data_1D[r])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

def graph_profiles():

	#1D profile
	#----------
	
	#filenames
	filename_svg = os.getcwd() + '/' + str(args.dxfilename[:-4]) + '_1D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Electrostatic profile along z")

	#plot data
	ax = fig.add_subplot(111)
	plt.plot(coords_z, data_1D, color = 'k', linewidth = 2)
	plt.vlines(-21, min(data_1D), max(data_1D), linestyles = 'dashed')
	plt.vlines(21, min(data_1D), max(data_1D), linestyles = 'dashed')
	plt.vlines(0, min(data_1D), max(data_1D), linestyles = 'dashdot')
	plt.hlines(0, min(coords_z), max(coords_z))
	plt.xlabel('z distance to bilayer center ($\AA$)')
	plt.ylabel('electrostatic potential (V)')
	
	#save figure
	ax.set_xlim(min(coords_z), max(coords_z))
	#ax.set_ylim(min_density_charges, max_density_charges)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_svg)
	plt.close()


	#2D profile
	#----------
	
	#filenames
	filename_svg = os.getcwd() + '/' + str(args.dxfilename[:-4]) + '_2D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Electrostatic profile slice")

	#plot data
	ax = fig.add_subplot(111)
	im = plt.imshow(data_2D)
	plt.vlines(-21, min(data_2D[:,0]), max(data_2D[:,0]), linestyles = 'dashed')
	plt.vlines(21, min(data_2D[:,0]), max(data_2D[:,0]), linestyles = 'dashed')
	plt.vlines(0, min(data_2D[:,0]), max(data_2D[:,0]), linestyles = 'dashdot')
	plt.xlabel('z distance to bilayer center ($\AA$)')
	plt.ylabel('x axis ($\AA$)')
	
	#color bar
	cax = fig.add_axes([0.85, 0.26, 0.025, 0.48])
	cbar = fig.colorbar(im, orientation='vertical', cax=cax)
	cbar.ax.tick_params(axis='y', direction='out')
	cbar.set_label(r'potential (V)')
		
	#save figure
	ax.set_xlim(0, dims[2])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.1, right = 0.8)
	fig.savefig(filename_svg)
	plt.close()


	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading file..."
load_dx()

print "\nCalculating profiles..."
calc_profiles()
write_xvg()
graph_profiles()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully!"
print ""
sys.exit(0)
