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

By default the content of the dx file is assumed to contain data in Volts (V) (or a 
potential expressed in PMEPot units, based on kT.e-1, which can be convered to Volts
via the --pmepot flag).

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
--xticks	[10]	: nb of ticks along the plot horizontal axis
--yticks	[7]	: nb of ticks along the plot vertical axis
--cticks	[10]	: nb of ticks on the colour bar
--pmepot		: use this flag to convert units from PMEPot to V

Volume to process
-----------------------------------------------------
--xmin		[0]	: position of lower delimiter on the x axis (as a %)
--ymin		[0]	: position of lower delimiter on the y axis (as a %)
--zmin		[0]	: position of lower delimiter on the z axis (as a %)
--xmax		[100]	: position of upper delimiter on the x axis (as a %)
--ymax		[100]	: position of upper delimiter on the y axis (as a %)
--zmax		[100]	: position of upper delimiter on the z axis (as a %)

Potential reference
-----------------------------------------------------
-r		[z]	: the potential will be set to 0 at the lower extremity of this axis (x,y or z)
--pad		[5]	: nb of slices used to calculate potential offset (set to 0 to not offset)

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
parser.add_argument('--xticks', nargs=1, dest='xticks', default=[10], type=int, help=argparse.SUPPRESS)
parser.add_argument('--yticks', nargs=1, dest='yticks', default=[7], type=int, help=argparse.SUPPRESS)
parser.add_argument('--cticks', nargs=1, dest='cticks', default=[10], type=int, help=argparse.SUPPRESS)
parser.add_argument('--pmepot', dest='pmepot', action='store_true', help=argparse.SUPPRESS)

#volume to process
parser.add_argument('--xmin', nargs=1, dest='xmin', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--ymin', nargs=1, dest='ymin', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--zmin', nargs=1, dest='zmin', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--xmax', nargs=1, dest='xmax', default=[100], type=float, help=argparse.SUPPRESS)
parser.add_argument('--ymax', nargs=1, dest='ymax', default=[100], type=float, help=argparse.SUPPRESS)
parser.add_argument('--zmax', nargs=1, dest='zmax', default=[100], type=float, help=argparse.SUPPRESS)

#potential reference
parser.add_argument('-r', dest='axisref', choices=['x','y','z'], default='z', help=argparse.SUPPRESS)
parser.add_argument('--pad', nargs=1, dest='pad', default=[5], type=int, help=argparse.SUPPRESS)

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
args.pad = args.pad[0]
args.xticks = args.xticks[0]
args.yticks = args.yticks[0]
args.cticks = args.cticks[0]

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

if args.pad < 0:
	print "Error: --pad cannot be negative. Set it to 0 to avoid offsetting the potential."
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
global nx_min
global ny_min
global nz_min
global nx_max
global ny_max
global nz_max

#=========================================================================================
# data loading
#=========================================================================================

def load_dx():

	global dims
	global data
	global coords_x
	global coords_y
	global coords_z
	global nx_min
	global ny_min
	global nz_min
	global nx_max
	global ny_max
	global nz_max

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
	
	#calculate delimiting indexes
	nx_min = int(args.xmin * dims[0] / float(100))
	ny_min = int(args.ymin * dims[1] / float(100))
	nz_min = int(args.zmin * dims[2] / float(100))
	nx_max = int(args.xmax * dims[0] / float(100))
	ny_max = int(args.ymax * dims[1] / float(100))
	nz_max = int(args.zmax * dims[2] / float(100))
	
	#center coordinates
	coords_x -= np.average(coords_x)
	coords_y -= np.average(coords_y)
	coords_z -= np.average(coords_z)
	coords_x = coords_x[nx_min:nx_max]
	coords_y = coords_y[ny_min:ny_max]
	coords_z = coords_z[nz_min:nz_max]
	
	return
	
#=========================================================================================
# averages
#=========================================================================================

def calc_profiles():
	
	global data_1D
	global data_2D

	
	# 1D average along chosen axis
	#-----------------------------
	if args.axis1D == "x":
		data_1D = np.zeros(nx_max-nx_min)
		for nx in range(nx_min,nx_max):
			data_1D[nx] = np.average(data[nx,ny_min:ny_max,nz_min:nz_max])
	elif args.axis1D == "y":
		data_1D = np.zeros(ny_max-ny_min)
		for ny in range(ny_min,ny_max):
			data_1D[ny] = np.average(data[nx_min:nx_max,ny,nz_min:nz_max])
	else:
		data_1D = np.zeros(nz_max-nz_min)
		for nz in range(nz_min,nz_max):
			data_1D[nz] = np.average(data[nx_min:nx_max,ny_min:ny_max,nz])

	# 2D average
	#-----------
	#case: xz
	if args.axis2D == "xz":
		data_2D = np.zeros((nx_max-nx_min,nz_max-nz_min))
		for nz in range(nz_min,nz_max):
			for nx in range(nx_min,nx_max):
				data_2D[nx,nz] = np.average(data[nx,ny_min:ny_max,nz])
	#case: yz
	elif args.axis2D == "yz":
		data_2D = np.zeros((ny_max-ny_min,nz_max-nz_min))
		for nz in range(nz_min,nz_max):
			for nx in range(ny_min,ny_max):
				data_2D[nx,nz] = np.average(data[nx_min:nx_max,ny,nz])
	#case: xy
	else:
		data_2D = np.zeros((nx_max-nx_min,ny_max-ny_min))
		for ny in range(ny_min,ny_max):
			for nx in range(ny_min,ny_max):
				data_2D[nx,nz] = np.average(data[nx,ny,nz_min:mz_max])
		
	#convert units to V
	#------------------
	if args.pmepot:
		#convert to kT to kJ (in PMEpot C++ code a temperature of 300 K is used to obtain kT)
		factor = 8.3144621 * 300 / float(1000) * 0.010364272
		data_1D *= factor
		data_2D *= factor

	#sets potential to 0 V at the lower extremity of the reference axis using specified padding
	#------------------------------------------------------------------------------------------
	if args.pad > 0:
		offset = np.average(data_1D[0:args.pad])
		data_1D -= offset
		data_2D -= offset

	#set upper an lower boundaries if need be
	#----------------------------------------
	if args.vmin == "auto":
		args.vmin = min(data_1D)
	else:
		args.vmin = float(args.vmin)
	if args.vmax == "auto":
		args.vmax = max(data_1D)
	else:
		args.vmax = float(args.vmax)

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
	output_xvg.write("# -> 1D axis: " + str(args.axis1D) + "\n")
	output_xvg.write("# -> 2D axis: " + str(args.axis2D) + "\n")
	output_xvg.write("# -> x axis: " + str(args.xmin) + "-" + str(args.xmax) + "\n")
	output_xvg.write("# -> y axis: " + str(args.ymin) + "-" + str(args.ymax) + "\n")
	output_xvg.write("# -> z axis: " + str(args.zmin) + "-" + str(args.zmax) + "\n")
	output_xvg.write("# -> pad: " + str(args.pad) + " slices on " + str(args.axisref) + " axis lower extremity\n")	
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label \"distance from box center along " + str(args.axis1D) + " (A)\"\n")
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
	if args.axis1D == "x":
		for r in range(0, len(data_1D)):
			results = str(round(coords_x[r],2)) + "	" + "{:.6e}".format(data_1D[r])
			output_xvg.write(results + "\n")		
	elif args.axis1D == "y":
		for r in range(0, len(data_1D)):
			results = str(round(coords_y[r],2)) + "	" + "{:.6e}".format(data_1D[r])
			output_xvg.write(results + "\n")		
	else:
		for r in range(0, len(data_1D)):
			results = str(round(coords_z[r],2)) + "	" + "{:.6e}".format(data_1D[r])
			output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

def graph_profile_1D():
	
	#filenames
	filename_svg = os.getcwd() + '/' + str(args.dxfilename[:-4]) + '_1D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Electrostatic profile along " + str(args.axis1D))
	
	#plot data
	ax = fig.add_subplot(111)
	if args.axis1D == "x":
		plt.plot(coords_x, data_1D, color = 'k', linewidth = 2)
		plt.hlines(0, min(coords_x), max(coords_x))
	elif args.axis1D == "y":
		plt.plot(coords_y, data_1D, color = 'k', linewidth = 2)
		plt.hlines(0, min(coords_y), max(coords_y))
	else:
		plt.plot(coords_z, data_1D, color = 'k', linewidth = 2)
		plt.hlines(0, min(coords_z), max(coords_z))	
	plt.vlines(-21, args.vmin, args.vmax, linestyles = 'dashed')
	plt.vlines(21, args.vmin, args.vmax, linestyles = 'dashed')
	plt.vlines(0, args.vmin, args.vmax, linestyles = 'dashdot')
	plt.xlabel(str(args.axis1D) + ' distance to box center ($\AA$)')
	plt.ylabel('electrostatic potential (V)')
	
	#save figure
	ax.set_ylim(args.vmin, args.vmax)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=args.xticks))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=args.yticks))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_svg)
	plt.close()

	return

def graph_profile_2D():	

	#filenames
	filename_svg = os.getcwd() + '/' + str(args.dxfilename[:-4]) + '_2D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Electrostatic profile slice")

	#plot data
	ax = fig.add_subplot(111)
	if args.axis2D == "xz":
		data_2D_oriented = np.zeros((np.shape(data_2D)[1],np.shape(data_2D)[0]))
		for nx in range(nx_min, nx_max):
			for nz in range(nz_min, nz_max):
				data_2D_oriented[nz,nx] = data_2D[nx,nz_max-1-nz]
		im = plt.imshow(data_2D_oriented, extent = [min(coords_x),max(coords_x),min(coords_z),max(coords_z)], vmin = args.vmin, vmax = args.vmax, cmap = matplotlib.cm.jet_r)
		ax.set_xlim(min(coords_x), max(coords_x))
		ax.set_ylim(min(coords_z), max(coords_z))
		plt.xlabel('x axis ($\AA$)')
		plt.ylabel('z axis ($\AA$)')
	elif args.axis2D == "yz":
		data_2D_oriented = np.zeros((np.shape(data_2D)[1],np.shape(data_2D)[0]))
		for ny in range(ny_min, ny_max):
			for nz in range(nz_min, nz_max):
				data_2D_oriented[nz,ny] = data_2D[ny,nz_max-1-nz]
		im = plt.imshow(data_2D_oriented, extent = [min(coords_y),max(coords_y),min(coords_z),max(coords_z)], vmin = args.vmin, vmax = args.vmax, cmap = matplotlib.cm.jet_r)
		ax.set_xlim(min(coords_y), max(coords_y))
		ax.set_ylim(min(coords_z), max(coords_z))
		plt.xlabel('y axis ($\AA$)')
		plt.ylabel('z axis ($\AA$)')
	else:
		data_2D_oriented = np.zeros((np.shape(data_2D)[1],np.shape(data_2D)[0]))
		for nx in range(nx_min, nx_max):
			for ny in range(ny_min, ny_max):
				data_2D_oriented[ny,nx] = data_2D[nx,ny_max-1-ny]
		im = plt.imshow(data_2D_oriented, extent = [min(coords_x),max(coords_x),min(coords_y),max(coords_y)], vmin = args.vmin, vmax = args.vmax, cmap = matplotlib.cm.jet_r)
		ax.set_xlim(min(coords_x), max(coords_x))
		ax.set_ylim(min(coords_y), max(coords_y))
		plt.xlabel('x axis ($\AA$)')
		plt.ylabel('y axis ($\AA$)')

	if args.axis2D != "xy":
		plt.vlines(-21, args.vmin, args.vmax, linestyles = 'dashed')
		plt.vlines(21, args.vmin, args.vmax, linestyles = 'dashed')
		plt.vlines(0, args.vmin, args.vmax, linestyles = 'dashdot')
	
	#color bar
	cax = fig.add_axes([0.83, 0.2, 0.025, 0.65])
	cbar = fig.colorbar(im, orientation='vertical', cax=cax)
	cbar.ax.tick_params(axis='y', direction='out')
	cbar.set_label(r'potential (V)')
	plt.setp(cbar.ax.yaxis.get_majorticklabels(), fontsize = "small")
	cbar.locator = MaxNLocator(nbins=args.cticks)
	cbar.update_ticks()
	cbar.ax.yaxis.labelpad = 10

	#save figure	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=args.xticks))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=args.yticks))
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
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
graph_profile_1D()
graph_profile_2D()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully!"
print ""
sys.exit(0)
