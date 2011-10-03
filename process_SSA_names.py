#!/usr/bin/env python
""" Parsing and analyzing the Social Security Administration's baby name
    database, found at http://www.ssa.gov/oact/babynames/limits.html"""

__author__ = "Vedant Misra"
__copyright__ = "Copyright 2011, Vedant Misra"
__license__ = "CC-BY-SA"
__version__ = "0.1"
__email__ = "vedant.misra@gmail.com"
__status__ = "Demo"

import matplotlib.pyplot as plt
import simplejson as json
import numpy as np
import matplotlib
import sys
import os

localPath = os.path.dirname(__file__)
natlDataPath = os.path.join(localPath, 'data', 'national')

def buildDataset(natlData):
    names = {}
    for file in natlData:
	print "Reading", file, "..."
	year = int(file[3:7])
	f = open(os.path.join(natlDataPath, file), 'r')
	for row in f.read().split("\n"):
	    row = row.strip()
	    try:
		name, gender, count = row.split(",")
	    except ValueError:
		# found blank line
		continue
	    nameGender = name + "_" + gender
	    count = int(count)
	    if nameGender in names:
		names[nameGender][year] = count
	    else:
		names[nameGender] = {year: count}
    return names

def buildMatrix(data):
    names = data.keys()
    years = range(1800, 2011)
    matrix = np.zeros(shape=(len(names), len(years)))
    for i in range(len(names)):
	if not i % 1000:
	    print "Building row", i, "..."
	nameData = data[names[i]]
	for j in range(len(years)):
	    try:
		matrix[i][j] = nameData[str(years[j])]
	    except KeyError:
		matrix[i][j] = 0
    return matrix
	

def buildJson():
    natlDataFiles = filter(
	lambda x: os.path.splitext(x)[1] == ".txt",
	os.listdir(natlDataPath))
    print "Writing data to names.json..."
    data = buildDataset(sorted(natlDataFiles))
    json.dump(data, open("names.json", 'w'))
    return data

def PCA(data):
    data = data[:2000] 
    numRows, numDims = data.shape
    # mean centering
    dataMean = data.mean(axis=0)
    for i in range(numRows):
	data[i] -= dataMean
    eigenVals, eigenVecs = np.linalg.eigh(np.dot(data, data.T))
    temp = np.dot(data.T, eigenVecs).T
    V = temp[::-1]
    S = np.sqrt(eigenVals)[::-1]
    return V, S, dataMean

def getValues(nameData, normed=None, padded=False):
    values = []
    if not padded:
	years = []
	for key in sorted(nameData.keys()):
	    years.append(int(key))
	    values.append(nameData[key])
    if padded:
	years = range(1800, 2011)
	for year in years:
	    try:
		values.append(nameData[str(year)])
	    except KeyError:
		values.append(0)
    if normed == 'sum':
	sumValues = float(sum(values))
	values = map(lambda x: x / sumValues, values)
    elif normed == 'max':
	maxValue = float(max(values))
	values = map(lambda x: x / maxValue, values)
    return years, values

def getMostPopular(data):
    mostPopular = []
    i = 0
    for name in sorted(data.keys()):
	i += 1
	if not i % 1000:
	    sys.stderr.write(name + "\n")
	appearances = sum([data[name][key] for key in data[name].keys()])
	mostPopular.append([name, appearances])
    mostPopular.sort(lambda x, y: -cmp(x[1], y[1]))
    return mostPopular

def plotNames(data, names, smooth=0):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_alpha = 0
    for i in range(len(names)):
        name =  names[i]
	try:
	    years, values = getValues(data[name], normed='max', padded=True)
	except KeyError:
	    continue
        if smooth != 0:
            window = np.ones(smooth, 'd')
            values = np.convolve(window/window.sum(),
                values, mode='valid')
            years = years[1:-1]
	#ax.plot(years[-30:], values[-30:], label=name, color=str(plot_color))
        if i < 50:
            plot_color = 'k'
            plot_alpha += (.8/50)
        else:
            plot_color = 'r'
            plot_alpha -= (.8/50)
        if i == 50:
            plot_alpha = .8
	ax.plot(years, values, plot_color, label=name, alpha=plot_alpha)
    #plt.legend()
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(.1))
    #ax.xaxis.grid(True, which='minor')
    #ax.yaxis.grid(True, which='major')
    ax.set_xlim([1870, 2010])
    ax.set_xlabel('Time')
    ax.set_ylabel('Correlation')
    matplotlib.rcParams['font.sans-serif'] = 'Helvetica Neue'
    plt.title('Highest and lowest correlates with the male given name "Brody"')
    plt.savefig("ssa-names_highest-and-lowest-correlates_brody-M.png", dpi=75)

def computeCorrelations(data, index_name):
    index_years, index_values = getValues(data[index_name], padded=True)
    correlations = []
    i = 0
    for name in sorted(data.keys()):
        i += 1
        if not i % 1000: print name
        years, values = getValues(data[name], padded=True)
        correlations.append((name, np.corrcoef(values, index_values)[0][1]))
    correlations.sort(lambda x, y: -cmp(x[1], y[1]))
    return correlations

def computePCA(datafile="names.npy"):
    print "Performing PCA..."
    V, S, dataMean = PCA(np.load(open(datafile, 'r')))
    matrix = np.load(open("names.npy", 'r'))
    print "Plotting..."
    plt.plot(V[0,:], V[1,:], "ob")
    plt.axis('equal')
    plt.savefig('PCA.pdf')

if __name__ == "__main__":
    # Generate "names.json" file that contains national data for all years
    buildJson()
    
    # Build a numpy matrix with names as rows, years as columns
    matrix_data = buildMatrix(json.load(open("names.json", 'r')))
    np.save(open("names.npy", 'w'), matrix_data)

    # Perform Principal Component Analysis with the input matrix above
    computePCA("names.npy")
   
    # Generate plots for certain names
    data = json.load(open("names.json", 'r'))
    plotNames(data, ['James_M', 'Jill_F'])
    
    # To compute and serialize correlations. Plots the top and bottom 50.
    data = json.load(open("names.json", 'r'))
    correlations = computeCorrelations(data, 'Brody_M')
    correlations = json.load(open('correlations.json', 'r'))
    top_correlations = correlations[:50]
    lowest_correlations = correlations[-50:]
    for c in top_correlations + lowest_correlations:
        name = c[0].split("_")
        print name[0] + " (" + name[1] + ")" + "\t" + str("%.3f" % c[1])
    json.dump(correlations, open('correlations.json', 'w'))
    plotNames(data, [n[0] for n in top_correlations + lowest_correlations], 
        smooth=3)
    
