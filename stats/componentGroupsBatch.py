import sys
import os
sys.path.insert(0, '.')
sys.path.insert(0, os.path.join('.','SBMLparser'))
import argparse
import fnmatch
import progressbar
import radarChart
import numpy as np
from SBMLparser.rulifier import componentGroups
from collections import defaultdict

from sklearn.cluster import AffinityPropagation
import cPickle as pickle
from collections import Counter

import stats
def getFiles(directory,extension):
    """
    Gets a list of bngl files that could be correctly translated in a given 'directory'

    Keyword arguments:
    directory -- The directory we will recurseviley get files from
    extension -- A file extension filter
    """
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.{0}'.format(extension)):
            filepath = os.path.abspath(os.path.join(root, filename))
            matches.append([filepath,os.path.getsize(os.path.join(root, filename))])

    #sort by size
    #matches.sort(key=lambda filename: filename[1], reverse=False)
    
    matches = [x[0] for x in matches]

    return matches


def getContextRequirements(files):
    componentGroupsDict = {}
    progress = progressbar.ProgressBar()
    for fileidx in progress(range(len(files))):
        fileStr = files[fileidx]
        frequencyCounter = defaultdict(float)
        moleculeDict =  componentGroups.getContextRequirements(fileStr)
        total = 0
        for molecule in moleculeDict:

            counter =  {x:len(moleculeDict[molecule][x]) for x in moleculeDict[molecule]}
            for element in counter:

                frequencyCounter[element] += counter[element]
                total += counter[element]

        for element in frequencyCounter:
            frequencyCounter[element]/= total
        componentGroupsDict[fileStr] = frequencyCounter
    return componentGroupsDict


def getClusters(dependencyDict):
    names = []
    datapoints = []
    dimensions = ['independent','requirement','nullrequirement','exclusion']
    for filename in dependencyDict:
        names.append(filename.split('/')[-1])
        datapoints.append([dependencyDict[filename][dimension] for dimension in dimensions])
    datapoints = np.array(datapoints)

    af = AffinityPropagation().fit(datapoints)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_

    #n_clusters_ = len(cluster_centers_indices)
    #centers = datapoints[cluster_centers_indices]
    clusterDict = defaultdict(list)
    for name,label in zip(names,labels):
        clusterDict[label].append(name)
    return clusterDict,datapoints[cluster_centers_indices]
    
def getAnnotationsPerCluster(clusters,modelAnnotations):
    annotationClusters = defaultdict(Counter)
    for cluster in clusters:
        for element in clusters[cluster]:
            for element in modelAnnotations['XMLExamples/curated/{0}'.format('.'.join(element.split('.')[:-1]))]:
                if 'taxonomy' in element or 'mamo' in element or 'doid' in element:
                    continue
                annotationClusters[cluster][element] +=1
    return annotationClusters

if __name__ == "__main__":

    print getContextRequirements(['curated/output34.xml'])
    raise Exception
    #bngxml = getFiles('curated/atomized','xml')
    #dependencyDict = getContextRequirements(bngxml)
    #with open('dependencyDict.dump','wb') as f:
    #    pickle.dump(dependencyDict,f)

    with open('dependencyDict.dump','rb') as f:
        dependencyDict = pickle.load(f)

    with open('XMLExamples/curated/modelAnnotationDictionary.dump','rb') as f:
        modelAnnotations = pickle.load(f)

    clusters, centers = getClusters(dependencyDict)
    clusterAnnotations = getAnnotationsPerCluster(clusters,modelAnnotations)
    topCandidates = [clusterAnnotations[x].most_common(5) for x in clusterAnnotations]
    print len(topCandidates)

    
    for x in clusters:

        if len(clusters[x]) > 200:
            print '++++++',clusters[x]

    for cluster,center,rawclusters in zip(topCandidates,centers,clusters):
        print '-----'
        print 'centers:',zip(center,['independent','requirement','nullrequirement','exclusion'])
        print 'num elements:',len(clusters[rawclusters])
        annotations = []
        for annotation in cluster:
            try:
                _,annotationStr =  stats.resolveAnnotation(annotation[0])
            except:
                print annotation[0]
            if annotationStr == '':
                print annotation
            else:
                annotations.extend(annotationStr)
        print ', '.join(annotations)
