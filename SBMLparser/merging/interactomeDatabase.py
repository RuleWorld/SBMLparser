from utils import readBNGXML
import pickle
import networkx as nx
import marshal
import functools
import progressbar
from collections import Counter


def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = marshal.dumps([args, kwargs])
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer


@memoize
def getBondsFromFileName(fileName):
    constructedFilename = '../curated/{0}.xml'.format(fileName)
    moleculeTypes, rules, _ = readBNGXML.parseXML(constructedFilename)
    bonds = set([])
    for molecule in moleculeTypes:
        for component in molecule.components:
            name = [x.name for x in moleculeTypes if x.name.lower() == component.name]
            if len(name) > 0:
                bonds.add(tuple(sorted([molecule.name, name[0]])))
    return bonds


def createGraph(organism, organismInformation):
    fileNameBonds = {}
    if organism != 'Homo sapiens':
        return
    for fileName in set([y for x in organismInformation for y in x['fileName']]):
        try:
            fileNameBonds[fileName] = getBondsFromFileName(fileName)
        except IOError:
            continue
    graph = nx.Graph()
    nameDict = {}
    annotationSet = [(x['name'], x['annotation']) for x in organismInformation]

    #composition/modification nodes
    compositionEdges = []
    catalysisEdges = []
    for species in organismInformation:
        
        annotationIntersection = [[y for y in annotationSet if len(set([z[0] for z in x]).intersection(y[1])) > 0]
                                  for x in species['otherAnnotation'] if [y for y in annotationSet if len(set([z[0] for z in x]).intersection(y[1])) > 0]]
        if len(annotationIntersection) > 0:
            #annotationCounter = [Counter([z[1] for z in x]) for x in species['otherAnnotation']]
            pairIntersection = [[x for x in annotationIntersection if any(w[1].intersection(set([k[0] for k in y])) for w in x)]
                                for y in species['otherAnnotation']]

            #pairIntersection = [x for x in pairIntersection if x]
            if any(len(x) > 0 for x in pairIntersection):
                for intersection,otherAnnotation in zip(pairIntersection, species['otherAnnotation']):
                    for other in otherAnnotation:
                        if other[1] in ['BQB_HAS_PART', 'BQB_HAS_VERSION', 'BQB_IS_VERSION_OF']:
                            for annotationInter in intersection:
                                for specificIntersection in annotationInter:
                                    if other[0] in specificIntersection[1]:
                                        if '/'.join(specificIntersection[0]) != '/'.join(species['name']):
                                            if other[1] in ['BQB_HAS_PART']:
                                                compositionEdges.append(['/'.join(species['name']), specificIntersection[0])])
                                            elif other[1] in ['BQB_HAS_VERSION', 'BQB_IS_VERSION_OF']:
                                                catalysisEdges.append(['/'.join(species['name']), '/'.join(specificIntersection[0])])

        nodeName = '/'.join(species['name'])
        #nodeName = '/'.join(species['annotationName']) if '/'.join(species['annotationName']) != '' else '/'.join(species['name'])

        graph.add_node(nodeName, LabelGraphics={'text': nodeName}, graphics={'width': len(species['fileName'])})
        for name in species['name']:
            nameDict[name] = nodeName

    edgesCounter = Counter()
    for filebonds in fileNameBonds:
        for bonds in fileNameBonds[filebonds]:
            if bonds[0] in nameDict and bonds[1] in nameDict:
                bondPair = tuple(sorted([nameDict[bonds[0]], nameDict[bonds[1]]]))
            else:
                bondName1 = bonds[0] if bonds[0] not in nameDict else nameDict[bonds[0]]
                bondName2 = bonds[1] if bonds[1] not in nameDict else nameDict[bonds[1]]
                bondPair = tuple(sorted([bondName1, bondName2]))

            edgesCounter[bondPair] += 1
            graph.add_edge(bondPair[0], bondPair[1],
                           graphics={'fill': "#000000", 'width': edgesCounter[bondPair]}, LabelGraphics={'text': edgesCounter[bondPair]})
    for compositionEdge in compositionEdges:
        graph.add_edge(compositionEdge[0], compositionEdge[1],
                       graphics={'fill': "#800080", 'targetArrow': "standard", 'style': "dashed"})
    for catalysisEdge in catalysisEdges:
        graph.add_edge(catalysisEdge[0], catalysisEdge[1],
                       graphics={'fill': "#008000", 'targetArrow': "standard", 'style': "dashed"})

    nx.write_gml(graph, organism + '.gml')


if __name__ == "__main__":
    with open('results2.dump', 'rb') as f:
        speciesDict = pickle.load(f)
    progress = progressbar.ProgressBar(maxval=len(speciesDict)).start()
    for idx, organism in enumerate(speciesDict):
        percent_done = idx * 1.0 / len(speciesDict) * 100.0
        progress.update(percent_done)
        createGraph(organism, speciesDict[organism])
