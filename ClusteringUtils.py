import re

class Utils:
    def __init__(self):
        self.ind2node = {}
        self.node2ind = {}

    def loadNodes(self, nodeFile):
        """ Load node indices in nodes.txt file """
        i = 1
        for line in nodeFile:
            if line != "":
                line_s = line.strip()
                line_s = line_s.strip('"')
                self.ind2node[i] = line_s
                self.node2ind[line_s] = i
            i += 1

    def convertToNodes(self, df):
        """ Convert string entries in data frame to indices """
        return df.applymap(lambda x: self.ind2node[x])

    def convertToInds(self, df):
        """ Convert indices in data frame to strings """
        return df.applymap(lambda x: self.node2ind[x])

    def sortedTuples(self, s):
        """ Sort the tuples in the given set """
        return map(lambda x: tuple(sorted(x)), s)

    def convertToTuples(self, df):
        """ Convert edge dataframe to edge tuples for intersection testing """
        return df.apply(lambda x: tuple(sorted(x)), axis = 1).values

    def nx2dot(self, infname, outfname):
        prog = re.compile(r'(\d*) (\d*) {\'weight\': (\d*\.\d*)}')

        with open(infname, 'r') as infile, open(outfname, 'w') as outfile:
            outfile.write('strict graph {\n')
            for line in infile:
                if line != "":
                    m = prog.match(line)
                    if m:
                        outfile.write('  ' + m.group(1) + ' -- ' + m.group(2) + ' [weight = ' + m.group(3) + ']\n')

            outfile.write('}\n')

    def edge2dot(self, infname, outfname):
        with open(infname, 'r') as infile, open(outfname, 'w') as outfile:
            outfile.write('strict graph {\n')
            for line in infile:
                if line != "":
                    m = line.strip('\n')
                    m = m.split('\t')
                    if m:
                        outfile.write('  ' + m[0] + ' -- ' + m[1] + ' [weight = ' + m[2] + ']\n')

            outfile.write('}\n')

    def writeEdgeFile(self, g, outfname):
        """ Write edge file from graph-tools graph g """
        weight = g.edge_properties['weight']
        names = g.vertex_properties['vertex_name']

        with open(outfname, 'w') as outfile:
            for e in g.edges():
                outfile.write('{} {} {{\'weight\': {}}}\n'.format(names[e.source()], names[e.target()], weight[e]))
