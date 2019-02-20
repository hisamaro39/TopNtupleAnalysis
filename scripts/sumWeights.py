
import ROOT

input_file = open('input.txt')
sw = {}
c = ROOT.TChain("sumWeights")
for line in input_file:
  line = line[0:-1]
  print "On file ", line
  c.Add(line)

for entry in xrange(c.GetEntries()):
    c.GetEntry(entry)
    print "Entry ", entry, " of ", c.GetEntries()-1
    if not c.dsid in sw:
      sw[c.dsid] = 0
    sw[c.dsid] = c.totalEventsWeighted + sw[c.dsid]
input_file.close()

for k in sw:
  print "{:10} -> {:10}".format(k, sw[k])
