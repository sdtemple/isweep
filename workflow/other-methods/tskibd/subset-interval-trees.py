
import sys
import tskit

treesin, treesout, left, right = sys.argv[1:]
left = int(float(left))
right = int(float(right))

# Load the tree sequence from file
ts = tskit.load(treesin)

# Define the interval to keep
interval = tskit.Interval(left=left, right=right)

# Subset the tree sequence by the interval
sub_ts = ts.keep_intervals([interval])

# Write the new tree sequence to file
sub_ts.dump(treesout)
