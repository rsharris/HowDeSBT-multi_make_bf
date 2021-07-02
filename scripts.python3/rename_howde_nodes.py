#!/usr/bin/env python3
"""
Rename nodes an howdesbt tree, so that higher nodes have lower node numbers.
"""

from sys              import argv,stdin,stdout,stderr,exit
from howde_tree_parse import read_howde_tree_file


def usage(s=None):
	message = """
usage: cat howde_tree_file | rename_howde_nodes <template> [options]
  <template>   template for names to assign to internal tree nodes; this must
               contain the substring {node}

Rename nodes an howdesbt tree file, so that higher nodes have lower node
numbers. This only gives the nodes new names in the tree file -- it doesn't
rename any bf files that have been created.

Note that this functionality is now performed by default in "howdesbt cluster"."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	nameTemplate = None
	debug        = {}

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if ("{node}" in arg):
			nameTemplate = arg
		elif (arg == "--debug"):
			debug["debug"] = True
		elif (arg.startswith("--debug=")):
			for name in argVal.split(","):
				debug[name] = True
		else:
			usage("unrecognized option: %s" % arg)

	if (nameTemplate == None):
		usage("you must give me a template for the new node names")

	# read the tree

	forest = read_howde_tree_file(stdin,keepFileExtension=True,debug=("parse" in debug))
	assert (len(forest) != 0), "input has no tree"

	# collect nodes by depth

	depthToNodes = {}

	for tree in forest:
		tree.compute_depth()
		preOrder = tree.pre_order()

		for node in preOrder:
			if (node.depth not in depthToNodes): depthToNodes[node.depth] = []
			depthToNodes[node.depth] += [node]

	# assign new names

	nodeNumber = 0
	for depth in depthToNodes:
		for node in depthToNodes[depth]:
			if (len(node.children) == 0): continue
			nodeNumber += 1
			node.name = nameTemplate.replace("{node}",str(nodeNumber))

	# output the tree topology, with the new names

	for tree in forest:
		tree.list_pre_order()


if __name__ == "__main__": main()
