#!/usr/bin/env python3
"""
Parsing support for howde tree hierarchy files.

The tree hierarchy file looks something like this (below).  It is a pre-order
traversal of the tree, and the asterisks indicate the depth of each node. See
BloomTree::read_topology() in bloom_tree.cc for a more complete description.
  node5169.bf
  *node5160.bf
  **SRR765.bf
  **node5146.bf
  ***SRR364.bf
  ***node5124.bf
  ****SRR767.bf
  ****SRR673.bf
  ****SRR211.bf
  *SRR769.bf
   ...
"""

from sys import argv,stdin,stdout,stderr,exit
from os  import path as os_path


def read_howde_tree_file(f,keepFileExtension=False,keepTags=False,debug=False):
	nodeStack = []
	topLevelSame = False
	topLevel = None
	forest = []

	for (lineNumber,level,name,tags) in read_howde_list(f,keepFileExtension=keepFileExtension):
		if (debug):
			print("read %d %s" % (level,name),file=stderr)
			for (dgbLevel,dbgNode) in nodeStack:
				print("  %d %s" % (dgbLevel,dbgNode.name),file=stderr)

		if (not keepTags) and (tags != None):
			exit("%s: line %d of tree input contains extra fields (\"%s\")"
			   % (os_path.basename(argv[0]),lineNumber,tags[0]))

		if (level == 0) and (topLevel == 0):
			assert (len(nodeStack) == 1)
			(treeLevel,tree) = nodeStack.pop()
			assert (treeLevel == 0)
			forest += [tree]
			if (debug):
				print("adding %s to forest" % tree.name,file=stderr)

		while (topLevelSame) and (level < topLevel):
			(sibLevel,sibling) = nodeStack.pop()
			siblings = [sibling]
			if (level < 0) and (nodeStack == []):
				assert (sibLevel == 0)
				forest += [sibling]
				if (debug):
					print("adding %s to forest" % sibling.name,file=stderr)
				break
			(peekLevel,_) = nodeStack[-1]
			while (peekLevel == sibLevel):
				(_,sibling) = nodeStack.pop()
				siblings += [sibling]
				(peekLevel,_) = nodeStack[-1]
			siblings.reverse()
			(parentLevel,parent) = nodeStack[-1]
			assert (parentLevel == topLevel-1)

			if (debug):
				print("  assigning %s children %s"
				    % (parent.name,",".join([child.name for child in siblings])),
				     file=stderr)
			parent.children = siblings
			for sibling in siblings:
				sibling.parent = parent

			topLevel = parentLevel
			if (len(nodeStack) > 1):
				(prevLevel,_) = nodeStack[-2]
				topLevelSame = (prevLevel == topLevel)
			else:
				topLevelSame = False
				(rootLevel,tree) = nodeStack.pop()
				assert (rootLevel == 0)
				forest += [tree]
				if (debug):
					print("adding %s to forest" % tree.name,file=stderr)

		if (level < 0): break   # (end-of-list)

		if (debug):
			print("pushing %d %s" % (level,name),file=stderr)
		node = TreeNode(name)
		if (keepTags): node.tags = tags
		nodeStack += [(level,node)]
		topLevelSame = (level == topLevel)
		topLevel = level

	if (debug):
		print("returning [%s]"
		    % ",".join([tree.name for tree in forest]),
		      file=stderr)

	return forest


def read_howde_list(f,keepFileExtension=False):
	lineNumber = 0
	for line in f:
		lineNumber +=1 
		line = line.strip()

		indent = 0
		while (line[indent] == "*"):
			indent += 1

		line = line[indent:].split()
		name = line[0]
		if ("/" in name): name = name.split("/")[-1]
		if (not keepFileExtension): name = name.split(".bf")[0]

		if (len(line) == 1): yield (lineNumber,indent,name,None)
		else:                yield (lineNumber,indent,name,line[1:])

	yield (-1,-1,"end-of-list",None)


class TreeNode(object):

	def __init__(self,name,parent=None,children=None):
		self.name = name
		if (children == None): self.children = []
		else:                  self.children = list(children)
		self.parent = parent
		self.depth  = None

	def compute_depth(self,depth=1):  # root has depth 1
		self.depth = depth
		for child in self.children:
			child.compute_depth(depth+1)

	def build_dict(self,nameToNode=None):
		if (nameToNode == None): nameToNode = {}
		nameToNode[self.name] = self
		for child in self.children:
			child.build_dict(nameToNode)
		return nameToNode

	def pre_order(self,order=None):
		if (order == None): order = [self]
		else:               order += [self]
		for child in self.children:
			child.pre_order(order)
		return order

	def list_pre_order(self,f=None,indent=0,fileSpec=None):
		if (f == None): f = stdout
		name = self.name
		if (fileSpec != None):
			name = fileSpec.replace("{name}",name)
		print("%s%s" % ("*"*indent,name),file=f)
		for child in self.children:
			child.list_pre_order (f,indent=indent+1,fileSpec=fileSpec)

	def post_order(self,order=None):
		for child in self.children:
			child.post_order(order)
		if (order == None): order = [self]
		else:               order += [self]
		return order

	def list_post_order(self,f=None,indent=0,fileSpec=None):
		if (f == None): f = stdout
		for child in self.children:
			child.list_post_order (f,indent=indent+1,fileSpec=fileSpec)
		name = self.name
		if (fileSpec != None):
			name = fileSpec.replace("{name}",name)
		print("%s%s" % ("*"*indent,name),file=f)

	def list_leaf_groups(self,f=None):
		if (f == None): f = stdout

		allChildrenAreLeafs = True
		for child in self.children:
			childIsLeaf = (child.children == [])
			if (not childIsLeaf):
				allChildrenAreLeafs = False
				break

		if (allChildrenAreLeafs) and (self.children != []):
			print(" ".join([child.name for child in self.children]),file=f)
		else:
			for child in self.children:
				child.list_leaf_groups(f)

	def compute_height_etc(self):
		self.height      = 1   # (longest path to a leaf, including self)
		self.subtreeSize = 1   # (nodes and leaves in subtreee, including self)
		for child in self.children:
			(childHeight,childSubtreeSize) = child.compute_height_etc()
			self.height = max(self.height,childHeight+1)
			self.subtreeSize += childSubtreeSize
		return (self.height,self.subtreeSize)

	def list_height_etc(self,f=None,indent=0):
		if (f == None): f = stdout
		print("%s%s\theight=%d\tdescendants=%d"
		    % ("*"*indent,self.name,self.height,self.subtreeSize-1),
		      file=f)
		for child in self.children:
			child.list_height_etc(f,indent=indent+1)

