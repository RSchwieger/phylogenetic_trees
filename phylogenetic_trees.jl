### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 44972213-3979-4392-a208-173a35dea09d
using Combinatorics, Plots, Graphs, GraphRecipes

# ╔═╡ f5cab10a-23b2-11ed-0b96-8fc0acb2704b
md"""
# Phylogenetic Trees
"""

# ╔═╡ fab623cf-b40d-48bf-90c2-4d76554b0420
begin
	#uncomment to install
	#import Pkg
	#Pkg.add("Combinatorics")
	#Pkg.add("Plots")
	#Pkg.add("Graphs")
	#Pkg.add("GraphRecipes")
end

# ╔═╡ 6fa9071a-4b5b-45ce-bba1-d9249dc58587
md"""
## Background
Notation and ideas are based on the book *Recombinatorics*. In this notebook is explained how a phylogenetic tree is constructed from a set of input taxa.

**Assumptions**
- there exists a common ancestor
- **infinite-site model:** Each character mutates from the zero state to the one state **exactly once**.
"""

# ╔═╡ a4593138-3f6f-4a99-8ea0-1df66434e9c7
md"""
**Description of the input:**

- The set of taxa for which we want to construct the phylogenetic tree is written in  the form of a binary matrix.
- Its rows correspond to taxa.
- Its columns correspond to sites/characters (e.g. SNPs).
- The entry in the $i$-th row and $j$-th column of the input matrix indicates if a taxon $i$ possesses the $j$s character.
- 1 stands for "has the character".
- The common ancestor is without loss of generality the zero vector $0...0$.

For example:
"""

# ╔═╡ 2a62ea17-f928-4357-9e96-03be3952b56b
Array{Bool,}([0 1 0 0; 1 1 0 1; 0 0 0 1])

# ╔═╡ 23032982-cebc-4299-b131-368652eac0c6
md"""
- The first taxon possesses character 2.
- The second taxon possesses character 1,2,4.
- The third taxon possesses character 4.
"""

# ╔═╡ 6c83d66f-b25a-4e18-bce0-e061d1f782a1
md"""
#### Notation

- I write $ij$ for (i,j) e.g. $01$ for (0,1) 
- T denotes the phylogenetic tree
- M denotes the input matrix (rows correspond to taxa, columns correspond to sites)
"""

# ╔═╡ dcde1358-00a0-4a7d-b9ac-864086dba3d9
md"""
## The Perfect Phylogeny Problem
"""

# ╔═╡ b9815030-6516-4a96-aa37-50be6e53e8ab
md"""
**The perfect phylogeny problem:** Given $n$ by $m$ binary matrix $M$, determine whether there is a perfect phologeny for $M$, and if so, build one.
"""

# ╔═╡ f24b674f-7db8-4794-ad1e-6222d366c85e
md"""
**Definition (perfect phylogeney):** Given an $n \times m$ binary matrix $M$ for $n$ taxa, a **perfect phylogeny** for $M$ is a rooted directed tree with $n$ leaves, which satisfies:

1. Each of the taxa labels exactly one leaf.
2. Each of the characters labels exactly one edge.
3. For any taxon $f$, the characters that label the edges along the path from the root to the leaf $f$, specify the characters of $f$.
"""

# ╔═╡ 2e79a12d-f9e1-4c7a-9c38-f01b5302bd90
md"""
For example:
"""

# ╔═╡ b6bb3b31-653a-48ab-bf90-5ba3a632fa25
#Array{Bool,}([0 1 0 0; 1 1 0 1; 0 0 1 1])

# ╔═╡ 73dd5022-873e-4a23-ab47-a0f6ee6db5c4
Array{Bool,}([1 0 0; 1 1 0; 1 1 1])

# ╔═╡ eef9bc48-e205-4c6d-ab89-bbdb31954b32
begin
	g = DiGraph(7)
	add_edge!(g, 1, 2)
	add_edge!(g, 2, 3)
	add_edge!(g, 3, 7)
	add_edge!(g, 2, 4)
	add_edge!(g, 3, 5)
	add_edge!(g, 5, 6)
	plot(g, names=["000", "100", "110", "100", "111", "111", "110"], curves=false, edgelabel=Dict((1,2) => "1", (2,3) => "2", (3,5) => "3"))
end

# ╔═╡ e0bf952c-76f8-4088-9fb5-36001c57ebb2
md"""
For simplicity, we make also the assumption that no rows and columns are identical.

We can test this with the following function.
"""

# ╔═╡ b00ac74f-d9b4-4c25-b14d-37befa4f786f
	"""
	Tests if the matrix has duplicate rows.
	
	**Example:**
	
	```
	m_dup = [0 1 1 0; 1 1 1 1; 1 0 0 0]
	m_no_dup = [0 1 0; 1 1 1; 1 0 0]
	has_no_duplicate_columns(m_dup), has_no_duplicate_columns(m_no_dup)
	```
	The first function call returns false and the second true.
	"""
function has_no_duplicate_columns(m::AbstractMatrix)::Bool
	# Is the cardinality of the set of column the same than its number of columns?
	length(Set([m[:,i] for i in 1:size(m,2)])) == size(m,2)
end

# ╔═╡ a0a3ba95-5b43-4dd1-9e75-cbcd84528fa3
function has_no_duplicate_rows(m::AbstractMatrix)::Bool
	# Is the cardinality of the set of column the same than its number of columns?
	length(Set([m[i,:] for i in 1:size(m,1)])) == size(m,1)
end

# ╔═╡ 51434ae7-5eda-447a-a01c-fa9b694e7fd6
m_dup = [0 1 1 0; 1 1 1 1; 1 0 0 0]

# ╔═╡ fbd915b6-97bf-402d-a10a-4adfeac658f0
m_no_dup = [0 1 0; 1 1 1; 1 0 0; 1 0 0]

# ╔═╡ 568e0ed7-ccfc-4190-ad3b-0517386c6bcc
begin
	# Testcases
	@assert ( has_no_duplicate_rows(m_dup) )
	@assert ( !has_no_duplicate_columns(m_dup) )
	@assert ( !has_no_duplicate_rows(m_no_dup) )
	@assert ( has_no_duplicate_columns(m_no_dup) )
end

# ╔═╡ 0e86b461-39f0-4097-a8f1-6630b246f5be
md"""
Let us assume in this notebook w.l.o.g. that the **ancestral sequence is the zero vector**. Then the four gametes theorem becomes:
"""

# ╔═╡ 24e7cee9-33e7-4392-8c70-44b8f00c6e30
md"""
**The perfect-phylogeny theorem (p. 60)**: Matrix $M$ has a perfect phylogeny (with all-zero ancestral sequence) iff no pair of columns c,d contains the three binary pairs $01$, $10$ and $11$.
"""

# ╔═╡ de8a56f9-3650-4263-9867-fb6fc5912933
md"""
Illustration:
"""

# ╔═╡ f17bd807-ada1-4209-acbe-f913751e3261
[1 1; 0 1; 1 0]

# ╔═╡ 87a56a18-55e7-45ce-91ff-512120b747fb
begin
	g2 = DiGraph(4)
	add_edge!(g2, 1, 2)
	add_edge!(g2, 1, 3)
	add_edge!(g2, 3, 4)
	add_edge!(g2, 2, 4)
	plot(g2, names=["00", "10", "01", "11"], curves=false, edgelabel=Dict((1,2) => "1", (1,3) => "2", (3,4) => "?", (2,4) => "?"), method=:directed_tree)
end

# ╔═╡ b08d93c6-72c9-45ee-aef7-2e496aad3bdc
md"""
If a pair of columns/characters contain the pairs $01$, $10$, $11$ one of the two characters needs to change at least two times to explain the occurrence of all three taxa.
"""

# ╔═╡ 0cbf7248-c762-4d1f-81b6-a71f2e743a5c
md"""
We can test if a matrix has a perfect phylogeny with the following two functions:
"""

# ╔═╡ 5806317c-f81a-4220-aa9c-ad8e39091b05
begin
	test_set_for_perfect_phylogeny = Set([(0,1), (1,0), (1,1)])
	"""
	Tests if two columns are compatible. That is they do not contain simultanously the
	pairs 01, 10, 11
	
	**Example:**
	
	```
	m = [0 1 1 0; 1 1 1 1; 1 0 0 0]
	contruct_phylogenetic_tree(m, (1, 2))
	contruct_phylogenetic_tree(m, (2, 4))
	```
	The first function call returns false and the second true.
	"""
	function are_columns_compatible(m::AbstractMatrix{Bool,}, pair_of_indices)::Bool
		Set(zip(m[:, pair_of_indices[1]], m[:, pair_of_indices[2]])) != test_set_for_perfect_phylogeny
	end
end

# ╔═╡ 5108ea79-f805-490c-a855-1b216f8a3e33
testcase4are_columns_compatible = Array{Bool,}([0 1 1 0; 1 1 1 1; 1 0 0 0])

# ╔═╡ 24f21ee2-2fe0-4f9c-90c7-6cfa8356c30b
md"""
- The pair of columns (1,2) contains the pairs $\{01, 11, 10\}$. Therefore, the function evaluates to false.
- The pair of columns (1,4) contains the pairs $\{00, 11, 10\} \not = \{01, 10, 11\}$. Therefore, the function evaluates to true.
"""

# ╔═╡ bee7850d-ade4-4bcc-b0c0-f92bbd8d2df3
begin
	@assert( !are_columns_compatible(testcase4are_columns_compatible, (1, 2)) )
	@assert( are_columns_compatible(testcase4are_columns_compatible, (1, 4)) )
end

# ╔═╡ a9b8a97b-aa72-45c9-8f8a-ed30c1f0b0f3
"""
	Tests if the input matrix has only perfect characters.
	
	**Example:**
	
	```
	m1 = Array{Bool,}([0 1 1 0; 1 1 1 1; 1 0 1 0])
	m2 = Array{Bool,}([0 1 0 1; 1 1 0 1; 0 0 0 1])
	has_only_perfect_characters(m1)
	has_only_perfect_characters(m2)
	```
	The first function call returns false and the second true.
	"""
function has_phylogentic_tree(m::AbstractMatrix{Bool,})::Bool
	pairs_of_columns = combinations(1:size(m, 2), 2) # All pairs of columns
	reduce(&, map(x -> are_columns_compatible(m, x), pairs_of_columns))
end

# ╔═╡ f83df4a0-01b1-47d8-a186-850af7538e06
testcase4has_phylogentic_tree_1 = Array{Bool,}([0 1 1 0; 1 1 1 1; 1 0 1 0])

# ╔═╡ 2d62a44a-1d0f-4342-a3d1-9095437d749f
md"""
This matrix has **no perfect phologeny** since the pair of columns (1,2) is not compatible.
"""

# ╔═╡ 472ee63f-3bae-44c2-bbe2-1192047dd17d
testcase4has_phylogentic_tree_2 = Array{Bool,}([0 1 1; 1 1 1; 0 0 1])

# ╔═╡ 1f0f7f6c-8119-499a-8608-632e99af4158
md"""
This matrix has **a perfect phologeny** since:

- The pair of columns (1,2) contains the pairs $\{01, 11, 00\} \not = \{01, 10, 11\}$
- The pair of columns (1,3) contains the pairs $\{00, 11\} \not = \{01, 10, 11\}$
- ...
"""

# ╔═╡ d8ada293-ca93-4ec5-b8aa-630e07bb068c
begin
	@assert( !has_phylogentic_tree(testcase4has_phylogentic_tree_1) )
	@assert( has_phylogentic_tree(testcase4has_phylogentic_tree_2) )
end

# ╔═╡ 457deb47-d57c-489d-8a47-360ce65bfcb3
md"""
## Construction of a Phylogenetic Tree
"""

# ╔═╡ 7ec6bd75-f906-40d5-b71d-86730cff8ed9
md"""
In a first step we sort the input matrix columns by the number of ones they contain (largest first).
A sorted version of the matrix $M$ will be denoted by $\overline{M}$.
"""

# ╔═╡ 5d38aef5-5dfe-4ed6-b788-27379b4b3262
"""
	Sorts the columns of a binary input matrix and returns a sorted version of the input matrix. The columns of the sorted matrix are ordered in descending order with respect to the number of ones they possess.
	
	**Example**:
	
	```
	m = Array{Bool,}([0 1 1 0; 1 1 1 1; 0 0 0 1; 0 0 0 1])
	sorted_m = Array{Bool,}([0 1 1 0; 1 0 0 0; 1 0 0 0; 1 0 0 0])
	has_only_perfect_characters(m) == sorted_m
	```
	"""
function sort_binary_matrix(m::AbstractMatrix{Bool,})::AbstractMatrix{Bool,}
	# Returns a sorted version of the matrix
	# Sorting is done by number of ones
	sortslices(m, dims=2, lt=(x,y)->isless(sum(x),sum(y)), rev=true )
end

# ╔═╡ b24bb206-a9b2-424b-8ec8-87f9b2c06d2b
"""
	Checks if a binary matrix is sorted with respect to columnwise sums.
	
	**Example**:
	
	```
	random_binary_matrix = rand(Bool, (10,5))
	@assert( is_sorted(sort_binary_matrix(random_binary_matrix)) )
	```
	"""
function is_sorted(m::AbstractMatrix{Bool,})::Bool
	@assert( size(m, 2) > 0 )
	last = sum(m[:, 1])
	for i in 2:size(m, 2)
		next = sum(m[:, i])
		if next > last
			return false
		end
		last = next
	end
	return true
end

# ╔═╡ f3ba36c3-c9e6-4c0d-aace-12f825d2e843
example_unsorted_matrix = Array{Bool,}([0 1 1 0 0; 1 1 1 1 0; 0 0 0 1 0; 0 0 0 1 0])

# ╔═╡ 43adfdc8-8ddf-4ae9-91a7-961772e9db16
begin
	# Testcase
	random_binary_matrix = rand(Bool, (10,5))
	@assert( is_sorted(sort_binary_matrix(example_unsorted_matrix)) )
end

# ╔═╡ bb123979-63f9-43a1-9b92-ea1db15be535
md"""
Its sorted version:
"""

# ╔═╡ d1d5a24e-aaf3-42d7-be06-62d1de791a04
sort_binary_matrix(example_unsorted_matrix)

# ╔═╡ 85bbd942-c696-4481-9a99-dd6ebf3705a7
example4shared_prefix_property = Array{Bool,}([1 1 1 0; 1 1 1 1; 1 1 0 0; 1 0 0 0])

# ╔═╡ 95bbc5b0-9f29-4128-874c-57bb97fb2288
md"""
**The shared-prefix property (p.62)**: For two taxa $f$ and $g$ let $d$ be the largest index character (rightmost in $\overline{M}$ that taxa $f$ and $g$ both possess (i.e. where both have state $1$). Then, assuming no pair of columns contain all three binary pairs $01$, $10$, $11$ rows $f$ and $g$ in $\overline M$ must be indentical from column one (at the left end of $\overline{M}$) to column $d$.)
"""

# ╔═╡ 865fd300-b2dd-4858-a3ae-f39506a5233d
md"""
**Algorithm (p. 43)**

Process the rows (i.e. taxa) of $\overline{M}$ in order:

(1) create the root node of $T$

(2) Create a path from the root to the taxon $r_1$ (if taxon $r_1$ contains $t$ characters that path contains $t$ edges.)

(3) Let $T_f$ denote the intermediate tree containing all first $f$ taxa of $\overline{M}$ 

Then $T_{f+1}$ is constructed from $T_f$ as follows:
- examine characters that taxon $f+1$ possesses (left to right) and in parallel walk from the root of $T_f$ down the path in the tree, as long as the successive characters on the path match the successive charactrers that taxon $f+1$ possesses.
- let $c$ denote the last matched character on the walk.
- the walk ends at a node $v_{f+1}$ containing all the characters to the right c that taxon $f+1$ possesses (in the order that they appear in $\overline{M}$)
- create a new path out of $v_{f+1}$ containing all the characters to the right $c$ that taxon $f+1$ possesses, followed by an unlabeled edge to a leaf labeled $f+1$.

"""

# ╔═╡ 666eb7bc-25a4-4c37-81d6-b4b80385aa17
"""
Traverses the phylogenetic tree until the first mismatch is detected. The last node of the traversal still existing in
the graph and the first character/column which does not correspond to an edge in the graph is returned.

"""
function traverse_phylogenetic_tree_until_mismatch(phylogenetic_tree::SimpleDiGraph{Int64}, row::Vector{Bool}, label2node::Dict{Vector{Bool}, Int64})::Tuple{Vector{Bool}, Int64}
	org = zeros(Bool, length(row))
	for col in 1:length(row)
		# if no character go to the next one
		if !row[col]
			continue
		end
		
		# construct the successor node in the graph
		dest = copy(org)
		dest[col] = true

		# check if the node is already in the graph
		if !haskey(label2node, dest)
			return org, col
		end
		org = dest # update the current node
	end
	return org, length(row)
end

# ╔═╡ 27fbffe2-5e97-4762-9588-12b751a311b1
"""
Constructs the phylogenetic tree of an input matrix `unsorted_m`. Returns a dictionary
of the form: character => edge, where edge is a tuple of Boolean vectors.

Note that in a phylogenetic tree in the infinite site model each character corresponds
to exactly one labeled edge.

**Example:**

```
m = Array{Bool,}([0 1 1 0; 1 1 1 1; 0 0 0 1; 0 0 0 1])
contruct_phylogenetic_tree(m)
```

"""
function contruct_phylogenetic_tree(m::AbstractMatrix{Bool,})
	# check if assumpotions are met
	@assert( has_phylogentic_tree(m) )
	@assert( has_no_duplicate_columns(m) )
	@assert( has_no_duplicate_rows(m) )
	@assert( is_sorted(m) )

	phylogenetic_tree = SimpleDiGraph()

	# insert root node
	root = zeros(Bool, size(m, 2))
	add_vertex!(phylogenetic_tree)
	inserted_node = 1 # nodes are integers
	node2label = [root] # we save the node and its label
	label2node = Dict(root => inserted_node)
	edge2label = Dict()
	
	# Iterate over the rows/taxa of m
	for row in 1:size(m, 1)
		# Invariant property: At the beginning of each iteration phylogenetic_tree is 
		# a perfect phylogeny for the rows 1 to (row-1)
		
		# traverse tree until first mismatch
		last_node_in_tree, col_of_mismatch = traverse_phylogenetic_tree_until_mismatch(phylogenetic_tree, m[row,:], label2node)

		# Due to the shared prefix property the remaining characters of
		# the row/taxon label no edges in the tree constructed until now
		# This guarantees property 2 of a perfect phylogeny
		# (2) each of the characters labels exactly one edge
		
		# iterate over the remaining columns of m
		org = last_node_in_tree
		# For each remaining character of the taxon we will insert a new edge
		# That guarantees property 3
		# (3) For any taxon 'row', the characters that label the edges along the path
		# from the root to the leaf row, specify the characters of 'row'
		for col in col_of_mismatch:size(m, 2)
			# invariant property: At the beginning of each iteration the phylogenetic_tree is a tree
			
			# if this taxon has no character 'col' go to the next one
			if !m[row, col]
				continue
			end
			
			# construct the successor node which needs to be inserted into
			# the existing graph
			dest = copy(org)
			dest[col] = true

			# we insert a new node into the the phylogenetic tree
			add_vertex!(phylogenetic_tree)
			inserted_node = inserted_node+1
			append!(node2label, [dest]) # we label the node
			label2node[dest] = inserted_node
			# add a new edge from the last node to the newly created one.
			# Since we just created the new vertex this edge
			# is guaranteed to not induce any cycle
			# => The new graph is as well a tree
			add_edge!(phylogenetic_tree, label2node[org], label2node[dest])
			edge2label[(label2node[org], label2node[dest])] = col

			# update org
			org = dest
		end
		# After we inserted the new path we insert the taxon as leaf to the last inserted leaf.
		# this guarantees property (3)
		# (3) Each taxon labels exactly one leaf
		add_vertex!(phylogenetic_tree) # we insert a new vertex (the leaf)
		inserted_node = inserted_node+1
		dest = copy(org) # we make a copy since we append that to the list of labels
		node2label = append!(node2label, [dest]) # relate node to label
		# We do not add the new node to label2node
		# this guarantees that this leaf will remain a leaf in each iteration
		add_edge!(phylogenetic_tree, label2node[org], inserted_node)
	end
	phylogenetic_tree, node2label, label2node, edge2label
end

# ╔═╡ 1a94bcac-c1d5-4aee-a970-a220a78d3635
"""
Small helper function to convert a Boolean vector into a string.

**Example:**

```
example = [0 1 1 0; 1 1 1 1; 0 0 0 1; 0 0 0 1]
labels2string(example)
```
This returns ["0110", "1111", "0001", "0001"]
"""
labels2string = labels -> map(join, map(x -> convert(Vector{Int}, x), labels))

# ╔═╡ 02d97d4c-2a07-4993-ac7d-2136ec22fa25
begin
	example4contruct_phylogenetic_tree = Array{Bool,}([1 1 1 0; 1 1 1 1; 1 1 0 0; 1 0 0 0])

	@assert( has_phylogentic_tree(example4contruct_phylogenetic_tree) )
	@assert( has_no_duplicate_columns(example4contruct_phylogenetic_tree) )
	@assert( has_no_duplicate_rows(example4contruct_phylogenetic_tree) )
	example4contruct_phylogenetic_tree
end

# ╔═╡ b8bbcb17-cdda-4464-bbdd-4e477b594272
begin
	phylogenetic_tree_1, labels_1, nodes_1, edge_labels_1 = contruct_phylogenetic_tree(example4contruct_phylogenetic_tree)

	plot(phylogenetic_tree_1, names=labels2string(labels_1), curves=false, edgelabel=edge_labels_1)
end

# ╔═╡ e3dbe35c-2fd5-4937-b3a2-ade39fd02db3
begin
	example4contruct_phylogenetic_tree_2 = sort_binary_matrix(Array{Bool,}([0 0 1 0; 0 1 1 0; 1 0 0 1; 1 0 0 0]))

	@assert( has_phylogentic_tree(example4contruct_phylogenetic_tree_2) )
	@assert( has_no_duplicate_columns(example4contruct_phylogenetic_tree_2) )
	@assert( has_no_duplicate_rows(example4contruct_phylogenetic_tree_2) )
	example4contruct_phylogenetic_tree_2
end

# ╔═╡ df59f942-37dc-4440-a311-75eaf9a79022
begin
	phylogenetic_tree_2, labels_2, nodes_2, edge_labels_2 = contruct_phylogenetic_tree(example4contruct_phylogenetic_tree_2)

	plot(phylogenetic_tree_2, names=labels2string(labels_2), curves=false, edgelabel=edge_labels_2)
end

# ╔═╡ ac9ba96a-d5aa-463c-af9b-9c70d35b6bf0
md"""
## Complexity
"""

# ╔═╡ e27e9c6d-4b8e-4562-92f7-08c2f8d8861f
md"""
Since we iterate over the rows and columns we have a complexity of $O(\text{number of character} \cdot \text{number of taxa})$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
GraphRecipes = "bd48cda9-67a9-57be-86fa-5b3c104eda73"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
Combinatorics = "~1.0.2"
GraphRecipes = "~0.5.9"
Graphs = "~1.7.2"
Plots = "~1.31.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "5856d3031cdb1f3b2b6340dfdc66b6d9a149a374"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.2.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccd479984c7838684b3ac204b716c89955c76623"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "cf0a9940f250dc3cb6cc6c6821b4bf8a4286cf9c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.66.2"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2d908286d120c584abbe7621756c341707096ba4"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.66.2+0"

[[GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "a7a97895780dab1085a97769316aa348830dc991"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.3"

[[GeometryTypes]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "d796f7be0383b5416cd403420ce0af083b0f9b28"
uuid = "4d00f742-c7ba-57c2-abde-4428a4b178cb"
version = "0.8.5"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[GraphRecipes]]
deps = ["AbstractTrees", "GeometryTypes", "Graphs", "InteractiveUtils", "Interpolations", "LinearAlgebra", "NaNMath", "NetworkLayout", "PlotUtils", "RecipesBase", "SparseArrays", "Statistics"]
git-tree-sha1 = "1735085e3a8dd0e14020bdcbf8da9893a5508a3f"
uuid = "bd48cda9-67a9-57be-86fa-5b3c104eda73"
version = "0.5.9"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "a6d30bdc378d340912f48abf01281aab68c0dec8"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.2"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "MbedTLS", "Sockets"]
git-tree-sha1 = "c7ec02c4c6a039a98a15f955462cd7aea5df4508"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.8.19"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "2f0be365951a88dfb084f754005177e6dfb00ed0"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.4"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[NetworkLayout]]
deps = ["GeometryBasics", "LinearAlgebra", "Random", "Requires", "SparseArrays"]
git-tree-sha1 = "cac8fc7ba64b699c678094fa630f49b80618f625"
uuid = "46757867-2c16-5918-afeb-47bfcb05e46a"
version = "0.4.4"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "a19652399f43938413340b2068e11e55caa46b65"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.31.7"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "e7eac76a958f8664f2718508435d058168c7953d"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.3"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "dfec37b90740e3b9aa5dc2613892a3fc155c3b42"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.6"

[[StaticArraysCore]]
git-tree-sha1 = "ec2bd695e905a3c755b33026954b119ea17f2d22"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.3.0"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "8c6ac65ec9ab781af05b08ff305ddc727c25f680"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.12"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─f5cab10a-23b2-11ed-0b96-8fc0acb2704b
# ╟─fab623cf-b40d-48bf-90c2-4d76554b0420
# ╟─44972213-3979-4392-a208-173a35dea09d
# ╟─6fa9071a-4b5b-45ce-bba1-d9249dc58587
# ╟─a4593138-3f6f-4a99-8ea0-1df66434e9c7
# ╟─2a62ea17-f928-4357-9e96-03be3952b56b
# ╟─23032982-cebc-4299-b131-368652eac0c6
# ╟─6c83d66f-b25a-4e18-bce0-e061d1f782a1
# ╟─dcde1358-00a0-4a7d-b9ac-864086dba3d9
# ╟─b9815030-6516-4a96-aa37-50be6e53e8ab
# ╟─f24b674f-7db8-4794-ad1e-6222d366c85e
# ╟─2e79a12d-f9e1-4c7a-9c38-f01b5302bd90
# ╟─b6bb3b31-653a-48ab-bf90-5ba3a632fa25
# ╟─73dd5022-873e-4a23-ab47-a0f6ee6db5c4
# ╟─eef9bc48-e205-4c6d-ab89-bbdb31954b32
# ╟─e0bf952c-76f8-4088-9fb5-36001c57ebb2
# ╟─b00ac74f-d9b4-4c25-b14d-37befa4f786f
# ╟─a0a3ba95-5b43-4dd1-9e75-cbcd84528fa3
# ╟─51434ae7-5eda-447a-a01c-fa9b694e7fd6
# ╟─fbd915b6-97bf-402d-a10a-4adfeac658f0
# ╠═568e0ed7-ccfc-4190-ad3b-0517386c6bcc
# ╟─0e86b461-39f0-4097-a8f1-6630b246f5be
# ╟─24e7cee9-33e7-4392-8c70-44b8f00c6e30
# ╟─de8a56f9-3650-4263-9867-fb6fc5912933
# ╟─f17bd807-ada1-4209-acbe-f913751e3261
# ╟─87a56a18-55e7-45ce-91ff-512120b747fb
# ╟─b08d93c6-72c9-45ee-aef7-2e496aad3bdc
# ╟─0cbf7248-c762-4d1f-81b6-a71f2e743a5c
# ╟─5806317c-f81a-4220-aa9c-ad8e39091b05
# ╟─5108ea79-f805-490c-a855-1b216f8a3e33
# ╟─24f21ee2-2fe0-4f9c-90c7-6cfa8356c30b
# ╠═bee7850d-ade4-4bcc-b0c0-f92bbd8d2df3
# ╟─a9b8a97b-aa72-45c9-8f8a-ed30c1f0b0f3
# ╟─f83df4a0-01b1-47d8-a186-850af7538e06
# ╟─2d62a44a-1d0f-4342-a3d1-9095437d749f
# ╟─472ee63f-3bae-44c2-bbe2-1192047dd17d
# ╟─1f0f7f6c-8119-499a-8608-632e99af4158
# ╠═d8ada293-ca93-4ec5-b8aa-630e07bb068c
# ╟─457deb47-d57c-489d-8a47-360ce65bfcb3
# ╟─7ec6bd75-f906-40d5-b71d-86730cff8ed9
# ╟─5d38aef5-5dfe-4ed6-b788-27379b4b3262
# ╟─b24bb206-a9b2-424b-8ec8-87f9b2c06d2b
# ╠═43adfdc8-8ddf-4ae9-91a7-961772e9db16
# ╟─f3ba36c3-c9e6-4c0d-aace-12f825d2e843
# ╟─bb123979-63f9-43a1-9b92-ea1db15be535
# ╠═d1d5a24e-aaf3-42d7-be06-62d1de791a04
# ╠═85bbd942-c696-4481-9a99-dd6ebf3705a7
# ╟─95bbc5b0-9f29-4128-874c-57bb97fb2288
# ╟─865fd300-b2dd-4858-a3ae-f39506a5233d
# ╟─666eb7bc-25a4-4c37-81d6-b4b80385aa17
# ╟─27fbffe2-5e97-4762-9588-12b751a311b1
# ╟─1a94bcac-c1d5-4aee-a970-a220a78d3635
# ╟─02d97d4c-2a07-4993-ac7d-2136ec22fa25
# ╟─b8bbcb17-cdda-4464-bbdd-4e477b594272
# ╟─e3dbe35c-2fd5-4937-b3a2-ade39fd02db3
# ╟─df59f942-37dc-4440-a311-75eaf9a79022
# ╟─ac9ba96a-d5aa-463c-af9b-9c70d35b6bf0
# ╟─e27e9c6d-4b8e-4562-92f7-08c2f8d8861f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
