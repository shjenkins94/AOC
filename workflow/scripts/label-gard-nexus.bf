RequireVersion ("2.5.74");

LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");


LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");

labeler.analysis_description = {
  terms.io.info : "Read a NEXUS file from GARD and a list of sequence names, and annotate trees in the tree block for subsequent analyses",
  terms.io.version :          "0.1",
  terms.io.reference :        "TBD",
  terms.io.authors :          "Scott H. Jenkins",
  terms.io.contact :          "shjenkins94@gmail.com",
  terms.io.requirements :     "A NEXUS file from GARD"
};

io.DisplayAnalysisBanner (labeler.analysis_description);

KeywordArgument ("msa", "The GARD output file to annotate (NEXUS)");
SetDialogPrompt ("The GARD output file to annotate (NEXUS)");

fscanf (PROMPT_FOR_FILE, "Raw", labeler.raw_string);

// Ends up being easier to extract the trees from the raw string than to try to reconstruct them from the tree objects, so we'll do that and then load the trees from the newick strings we extract
labeler.tree_strings = extract_tree_strings(labeler.raw_string);

labeler.partition_count = Abs(labeler.tree_strings);

labeler.trees = utility.Map(labeler.tree_strings, "_tree_string_", 'trees.LoadAnnotatedTopology(_tree_string_)');
labeler.ts = utility.Map(labeler.trees, "_tree_", 'set_topology(_tree_ [^"terms.trees.newick_with_lengths"])');
labeler.existing = utility.Map(labeler.trees, "_tree_", '_tree_ [terms.trees.model_map]');

console.log ("Loaded dataset with " + labeler.partition_count + " partitions.");

KeywordArgument  ("regexp", "Use the following regular expression to select a subset of leaves", "()");
labeler.regexp = io.PromptUserForString ("Use the following regular expression to select a subset of leaves");

KeywordArgument  ("label", "Use the following label for annotation", "Foreground");
labeler.tag = io.PromptUserForString ("Use the following label for annotation");

KeywordArgument  ("reroot", "Reroot the tree on this node ('None' to skip rerooting)", "None");
labeler.reroot = io.PromptUserForString ("Reroot the tree on this node ('None' to skip rerooting)");

// Reroot the trees if requested
// Test later
if (labeler.reroot != "None") {
  labeler.ts = utility.Map(labeler.ts, "_topology_", 'set_topology(RerootTree (_topology_, labeler.reroot))');  
}

KeywordArgument  ("invert", "Invert selection", "No");
labeler.inverse  = io.SelectAnOption (
  {
    {"No","Matching list/regexp"}
    {"Yes","NOT matching list/regexp"}
    
  },
  "Invert selection"
) == "Yes";

if (labeler.regexp == "()") {
  KeywordArgument  ("list", "Line list of sequences to include in the set (required if --regexp is not supplied)");
  SetDialogPrompt  ("Line list of sequences to include in the set");
  labeler.list_dict = io.ReadDelimitedFile (null,",",FALSE);
}


KeywordArgument  ("internal-nodes", "Strategy for labeling internal nodes", "All descendants");
labeler.kind  = io.SelectAnOption (
  {
    {"None","Only assign labels to selected nodes"}
    {"All descendants","Only label an internal node if all its descendants are labeled"}
    {"All descendants, no MRCA","Only label an internal node if all its descendants are labeled"}
    {"Some descendants","Only label an internal node if some of its descendants are labeled"}
    {"Parsimony","Use maximum parsimony to label internal nodes"}
  },
  "Strategy for labeling internal nodes"
);

KeywordArgument  ("leaf-nodes", "Strategy for labeling selected leaves", "Label");

labeler.tips  = io.SelectAnOption (
  {
    {"Label","Include tips"}
    {"Skip","Only label internal nodes using the selected strategy"}
    
  },
  "Include tips in labeling"
) == "Label";

labeler.lists = utility.Map(labeler.ts, "_topology_", 'get_labels(_topology_)');

// let's get a version of the original nexus that cuts off before the trees block, so we can append our new trees block to it
labeler.header_string = regexp.Replace (labeler.raw_string, "(BEGIN TREES;.*)", "");

KeywordArgument ("output", "Write labeled NEXUS tree to");
labeler.outpath = io.PromptUserForFilePath ("Write labeled Newick tree to");

// now let's print the header string to the output file, and then append the new trees block to it
fprintf (labeler.outpath, CLEAR_FILE, labeler.header_string);
fprintf (labeler.outpath, "BEGIN TREES;\n");
for (p = 0; p < labeler.partition_count; p += 1) {
  labeler.current_topology = labeler.ts[p];
  labeler.current_labels = labeler.lists[p];
  labeler.current_existing = labeler.existing[p];

  labeled_topology = tree.Annotate("labeler.current_topology", "relabel_and_annotate", "{}", TRUE);

  fprintf (labeler.outpath, "\tTREE tree_", p + 1, " = ", labeled_topology, ";\n");
}

fprintf (labeler.outpath, "END;\n");


/**
 * Extracts the tree strings from a NEXUS file as a list of Newick strings. Assumes that the trees are in a block that starts with "BEGIN TREES;" and ends with "END;", and that each tree is on its own line and starts with "\tTREE (tree name) = ".
 * @name extract_tree_strings
 * @param {String} nexus_string - the NEXUS file as a string
 * @returns {AssociativeList} - a list of Newick strings
 */
function extract_tree_strings(nexus_string) {
  tree_block = regexp.FindSubexpressions (nexus_string, "(\tTREE .*?)\nEND;")[1];
  tree_lines = regexp.Split (tree_block, "\n");
  tree_strings = utility.Map(tree_lines, "_line_", 'regexp.Replace(_line_, "\tTREE (.*?) = ", "")');
  return tree_strings;
}

/**
 * helper function for setting tree topologies.
 * @name set_topology
 * @param {String} tree_string - a Newick string representing the tree topology
 * @returns {String} Topology - a Topology object representing the tree topology
 */
function set_topology(tree_string) {
  Topology T = tree_string;
  return T;
}

/**
 * helper function for setting tree topologies.
 * @name get_labels
 * @param {String} Topology - a Topology object representing the tree topology
 * @returns {AssociativeList} - a list of nodes to be labeled
 */
function get_labels(topology) {
  labels = {};

  if (labeler.regexp == "()") {
    label_list = {};

    for (k, v; in; labeler.list_dict["rows"]) {
      label_list[regexp.Replace (v[0],"\\ +$", "")] = 1;
    }

    if (labeler.inverse) {
      for (n; in; topology) {
        if (label_list[n] == 0) {
          labels[n] = labeler.tag;
        }
      }
    } else {
      for (n; in; topology) {
        if (label_list[n]) {
          labels[n] = labeler.tag;
        }
      }
    }
  } else {
    if (labeler.inverse) {
      for (n; in; topology) {
        if (None == regexp.Find (n,labeler.regexp)) {
          labels[n] = labeler.tag;
        }
      }
    } else {
      for (n; in; topology) {
        if (regexp.Find (n,labeler.regexp)) {
          labels[n] = labeler.tag;
        }
      }
    }
  }

  // might be fine if some trees don't have labels.
  if (utility.Array1D (labels) == 0) {
    console.log ("No labels found for tree " + topology);
  } else {
      label_num = utility.Array1D (labels);

      console.log ("\nSelected " + label_num + " branches to label for " + topology + "\n");
    
      if (labeler.tips) {

        internal_node_labels = {};
      
        if (labeler.kind == "All descendants" || labeler.kind == "All descendants, no MRCA" ) {
          internal_node_labels  =  ((trees.ConjunctionLabel ("topology", labels))["labels"]);
        }

        if (labeler.kind == "Some descendants") {
          internal_node_labels = ((trees.DisjunctionLabel ("topology", labels))["labels"]);
        }

        if (labeler.kind == "Parsimony") {
          internal_node_labels = ((trees.ParsimonyLabel ("topology", labels))["labels"]);
        }
        
        if (labeler.kind == "All descendants, no MRCA") {
            /* 
                here, we are going to perform a POST-order traversal of the tree, and remove labels for the internal nodes that are closest to the root,
                i.e. nodes that are roots of subtrees
            */
            
            rem_label = {};
            
            for (n,p; in; trees.ParentMap ("topology")) {
                if (p) {
                    if (internal_node_labels / n) {
                        if (internal_node_labels[n] != internal_node_labels[p]) {
                            rem_label + n;
                        }
                    }
                }
            }   
            
            for (n; in;  rem_label) {
                internal_node_labels - n;
            }
        }
        labels * internal_node_labels;
      
        console.log ("\nLabeled " + (utility.Array1D (labels) - label_num) + " additional internal branches for " + topology + "\n");
      } else {
          if (labeler.kind == "All descendants") {
            labels = ((trees.ConjunctionLabel ("topology", labels))["labels"]);
          }

          if (labeler.kind == "Some descendants") {
            labels = ((trees.DisjunctionLabel ("topology", labels))["labels"]);
          }

          if (labeler.kind == "Parsimony") {
            labels = ((trees.ParsimonyLabel ("topology", labels))["labels"]);
          }
          console.log ("\nLabeled " + (utility.Array1D (labels) ) + " internal branches for " + topology + "\n");
      }
  }
  return labels;
}


function relabel_and_annotate (node_name) {
    _label = "";
    if (Abs(labeler.current_existing [node_name]) > 0 ) {
        _label = "{" + labeler.current_existing [node_name] + "}";
    } else {
        if (labeler.current_labels / node_name) {
            _label = "{" + labeler.current_labels[node_name] + "}";
        }
    }
    return node_name + _label;
}