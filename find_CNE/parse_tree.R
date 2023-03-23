argv<-commandArgs();
library(phytools);
tree1<-read.tree(file=argv[6]);# modified read.newick to read.tree function for phytools 1.0.1 2022.3.14
parents_nodes_lab<-tree1$node.label;
child_nodes_lab<-tree1$tip.label;
edge<-tree1$edge;
file_parents<-paste(argv[6],".pnodes",sep="");
file_child<-paste(argv[6],".cnodes",sep="");
file_edge<-paste(argv[6],".edge",sep="");
write.table(parents_nodes_lab,file=file_parents,sep='\t',quote = F,row.names = F,col.names = F);
write.table(child_nodes_lab,file=file_child,sep='\t',quote = F,row.names = F,col.names = F);
write.table(edge,file=file_edge,sep='\t',quote = F,row.names = F,col.names = F);

