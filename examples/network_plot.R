library(igraph)
devtools::load_all("~/projects/NeuronModel/mnfa")
res <- readRDS("~/projects/NeuronModel/serverjobs/Apr9/TH2015results_cluster.rds")
fit <- res$fit
n_cell <- 46
n_factor <- 6
L <- varimax(get_FA_estim(fit, n_cell=n_cell, n_factor=n_factor)$L)$loadings[1:n_cell,]
nodes <- data.frame(id=c(paste0("N", 1:n_cell), paste0("F", 1:n_factor)),
                    type=c(rep("neuron", n_cell), rep("factor", n_factor)))
edges <- data.frame(neuron=paste0("N", rep(1:n_cell, each=n_factor)), 
		    factor=paste0("F", rep(1:n_factor, n_cell)),
                    loading=as.vector(t(L)))
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)

# Adjust edge width based on loading value
E(net)$width <- E(net)$loading*10

# Set different colors for different nodes
colors <- c("blue", "gold")
V(net)$color <- colors[as.numeric(as.factor(V(net)$type))]

# Specify layout
mylayout <- layout_as_tree(net, root=which(nodes$type=="factor"))
plot(net, layout=mylayout)
