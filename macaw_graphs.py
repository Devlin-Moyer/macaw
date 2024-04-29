# macaw_graphs.py
'''
Graph theory/NetworkX stuff
'''

import numpy as np
import networkx as nx
from macaw_utils import add_reaction_equations
import itertools as it
from multiprocessing import Pool
from macaw_utils import sigfig_round

def make_pathway_graph(
    model, exclude = list(), use_names = False, add_suffixes = False
):
    '''
    Group all reactions in the given Cobra.Model except those whose IDs are in
    `exclude` into a set of "pathways" by connecting reactions that both use at
    least one metabolite that is produced by one and consumed by the other,
    except for any metabolites whose IDs are in `exclude`.
    First, create a bipartite graph with all reactions and metabolites in the
    Model, where each node represents a reaction or a metabolite, and edges
    connect metabolites to the reactions that produce or consume them. Direct
    edges to distinguish consumption from production (reversible reactions have
    two edges connecting them to each of their metabolites, one in each
    direction)
    Then remove all nodes corresponding to metabolites or reactions whose IDs
    are in `exclude`
    Project the remaining bipartite graph onto the reaction nodes to form a
    monopartite reaction-reaction graph where two reactions are connected if
    there is a two-step path between them (i.e. they share at least one
    metabolite and one produces it and one consumes it or at least one reaction
    is reversible) or they're both upstream or downstream of another reaction
    that wasn't in `exclude` with which they both share at least one metabolite
    that wasn't in `exclude`
    '''
    bip_graph = make_bipartite_graph(model)
    bip_graph.remove_nodes_from(exclude)
    # remove any nodes that are isolated after that
    bip_graph.remove_nodes_from([
        n for n in bip_graph if bip_graph.degree(n) == 0
    ])
    # project onto reaction nodes to get a monopartite reaction-reaction graph
    rxn_ids = [r.id for r in model.reactions]
    rxn_graph = nx.algorithms.bipartite.projected_graph(
        bip_graph, [n for n in bip_graph if n in rxn_ids]
    )
    # remove any reactions that are isolated after that
    rxn_graph.remove_nodes_from([
        r for r in rxn_graph if rxn_graph.degree(r) == 0
    ])
    # make sure we have some edges before trying to make an edge list
    if len(rxn_graph) == 0:
        msg = 'Did not connect any reactions into pathways; maybe try a higher '
        msg += 'degree_cutoff and verify that the given dataframe of test '
        msg += 'results has a decent number of reactions flagged by the '
        msg += 'dead-end, reversibility, and/or dilution tests.'
        raise ValueError(msg)
    # make the graph undirected cuz we don't need two edges between each
    # reversible reaction and all of its metabolites in any visualizations we
    # make of the pathways
    rxn_graph = rxn_graph.to_undirected()
    elist = nx.to_pandas_edgelist(rxn_graph)
    # add the reaction equations to the reaction IDs to make the node labels for
    # the edge list
    elist = add_reaction_equations(
        elist, model, id_col = 'source', use_names = use_names,
        add_suffixes = add_suffixes
    ).rename(columns = {'reaction_equation' : 'r1'})
    elist = add_reaction_equations(
        elist, model, id_col = 'target', use_names = use_names,
        add_suffixes = add_suffixes
    ).rename(columns = {'reaction_equation' : 'r2'})
    elist['reaction_1'] = elist[['source', 'r1']].agg(': '.join, axis = 1)
    elist['reaction_2'] = elist[['target', 'r2']].agg(': '.join, axis = 1)
    elist = elist[['reaction_1', 'reaction_2']]
    # make a dict of reaction IDs to pathway indices (use "pathway" 0 for all
    # reactions that didn't wind up in a pathway)
    pathway_dict = {r.id : 0 for r in model.reactions}
    i = 0
    for pathway in nx.connected_components(rxn_graph):
        i += 1
        for reaction in pathway:
            pathway_dict[reaction] = i
    return((pathway_dict, elist))

def add_reaction_btwns(model, all_test_results, threads = 1):
    '''
    Compute the betweenness for each reaction in the given Cobrapy Model in a
    reaction-reaction graph derived from the model that has one node for each
    reaction that is connected to all other reactions it shares any metabolites
    with. Edges are weighted by 1/the sum of 1/the number of reactions each
    shared metabolite participates in such that reactions that only share common
    metabolites (e.g. ATP, water) have large edge weights and reactions that
    share uncommon metabolites (e.g. intermediates in linear pathways) have
    small edge weights.
    The idea is that you can rank reactions by these betweennesses to help
    prioritize reactions when doing manual curation.
    I highly recommend using more than one thread; it can take multiple days
    with larger (~10,000 reactions) models on a single thread.
    Returns a dict with reaction IDs as keys and betweennesses as values.
    '''
    bip_graph = make_bipartite_graph(model)
    rxn_graph = make_weighted_rxn_graph(bip_graph, model)
    rxn_btwns = parallel_betweenness(rxn_graph, threads)
    # add columns to the dataframe of all test results for the betweenness of
    # each reaction and the largest betweenness of all reactions in each pathway
    all_test_results[
        'reaction_betweenness'
    ] = all_test_results['reaction_id'].map(rxn_btwns).apply(
        # round to 3 significant figures so that rounding errors don't cause
        # weirdness when sorting reactions by their betweennesses
        sigfig_round, sfs = 3
    )
    # in theory, we could compute the betweenness for each pathway by computing
    # the proportion of all shortest paths that pass through the subgraph of the
    # whole network that each pathway constitutes, but that'd be a lot more work
    # and it's not clear to me that it'd significantly change the ordering of
    # pathways when sorted by betweenness, which is all we're using it for
    all_test_results['pathway_betweenness'] = all_test_results.groupby(
        'pathway'
    )['reaction_betweenness'].transform(max)
    all_test_results = all_test_results.sort_values([
        # in case multiple pathways have the same betweenness
        'pathway_betweenness', 'pathway'
    ])
    return(all_test_results)

def make_bipartite_graph(model):
    '''
    Make a directed bipartite graph from a Cobrapy Model object where nodes are
    either reactions or metabolites and edges connect metabolites to the
    reactions they participate in and the direction indicates consumption or
    production (two edges connect each reversible reaction to each metabolite
    that participates in it, one in each direction)
    '''
    graph = nx.DiGraph()
    for reaction in model.reactions:
        # loop over products and reactants separately to direct edges properly
        for met in reaction.reactants:
            graph.add_edge(met.id, reaction.id)
            # if this is a reversible reaction, also add another edge in the
            # opposite direction
            if reaction.reversibility:
                graph.add_edge(reaction.id, met.id)
        for met in reaction.products:
            graph.add_edge(reaction.id, met.id)
            if reaction.reversibility:
                graph.add_edge(met.id, reaction.id)
    return(graph)

def make_weighted_rxn_graph(bip_graph, model):
    '''
    Given a bipartite graph representing a metabolic network, create a
    monopartite graph with one node per reaction in which two reactions are
    connected if they share any metabolites. Set the weights of these edges by
    getting the degree of each shared metabolite in the bipartite network, then
    doing 1/the sum of 1/each metabolite degree, so that edge weights represent
    "distance" between each pair of reactions, and reactions that share highly-
    connected metabolites are "further" from each other than reactions that
    share metabolites that participate in few reactions.
    '''
    rxn_ids = [r.id for r in model.reactions]
    rxn_graph = nx.algorithms.bipartite.projected_graph(
        bip_graph, [n for n in bip_graph if n in rxn_ids]
    )
    # set edge weights by adding the sum of 1/the degree of each shared
    # metabolite in the bipartite graph (also remove edges between pairs of
    # reactions that don't actually share metabolites that networkx added for
    # no discernible reason)
    to_remove = list()
    for (r1, r2) in rxn_graph.edges():
        # figure out which metabolites were shared
        mets1 = bip_graph.neighbors(r1)
        mets2 = bip_graph.neighbors(r2)
        shared = set(mets1).intersection(set(mets2))
        # remove this edge if no metabolites were shared
        if len(shared) == 0:
            to_remove.append((r1, r2))
        else:
            # get 1/degree for each shared metabolite in the bipartite graph
            weight = sum([1/bip_graph.degree(m) for m in shared])
            # use 1/that sum as the edge weight so that they can be interpreted
            # as distances between reactions and reactions only sharing common
            # metabolites are more distant than ones sharing rare ones
            rxn_graph[r1][r2]['weight'] = 1/weight
    # drop the edges that shouldn't have existed in the first place
    rxn_graph.remove_edges_from(to_remove)
    return(rxn_graph)

def parallel_betweenness(G, threads = None):
    '''
    Parallel betweenness centrality function
    (From NetworkX example code, lightly edited)
    '''
    p = Pool(processes = threads)
    node_chunks = list(chunks(G.nodes(), G.order() // (len(p._pool) * 4)))
    num_chunks = len(node_chunks)
    # compute betweennesses using shortest paths between all pairs of many
    # disjoint subsets of nodes in the network that collectively span all nodes
    partial_btwn_lists = p.starmap(
        nx.betweenness_centrality_subset,
        zip(
            [G] * num_chunks,
            node_chunks,
            [list(G)] * num_chunks,
            [True] * num_chunks,
            ['weight'] * num_chunks,
        ),
    )
    # concatenate all of the partial betweennesses obtained with each subset
    full_btwns = partial_btwn_lists[0]
    for some_partial_btwns in partial_btwn_lists[1:]:
        for node in some_partial_btwns:
            full_btwns[node] += some_partial_btwns[node]
    return(full_btwns)

def chunks(l, n):
    '''
    Divide a list of nodes `l` into `n` chunks
    (From NetworkX example code)
    '''
    l_c = iter(l)
    while 1:
        x = tuple(it.islice(l_c, n))
        if not x:
            return
        yield x
