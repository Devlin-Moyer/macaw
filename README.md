# Metabolic Accuracy Checks And Workflow (MACAW)

![MACAW logo](logo.png)

MACAW is a collection of tests for automatically highlighting reactions in an arbitrary Genome-Scale Metabolic Model (GSMM) that are likely to contain or be near an error of some kind. It also provides some information that suggests what changes would be necessary to correct the error. MACAW should work with any GSMM that can be parsed by [Cobrapy](https://opencobra.github.io/cobrapy/) (which currently is currently capable of reading GSMMs saved as SBML, Matlab, JSON, or YAML files).

If you use MACAW, cite [this paper](https://doi.org/10.1186/s13059-025-03533-6).

## Quick Start

0. Requires Python 3.8 or newer

1. Install with pip: `pip install git+https://github.com/Devlin-Moyer/macaw.git@main`

2. Read in a GSMM and run the tests:

```python
import cobra
from macaw.main import run_all_tests

model = cobra.io.read_sbml_model('GSMMs/iML1515.xml')
(test_results, edge_list) = run_all_tests(model)
```

`test_results` will be a `Pandas.DataFrame` with one row for each reaction in the model and columns containing the results of each test (described in detail below).
`edge_list` will be a list of lists where each sub-list has two reaction IDs; you can use it to visualize a network connecting all reactions flagged by any of the tests (with e.g. Cytoscape, Networkx, etc.)

All of the code in `figure_scripts` is only necessary to reproduce the figures used in the paper; the only code that actually constitutes MACAW is in `macaw`.

## The Tests

<details>
  <summary><code>dead_end_test</code></summary>

  Looks for metabolites in the given GSMM that can only be produced by all reactions they participate in or only consumed, then identifies all reactions that are prevented from sustaining steady-state fluxes because of each of these dead-end metabolites. The simplest case of a dead-end metabolite is one that only participates in a single reaction. Also flags all reversible reactions that can only carry fluxes in a single direction because one of their metabolites can either only be consumed or only be produced by all other reactions it participates in.

  Arguments:

  - `given_model`: the Cobrapy Model object containing the GSMM to be tested.
  - `use_names`: (optional) whether or not to use the "name" attributes of the metabolites in `given_model` in the reaction equation column in the output DataFrame instead of the metabolite IDs. False by default.
  - `add_suffixes`: (optional) whether or not to add suffixes indicating which compartment each metabolite is in to the names or IDs of metabolites in the reaction equation column in the output. Setting `add_suffixes` to True and `use_names` to False is generally not recommended, as most GSMMs with multiple compartments already encode the compartment each metabolite is in in the metabolite's ID. False by default.
  - `verbose`: (optional) controls how many messages are printed when the test runs. Set to 0 to print no messages. 1 by default.

  Returns an edge list defining a network that connects each dead-end metabolite to all the reactions it blocks fluxes through and a `Pandas.DataFrame` with one row for each reaction in the given GSMM, columns for the IDs and equations of each reaction, and the result of the dead-end test:

  - "ok" if the reaction was not a dead-end.
  - one or more metabolite IDs separated by semicolons indicating which dead-end metabolites participate in that reaction.
  - "only when going forwards" or "only when going backwards" if it was a reversible reaction that was not a dead-end but had a reactant or product that could only be consumed by all other reactions it participates in or only be produced by all other reactions it participates in.

</details>

<details>
  <summary><code>dilution_test</code></summary>

  Separately tests each metabolite in the given GSMM to see if adding a dilution reaction and dilution constraint for that metabolite renders all reactions that it participates in incapable of non-zero steady-state fluxes. A dilution reaction just consumes a single metabolite and produces nothing, and dilution constraint sets the flux through a particular metabolite's dilution reaction equal to some fraction of the sum of the absolute values of the fluxes through all other reactions that that metabolite participates in. Dilution constraints generally only block fluxes through metabolites that can only be recycled within a GSMM and lack a biosynthesis or uptake pathway. If you run the dilution test multiple times on the same GSMM, it will sometimes flag 1-3 more or fewer reactions on different runs. I haven't figured out why.

  Arguments:

  - `given_model`: the Cobrapy Model object containing the GSMM to be tested.
  - `dead_end_results`: (optional) the `Pandas.DataFrame` returned by `dead_end_test` on `given_model`. Will not verify that the DataFrame contains results from running `dead_end_test` on `given_model`, so providing results from other GSMMs may produce unusual errors. If not provided, will run the dead-end test automatically before beginning the dilution test.
  - `media_mets`: (optional) list of IDs of metabolites in `given_model` that you want to allow uptake of through exchange reactions. If empty, will not alter bounds on exchange reactions. Otherwise, will set the lower bounds on all exchange reactions except those involving metabolites in `media_mets` to 0 to prevent their uptake. This can significantly increase the number of reactions flagged by the dilution test if there are exchange reactions for metabolites that the cell being modeled should be capable of producing on its own and is unlikely to encounter in its surroundings. Consider using the list of ingredients in a defined culture medium for the cell(s) the GSMM represents, if one exists.
  - `zero_thresh`: (optional) how close to zero is close enough to consider a reaction incapable of sustaining flux? 10^-8 by default.
  - `timeout`: (optional) sometimes, Cobrapy/the underlying linear programming optimizer will hang when optimizing models with dilution constraints. How long should the script wait on results for a single metabolite before giving up and starting over? Keep in mind that it takes much longer to test every single reaction that e.g. water or ATP participates in than most other metabolites in most GSMMs. 1800 seconds (30 minutes) by default.
  - `max_attempts`: (optional) if it takes longer than `timeout` to test a single metabolite, how many total times should that metabolite be tested before giving up and assuming that it is probably dilution-blocked? 3 by default.
  - `use_names`: (optional) whether or not to use the "name" attributes of the metabolites in `model` in the reaction equation column in the output DataFrame instead of the metabolite IDs. False by default.
  - `add_suffixes`: (optional) whether or not to add suffixes indicating which compartment each metabolite is in to the names or IDs of metabolites in the reaction equation column in the output. Setting `add_suffixes` to True and `use_names` to False is generally not recommended, as most GSMMs with multiple compartments already encode the compartment each metabolite is in in the metabolite's ID. False by default.
  - `verbose`: (optional) controls how many messages are printed when the test runs. Set to 0 to print no messages. 1 by default.

  Returns an edge list defining a network that connects each dilution-blocked metabolite to all the reactions its dilution constraint blocks and a `Pandas.DataFrame` with one row for each reaction in the given GSMM, columns for the IDs and equations of each reaction, and the result of the dilution test:

  - "ok" if the reaction was always capable of non-zero fluxes when any individual metabolite's dilution constraint was imposed on the model.
  - "always blocked" if the reaction was incapable of non-zero fluxes regardless of whether or not any dilution constraints were imposed on the model.
  - "blocked by dilution" if the reaction was capable of non-zero fluxes when no dilution constraints were imposed on the model but became incapable of non-zero fluxes when one or more metabolites' dilution constraints were imposed.
  - "unblocked by dilution" if the reaction was only capable of non-zero fluxes when at least one metabolite's dilution constraint was imposed on the model. This is rare, and ideally all such reactions would be flagged by the dead-end test and blocked in the pre-processing step of the dilution test where it sets both bounds of all reactions flagged by the dead-end test to 0.

</details>

<details>
  <summary><code>diphosphate_test</code></summary>

  Identifes all reversible reactions that involve diphosphate that aren't transporting it between compartments. Requires the IDs of the metabolites in the GSMM that represent the diphosphate and inorganic (mono)phosphate ions. Most reactions involving the diphosphate ion should be irreversible in the direction that produces diphosphate, since most cells express a variety of highly active diphosphatases that quickly turn diphosphate ions into two separate inorganic phosphate ions. While most other reactions involving diphosphate, specifically those that involve separating a (d)NTP into a (d)NMP + a diphosphate, have Gibbs free energy changes of approximately zero and are thus readily reversible, the diphosphate reaction is quite exergonic, so the diphosphatases present in most cells generally drive all other diphosphate-producing reactions in the direction of diphosphate production. Leaving these reactions as reversible when predicting steady-state fluxes from a GSMM can result in unrealistic predictions about how ATP is synthesized and create loops involving chains of reversible diphosphate reactions. This test generally flags very very few reactions in most GSMMs.

  Arguments:

  - `given_model`: the Cobrapy Model object containing the GSMM to be tested.
  - `ppi_ids`: the IDs of metabolites in `given_model` that represent diphosphate ions.
  - `pi_ids`: the IDs of metabolites in `given_model` that represent inorganic (mono)phosphate ions.
  - `use_names`: (optional) whether or not to use the "name" attributes of the metabolites in `model` in the reaction equation column in the output DataFrame instead of the metabolite IDs. False by default.
  - `add_suffixes`: (optional) whether or not to add suffixes indicating which compartment each metabolite is in to the names or IDs of metabolites in the reaction equation column in the output. Setting `add_suffixes` to True and `use_names` to False is generally not recommended, as most GSMMs with multiple compartments already encode the compartment each metabolite is in in the metabolite's ID. False by default.
  - `verbose`: (optional) controls how many messages are printed when the test runs. Set to 0 to print no messages. 1 by default.

  Returns a `Pandas.DataFrame` with one row for each reaction in the given GSMM, columns for the IDs and equations of each reaction, and the result of the diphosphate test:

  - "ok" if the reaction either does not involve diphosphate or is already irreversible.
  - "should be irreversible" if the reaction is reversible and diphosphate is a product.
  - "should be flipped and made irreversible" if the reaction is reversible and diphosphate is a reactant (the suggestion to flip such reactions is to ensure that no reactions in the GSMM are only capable of sustaining non-positive fluxes, which won't break any math one might want to do with a GSMM but might be confusing or aesthetically unappealing).

  Unlike the other tests, the diphosphate test does not return an edge list connecting the reactions it flags.

</details>

<details>
  <summary><code>duplicate_test</code></summary>

  Identifies sets of reactions that may be duplicates of each other because they:
  - Involve exactly the same metabolites with exactly the same stoichiometric coefficients (but potentially different associated genes).
  - Involve exactly the same metabolites, but go in different directions and/or some are reversible and some are not.
  - Involve exactly the same metabolites, but with different stoichiometric coefficients.
  - Represent the oxidation and/or reduction of the same metabolite, but use different electron acceptors/donors from the given list of pairs of oxidized and reduced forms of various electron carriers (e.g. NAD(H), NADP(H), FAD(H2), ubiquinone/ubiquinol, cytochromes).

  It is possible for a single reaction to fit in multiple of the above categories. There are sometimes cases where sets of reactions that fall into one of the above categories are completely legitimate representations of real biochemistry (e.g. separate irreversible reactions for importing vs exporting the same metabolite because two different transporters encoded by different genes are each responsible for transporting that metabolite in only one direction, enzymes that can use NAD(H) or NADP(H) interchangeably to catalyze the same redox reaction), but reactions that meet these criteria are generally worth close examination to ensure that they should actually all exist as separate reactions.

  Arguments:

  - `model`: the Cobrapy Model object containing the GSMM to be tested.
  - `redox_pairs`: (optional) a list of lists or tuples that each have exactly two strings corresponding to the IDs of metabolites in `model` that represent the oxidized and reduced forms of the same metabolite. For example, in a model that uses BiGG IDs for all metabolites, `redox_pairs` might look like `[('nad_c', 'nadh_c'), ('nadp_c', 'nadph_c')]`. Providing more pairs of IDs will generally lead to more reactions being flagged as redox duplicates, but the test does nothing to ensure that the provided pairs of metabolites actually represent oxidized and reduced forms of the same metabolite. Ignored if `proton_ids` is not also provided.
  - `proton_ids`: (optional) a list of strings containins IDs of metabolites in `model` that represent protons. For example, in a model that uses BiGG IDs for all metabolites, `proton_ids` might look like `['h_c', 'h_p', 'h_e']`. Ignored if `redox_pairs` is not also provided.
  - `use_names`: (optional) whether or not to use the "name" attributes of the metabolites in `model` in the reaction equation column in the output DataFrame instead of the metabolite IDs. False by default.
  - `add_suffixes`: (optional) whether or not to add suffixes indicating which compartment each metabolite is in to the names or IDs of metabolites in the reaction equation column in the output DataFrame. Setting `add_suffixes` to True and `use_names` to False is generally not recommended, as most GSMMs with multiple compartments already encode the compartment each metabolite is in in the metabolite's ID. False by default.
  - `verbose`: (optional) controls how many messages are printed when the test runs. Set to 0 to print no messages. 1 by default.

  Returns an edge list describing a network with one node for each reaction flagged as a potential duplicate where reactions are connected to the other reactions that they are potentially duplciates of, as well as a `Pandas.DataFrame` with one row for each reaction in the given GSMM, columns for the IDs and equations of each reaction, and several columns indicating the results of the duplicate test:

  - `duplicate_test_exact`: "ok" if the reaction had no exact duplicates or a semicolon-delimited list of the IDs of other reactions that were exact duplicates.
  - `duplicate_test_directions`: "ok" if there were no other reactions that involved the same metabolites but went in the opposite direction or had the opposite reversibility or a semicolon-delimited list of the IDs of those other reactions.
  - `duplicate_test_coefficients`: "ok" if there were no other reactions that involved the same metabolites but with different stoichiometric coefficients or a semicolon-delimited list of the IDs of those other reactions.
  - `duplicate_test_redox`: "ok" if there were no other reactions that involved the same metabolites aside from the ones provided in `redox_pairs` and `proton_ids` or a semicolon-delimited list of the IDs of those other reactions. N/A if `redox_pairs` or `proton_ids` were not provided.

</details>

<details>
  <summary><code>loop_test</code></summary>

  Identifies all reactions that are capable of sustaining non-zero fluxes when all exchange reactions (i.e. reactions representing the uptake and/or secretion of individual metabolites) are blocked. Also attempts to determine which "loop" each such reaction is a member of by generating 1,000 possible solutions to the GSMM, getting pairwise correlations between the distributions of 1,000 possible fluxes for each reaction, and identifying groups of reactions whose fluxes were highly correlated. Removes any objective functions from the given GSMM and sets all non-zero lower bounds (e.g. lower bounds on ATP maintenance reactions) to zero before starting.

  Arguments:

  - `model`: the Cobrapy Model object containing the GSMM to be tested.
  - `zero_thresh`: (optional) how close to zero is close enough to consider a reaction incapable of sustaining flux? 10^-8 by default.
  - `corr_thresh`: (optional) how correlated do the distributions of possible fluxes for two reactions have to be in order to consider them members of the same loop? Default is 0.9, which corresponds to correlations above 0.9 or below -0.9.
  - `use_names`: (optional) whether or not to use the "name" attributes of the metabolites in `model` in the reaction equation column in the output DataFrame instead of the metabolite IDs. False by default.
  - `add_suffixes`: (optional) whether or not to add suffixes indicating which compartment each metabolite is in to the names or IDs of metabolites in the reaction equation column in the output. Setting `add_suffixes` to True and `use_names` to False is generally not recommended, as most GSMMs with multiple compartments already encode the compartment each metabolite is in in the metabolite's ID. False by default.
  - `verbose`: (optional) controls how many messages are printed when the test runs. Set to 0 to print no messages. 1 by default.

  Returns an edge list defining a network connecting reactions that have at least one metabolite in common and have highly-correlated, as well as a `Pandas.DataFrame` with one row for each reaction in the given GSMM, columns for the IDs and equations of each reaction, and the result of the loop test:

  - "ok" for reactions that were not capable of non-zero fluxes when all exchange reactions were blocked.
  - "in loop" for reactions that were capable of non-zero fluxes when all exchange reactions were blocked.

</details>

## Other Useful Functions

<details>
  <summary><code>form_pathways</code></summary>

  Combines the edge lists produced by multiple of the above tests into a single comprehensive network. This is non-trivial because the dead-end and dilution tests produce edge lists that describe bipartite networks in which some nodes represent reactions and others represent metabolites, while the duplicate and loop tests produce edge lists that describe monopartite networks in which all nodes represent reactions. The resulting network generally contains many connected components. `form_pathways` will assign a unique integer to each component and add a column to the Pandas Dataframe of results from all tests indicating which connected component each reaction is in. Reactions that were not flagged by any tests or not connected to any other reactions that were flagged by any tests (this only happens with reactions flagged by the dead-end or diphosphate tests, and is generally uncommon) are always assigned a "pathway" of 0.

</details>

<details>
  <summary><code>run_all_tests</code></summary>

  Runs all four tests on the given model and calls `form_pathways` to combine the edge lists into one.

</details>

<details>
  <summary><code>simplify_test_results</code></summary>

  Makes each column in the `Pandas.DataFrame` produced by any test just say "ok" or "bad" for each reaction (most tests have more complicated/specific/variable information for each reaction in their column). Also merges the 4 duplicate test columns into a single column.

</details>
