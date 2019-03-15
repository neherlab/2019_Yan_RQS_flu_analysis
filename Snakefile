include: "Snakefile_base"

segments = ['ha'] #, 'na']
lineages = ['h3n2'] #, 'h1n1pdm', 'vic', 'yam', 'Ball']
resolutions = ['all']

passages = ['cell']
centers = ['cdc']
assays = ['hi']

rule all_lineages:
    input:
        auspice_tree = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tree.json",
                              lineage=lineages, segment=segments, resolution=resolutions),
        auspice_meta = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_meta.json",
                              lineage=lineages, segment=segments, resolution=resolutions),
        auspice_tip_frequencies = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json",
                              lineage=lineages, segment=segments, resolution=resolutions)

rule split:
    input:
        tree = "results/tree_cdc_B_{segment}_all_cell_hi.nwk",
    output:
        yam = "results/tree_cdc_yam_{segment}_all_cell_hi.nwk",
        vic = "results/tree_cdc_vic_{segment}_all_cell_hi.nwk"
    script:
        "scripts/split_trees.py"


def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model,
        rules.traits.output.node_data,
    ]

    # Only request a distance file for builds that have mask configurations
    # defined.
    # if _get_build_mask_config(wildcards) is not None:
    #     inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

# for lineage in ['h1n1', 'h1n1pdm', 'h3n2', 'B']:
#     rule:
#         input:
#             tree = rules.refine.output.tree,
#             metadata = rules.parse.output.metadata,
#             auspice_config = auspice_config,
#             node_data = _get_node_data_for_export,
#             colors=colors
#         output:
#             auspice_tree = expand("auspice/flu_{{center}}_{lineage}_{{segment}}_{{resolution}}_{{passage}}_{{assay}}_tree.json", lineage=[lineage]),
#             auspice_meta = expand("auspice/flu_{{center}}_{lineage}_{{segment}}_{{resolution}}_{{passage}}_{{assay}}_meta.json", lineage=[lineage])
#         shell:
#             """
#             augur export \
#                 --tree {input.tree} \
#                 --metadata {input.metadata} \
#                 --colors {input.colors} \
#                 --node-data {input.node_data} \
#                 --auspice-config {input.auspice_config} \
#                 --output-tree {output.auspice_tree} \
#                 --output-meta {output.auspice_meta}
#             """

rule:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = auspice_config,
        node_data = _get_node_data_for_export,
        colors=colors
    output:
        auspice_tree = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_tree.json",
        auspice_meta = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """

rule export_Bsub:
    input:
        tree_vic = rules.split.output.vic,
        tree_yam = rules.split.output.yam,
        node_data = "results/branch-lengths_cdc_B_{segment}_all_cell_hi.json",
        metadata = "results/metadata_B_{segment}.tsv",
        aa_muts = "results/aa-muts_cdc_B_{segment}_all_cell_hi.json",
        config = "config/auspice_config_B.json",
    output:
        tree_vic_json = "auspice/flu_seasonal_vic_{segment}_all_tree.json",
        meta_vic_json = "auspice/flu_seasonal_vic_{segment}_all_meta.json",
        tree_yam_json = "auspice/flu_seasonal_yam_{segment}_all_tree.json",
        meta_yam_json = "auspice/flu_seasonal_yam_{segment}_all_meta.json",
        node_data_yam = "results/branch-lengths_cdc_yam_{segment}_all_cell_hi.json",
        node_data_vic = "results/branch-lengths_cdc_vic_{segment}_all_cell_hi.json"
    shell:
        """
        augur export \
            --tree {input.tree_vic} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.config} \
            --output-tree {output.tree_vic_json} \
            --output-meta {output.meta_vic_json} &\
        augur export \
            --tree {input.tree_yam} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.config} \
            --output-tree {output.tree_yam_json} \
            --output-meta {output.meta_yam_json} &\
        cp {input.node_data} {output.node_data_vic} &\
        cp {input.node_data} {output.node_data_yam}
        """

rule simplify_auspice_names:
    input:
        tree = "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi_tree.json",
        meta = "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi_meta.json",
        frequencies = "auspice/flu_cdc_{lineage}_{segment}_{resolution}_cell_hi_tip-frequencies.json"
    output:
        tree = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tree.json",
        meta = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_meta.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_tip-frequencies.json"
    shell:
        '''
        cp {input.tree} {output.tree} &
        cp {input.meta} {output.meta} &
        cp {input.frequencies} {output.frequencies} &
        '''

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice ",
    shell:
        "rm -rfv {params}"


rule TMRCA:
    input:
        tree = "results/tree_cdc_{lineage}_{segment}_all_cell_hi.nwk",
        node_data = "results/branch-lengths_cdc_{lineage}_{segment}_all_cell_hi.json",
    output:
        "analysis/{lineage}_{segment}_Tmrca_trajectory.dat"
    shell:
        """
        python scripts/plot_Tmrca.py --tree {input.tree} --node-data {input.node_data} --output {output}
        """

rule allTMRCA:
    input:
        expand("analysis/{lineage}_{segment}_Tmrca_trajectory.dat", segment=['ha'],
                lineage=['h1n1', 'h1n1pdm', 'h3n2', 'vic', 'yam'])
    output:
        'figures/all_Tmrca.pdf'
    shell:
        "python scripts/plot_all_Tmrca.py"
