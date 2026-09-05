#!/usr/bin/env python3
"""Generate the nf-core/chipseq MultiQC Workflow Summary custom content.

Port of the nf-core template paramsSummaryMap (nfc_schema.groovy) +
paramsSummaryMultiqc (funcs.txt) pair, rendered against the oxo-flow engine:
per schema group, keep a param when its value differs from the schema default
(or, without a default, when it is non-empty and not false), then emit the
HTML <dl> section as a MultiQC custom-content YAML (workflow_summary_mqc.yaml,
auto-discovered by the `multiqc -f ... .` tree scan from the run root).
"""
import argparse
import os

N_A_SPAN = '<span style="color:#999999;">N/A</a>'


def parse_bool(text):
    return text == "true"


def groovy_to_string(value):
    """Mirror Groovy GString coercion for the values we render."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


# TEMPLATE_PARAMS TAIL (appended below)

# Schema groups from nf-core/chipseq 2.1.0 nextflow_schema.json
# (params_summary_map compares against these defaults). Only groups/params the
# port exposes in its [config] are listed; port-only keys live in "Port options".
SCHEMA_GROUPS = [
    ("Input/output options", [
        ("fragment_size", "integer", 200),
        ("seq_center", "string", None),
        ("read_length", "integer", None),
        ("multiqc_title", "string", None),
    ]),
    ("Reference genome options", [
        ("fasta", "string", None),
        ("gtf", "string", None),
        ("bwa_index", "string", None),
        ("bowtie2_index", "string", None),
        ("chromap_index", "string", None),
        ("star_index", "string", None),
        ("gene_bed", "string", None),
        ("macs_gsize", "number", None),
        ("blacklist", "string", None),
    ]),
    ("Adapter trimming options", [
        ("clip_r1", "integer", None),
        ("clip_r2", "integer", None),
        ("three_prime_clip_r1", "integer", None),
        ("three_prime_clip_r2", "integer", None),
        ("trim_nextseq", "integer", None),
        ("skip_trimming", "boolean", None),
    ]),
    ("Alignment options", [
        ("aligner", "string", "bwa"),
        ("keep_dups", "boolean", None),
        ("keep_multi_map", "boolean", None),
        ("bwa_min_score", "integer", None),
    ]),
    ("Peak calling options", [
        ("narrow_peak", "boolean", None),
        ("broad_cutoff", "number", 0.1),
        ("macs_fdr", "number", None),
        ("macs_pvalue", "number", None),
        ("min_reps_consensus", "integer", 1),
        ("save_macs_pileup", "boolean", None),
        ("skip_peak_qc", "boolean", None),
        ("skip_peak_annotation", "boolean", None),
        ("skip_consensus_peaks", "boolean", None),
    ]),
    ("Process skipping options", [
        ("skip_fastqc", "boolean", None),
        ("skip_picard_metrics", "boolean", None),
        ("skip_preseq", "boolean", None),
        ("skip_plot_profile", "boolean", None),
        ("skip_plot_fingerprint", "boolean", None),
        ("skip_spp", "boolean", None),
        ("skip_deseq2_qc", "boolean", None),
        ("skip_igv", "boolean", None),
        ("skip_multiqc", "boolean", None),
        ("skip_qc", "boolean", None),
    ]),
    ("Generic options", [
        ("fingerprint_bins", "integer", 500000),
    ]),
    ("Port options", [
        ("library_ids", "string", None),
        ("raw_dir", "string", None),
        ("chrom_sizes", "string", None),
        ("fai", "string", None),
        ("antibody", "string", None),
        ("multi_antibody", "boolean", None),
        ("replicates_exist", "boolean", None),
        ("multiple_groups", "boolean", None),
    ]),
]


def params_summary_map(config, engine_version, launch_dir, command_line):
    """Port of paramsSummaryMap: keep a param when the port config has the key
    AND (no schema default → value non-empty/non-false; else value != default).
    The 'Core Nextflow options' group is adapted to the oxo-flow engine:
    engine version, launch dir, invocation command line (upstream reads
    workflow.revision/runName/containerEngine/... which do not exist here).
    Config values arrive as CLI strings: "true"/"false" parse to bool so the
    false-vs-default comparison matches Groovy semantics."""
    core = {
        "version": engine_version,
        "commandLine": command_line,
        "launchDir": launch_dir,
        "projectDir": os.path.dirname(os.path.abspath(__file__)),
    }
    params_summary = {"Core Nextflow options": core}
    for group_name, group_params in SCHEMA_GROUPS:
        sub_params = {}
        for param, param_type, schema_value in group_params:
            if param not in config:
                continue
            params_value = config[param]
            if isinstance(params_value, str) and param_type == "boolean":
                if params_value in ("true", "false"):
                    params_value = parse_bool(params_value)
            elif param_type in ("integer", "number") and isinstance(params_value, str):
                # CLI numerics arrive as strings; coerce so falsiness and
                # rendering match upstream Groovy (0 → N/A span at render).
                try:
                    params_value = int(params_value)
                except ValueError:
                    try:
                        params_value = float(params_value)
                    except ValueError:
                        pass
            params_text = groovy_to_string(params_value)
            if schema_value is not None:
                schema_text = groovy_to_string(schema_value)
                if params_value != schema_text:
                    sub_params[param] = params_value
            else:
                # No schema default: keep unless Groovy-falsy ""/null/false.
                # Identity checks keep numeric 0 (upstream 0 != false); the
                # render-time ?: quirk then shows 0 as the N/A span.
                if params_value != "" and params_value is not None and params_value is not False:
                    sub_params[param] = params_value
        # Groovy renders booleans as lowercase true/false in GStrings.
        for param, value in list(sub_params.items()):
            if isinstance(value, bool):
                sub_params[param] = groovy_to_string(value)
        params_summary[group_name] = sub_params
    return params_summary


def workflow_summary_yaml(summary_params):
    """Port of paramsSummaryMultiqc: HTML dl rows per group, N/A span for
    Groovy-falsy values, params sorted alphabetically within each group."""
    summary_section = ""
    for group, group_params in summary_params.items():
        if not group_params:
            continue
        summary_section += f'    <p style="font-size:110%"><b>{group}</b></p>\n'
        summary_section += '    <dl class="dl-horizontal">\n'
        for param in sorted(group_params):
            rendered = group_params[param]
            if rendered is False:
                rendered = "false"
            elif not rendered:
                rendered = N_A_SPAN
            summary_section += f"        <dt>{param}</dt><dd><samp>{rendered}</samp></dd>\n"
        summary_section += "    </dl>\n"
    yaml_text = "id: 'nf-core-chipseq-summary'\n"
    yaml_text += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_text += "section_name: 'nf-core/chipseq Workflow Summary'\n"
    yaml_text += "section_href: 'https://github.com/nf-core/chipseq'\n"
    yaml_text += "plot_type: 'html'\n"
    yaml_text += "data: |\n"
    yaml_text += summary_section
    return yaml_text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--engine-version", required=True, help="oxo-flow engine version (upstream: workflow.nextflow.version)")
    parser.add_argument("--engine-command", default="oxo-flow run main.oxoflow", help="invocation command line (upstream: workflow.commandLine)")
    parser.add_argument("--config", action="append", default=[], metavar="KEY=VALUE", help="pipeline config key=value pairs for the Workflow Summary section (repeatable)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    config = dict(kv.split("=", 1) for kv in args.config if "=" in kv)
    summary_params = params_summary_map(
        config, args.engine_version, os.path.abspath(os.path.join(args.out_dir, os.pardir)), args.engine_command
    )
    with open(os.path.join(args.out_dir, "workflow_summary_mqc.yaml"), "w") as f:
        f.write(workflow_summary_yaml(summary_params))


if __name__ == "__main__":
    main()
