from argparse import ArgumentParser
import sys

from bioseq_pipeline.organisms.sars2.pipeline import (
    add_sars2_arguments,
    parse_genbank_to_excel,
    run_from_args as run_sars2_from_args,
)


def run_tuberculosis_placeholder(_args):
    raise SystemExit(
        "Error: tuberculosis pipeline is not implemented yet. "
        "The directory is reserved for the future TB-specific workflow."
    )


def parse_args(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or (argv[0].startswith("-") and argv[0] not in ("-h", "--help")):
        argv.insert(0, "sars2")

    parser = ArgumentParser(description="Run organism-specific bioinformatics pipelines.")
    subparsers = parser.add_subparsers(dest="organism", required=True)

    sars2_parser = subparsers.add_parser(
        "sars2",
        aliases=["sars-cov-2", "covid"],
        help="Run the SARS-CoV-2 pipeline.",
    )
    add_sars2_arguments(sars2_parser)
    sars2_parser.set_defaults(handler=run_sars2_from_args)

    tb_parser = subparsers.add_parser(
        "tuberculosis",
        aliases=["tb"],
        help="Reserved entry point for a future tuberculosis pipeline.",
    )
    tb_parser.set_defaults(handler=run_tuberculosis_placeholder)

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        args.handler(args)
    except (PermissionError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}") from None
