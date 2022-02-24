import argparse

import meth5.main_create_m5
import meth5.main_merge_m5
import meth5.main_annotate_reads
import meth5.main_list_chunks
import meth5.main_bedgraph


def main():
    parser = argparse.ArgumentParser(description="MetH5 file tools")
    
    parser.add_argument(
        "--chunk_size",
        type=int,
        required=False,
        default=int(5e4),
        help="Number of llrs per chunk",
    )
    
    # New command: merge_meth5
    subparsers = parser.add_subparsers(description="Subcommand: ", dest="subcommand")
    subparsers.required = True
    
    subcommands = {
        "create_m5": meth5.main_create_m5,
        "merge_m5": meth5.main_merge_m5,
        "annotate_reads": meth5.main_annotate_reads,
        "list_chunks": meth5.main_list_chunks,
        "bedgraph": meth5.main_bedgraph,
    }
    
    for subcommand, command_module in subcommands.items():
        sc_args = subparsers.add_parser(
            subcommand,
            description=command_module.__description__,
        )
        sc_args.set_defaults(func=command_module.main)
        command_module.set_arguments(sc_args)
    
    args = parser.parse_args()
    args_dict = vars(args)
    # Remove arguments that the subcommand doesn't take
    subcommand = args.func
    del args_dict["subcommand"]
    del args_dict["func"]
    
    subcommand(**args_dict)
