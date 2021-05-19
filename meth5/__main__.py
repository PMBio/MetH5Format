import argparse
from pathlib import Path

from meth5.meth5 import MetH5File
import meth5.main_create_h5
import meth5.main_annotate_reads

def argtype_M5File(value):
    try:
        MetH5File(value, "r").get_chromosomes()
    except:
        raise argparse.ArgumentTypeError(f"Failed to read '{value}'. Is it a valid MetH5 file?")
    return Path(value)

def main():
    parser = argparse.ArgumentParser(description="MetH5 file tools")
    
    parser.add_argument(
        "--chunk_size",
        type=int,
        required=False,
        default=int(1e6),
        help="Number of llrs per chunk",
    )
    
    subparsers = parser.add_subparsers(description="Subcommand: ", dest="subcommand")
    subparsers.required = True
    
    sc_args = subparsers.add_parser(
        "create_h5",
        description="Create H5 file from Nanopolish result files",
    )
    sc_args.set_defaults(func=meth5.main_create_h5.main)
    
    sc_args.add_argument(
        "--input_dir",
        type=Path,
        required=True,
        help="Input directory containing Nanopolish result files",
    )
    
    sc_args.add_argument(
        "--output_file",
        type=Path,
        required=True,
        help="Output MetH5 file",
    )
    
    sc_args.add_argument(
        "--quiet",
        action="store_true",
        help="No progress bar or warnings will be displayed",
    )
    
    sc_args.add_argument(
        "--compression",
        type=str,
        required=False,
        default="gzip",
        choices=["gzip", "None"],
        help="Compression method for the MetH5 data structures. Use 'gzip' for smaller file size, or 'None' for faster read and write speeds",
    )
    
    sc_args.add_argument(
        "--allowed_chromosomes",
        type=str,
        nargs="+",
        required=False,
        default=None,
        help="Only include these chromosomes",
    )


    # New command: Annotating reads
    sc_args = subparsers.add_parser("annotate_reads", description="Annotate reads with read group")
    sc_args.set_defaults(func=meth5.main_annotate_reads.main)
    
    sc_args.add_argument("--m5file", required=True, type=argtype_M5File,
        help="MetH5 file containing methylation calls", )
    
    sc_args.add_argument("--read_groups_key", type=str, required=True,
        help="Read groups key under which the groups should be stored", )
    
    sc_args.add_argument(
        "--read_group_file",
        type=Path,
        required=True,
        help="Tab-delimited file containing columns read_name and numeric group",
    )
    
    

    
    args = parser.parse_args()
    args_dict = vars(args)
    # Remove arguments that the subcommand doesn't take
    subcommand = args.func
    del args_dict["subcommand"]
    del args_dict["func"]
    
    subcommand(**args_dict)
