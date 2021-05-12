import argparse
import pathlib

import meth5.main_create_h5

def main():
    parser = argparse.ArgumentParser(
        description="MetH5 file tools")
    
    parser.add_argument("--chunk_size", type=int, required=False, default=int(1e6),
        help="Number of llrs per chunk", )
    
    subparsers = parser.add_subparsers(description="Subcommand: ", dest="subcommand")
    subparsers.required = True
    
    sc_args = subparsers.add_parser("create_h5",
        description="Create H5 file from Nanopolish result files", )
    sc_args.set_defaults(func=meth5.main_create_h5.main)
    
    sc_args.add_argument("--input_dir", type=pathlib.Path, required=True,
        help="Input directory containing Nanopolish result files", )

    sc_args.add_argument("--output_file", type=pathlib.Path, required=True,
                         help="Output MetH5 file", )

    sc_args.add_argument("--quiet", action="store_true", help="No progress bar or warnings will be displayed", )

    sc_args.add_argument("--compression", type=str, required=False, default="gzip",
                         choices = ["gzip", "None"],
                         help="Compression method for the MetH5 data structures. Use 'gzip' for smaller file size, or 'None' for faster read and write speeds")
    
    args = parser.parse_args()
    args_dict = vars(args)
    # Remove arguments that the subcommand doesn't take
    subcommand = args.func
    del args_dict["subcommand"]
    del args_dict["func"]
    
    subcommand(**args_dict)
