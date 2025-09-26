import argparse
import argparse

def create_plot_arguments(subparser: argparse._SubParsersAction) -> None:
    """
    This function refactors common arguments into a reusable block.
    """
    plot_args = subparser.add_parser(
        "plot_popgen_stats", 
        help="Plotting the SNP density, Nucleotide diversity, and Tajimas D over chromosome length"
    )

    plot_args.add_argument('-i', '--input_directory', 
                        required=True, 
                        type=str, 
                        help='Input directory containing files of type snpdens, tajimaD, or window.pi.')
    plot_args.add_argument('-p', '--file_pattern', 
                        type=str, 
                        required=True,
                        help='File pattern for the files to plot (e.g., *_SNP_Density_10kb*)')
    plot_args.add_argument('-o', '--output_file', 
                        required=True, 
                        type=str,
                        help='Name of the output plot file (without extension)')
    plot_args.add_argument('-f', '--output_file_format', 
                        required=True, 
                        type=str, 
                        choices=['svg', 'png', 'jpeg'], 
                        help="Format of the output figure (e.g., svg, png, jpeg).")
    plot_args.add_argument('-X', '--X_axis_values', 
                        type=str,
                        required=True,
                        help='Column name from input files to be plotted on the X-axis')
    plot_args.add_argument('-Y', '--Y_axis_values',
                        type=str,
                        required=True,
                        help='Column name from input files to be plotted on the Y-axis')
    plot_args.add_argument('-y', '--y_axis_title',
                        type=str,
                        required=True,
                        help='Title for the Y-axis')
    plot_args.add_argument('-c', '--plot_line_color',
                        required=True,
                        type=str,
                        help='Color of the plot line')

def create_distribution_arguments(subparser: argparse._SubParsersAction) -> None:
    """
    This function defines arguments for the plot_data_dist command.
    """
    dist_args = subparser.add_parser(
        "plot_data_dist", 
        help="Plot the distribution of the data as histograms"
    )

    dist_args.add_argument('-i', '--input_directory', 
                        required=True, 
                        type=str, 
                        help='Input directory containing files of type snpdens, tajimaD, or window.pi.')
    dist_args.add_argument('-p', '--file_pattern', 
                        type=str, 
                        required=True,
                        help='File pattern for the files to plot (e.g., *_SNP_Density_10kb*)')
    dist_args.add_argument('-o', '--output_file', 
                        required=True, 
                        type=str,
                        help='Name of the output plot file (without extension)')
    dist_args.add_argument('-f', '--output_file_format', 
                        required=True, 
                        type=str, 
                        choices=['svg', 'png', 'jpeg'], 
                        help="Format of the output figure (e.g., svg, png, jpeg).")
    dist_args.add_argument('-X', '--X_axis_values', 
                        type=str,
                        required=True,
                        help='Column name from input files to be plotted on the X-axis')
    dist_args.add_argument('-Y', '--Y_axis_values',
                        type=str,
                        required=True,
                        help='Column name from input files to be plotted on the Y-axis')
    dist_args.add_argument('-y', '--y_axis_title',
                        type=str,
                        required=True,
                        help='Title for the Y-axis')
    dist_args.add_argument('-c', '--plot_line_color',
                        required=True,
                        type=str,
                        help='Color of the plot line')

def CommandLineArguments():
    parser = argparse.ArgumentParser(prog="pgsp",
                                     usage='%(prog)s [options]', 
                                     description="Population Genomics Statistics Plotter")
    subparser: argparse._SubParsersAction[argparse.ArgumentParser] = parser.add_subparsers(dest="command", required=True)
    
    # Create arguments for both commands
    create_plot_arguments(subparser)
    create_distribution_arguments(subparser)

    cl_arguments = parser.parse_args()
    return cl_arguments