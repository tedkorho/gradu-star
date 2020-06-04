import argparse as ap

def get_args():
        
    """
    Parses the  arguments
    """    

    parser = ap.ArgumentParser( \
        description="Find flares from a lightcurve FITS file")
    parser.add_argument("inputfile",
        metavar="if",type=str,help="Input file")
    parser.add_argument("-o",metavar="of",type=str,help="Output file")
    parser.add_argument("--sens",type=float,help="Sensitivity ()")
    parser.add_argument("--stype",type=str,help="")
    parser.add_argument("-p",type=float,
                            help="Periodicity of the signal in days")
    
    args = parser.parse_args()
    return args

def write_flares(filename, data):
    
    

def lcplot():
    
