#!/usr/bin/env python3

def flexrfindmain(filename,sigmathreshold,plot,peakheight,peakprominence,peakwidth,peakdistance,step,mode,geotolerance,build_limit,ligand,pdb_in,distance,singleconfs):
    from src.tools import ringer_tools
    from src.tools.ringer_tools import ringer_parser, match_and_build,\
                                       parse_peak_find,assemble_matches,\
                                       output_build_list
    from src.flexrpkg import top_level
    from src.flexrpkg.top_level import check_library, test_input_file

    library = check_library()

    test_input_file(filename)

    df = parse_peak_find(filename,sigmathreshold,plot,peakheight,peakprominence,peakwidth,peakdistance,step,mode)
    alt_confs = assemble_matches(df,library,geotolerance)
    build_list = output_build_list(filename,sigmathreshold,alt_confs,build_limit,ligand,pdb_in,distance,singleconfs)
    print('')
    print('Peak finding results ----> peak_finder_output...csv')
    print('Alts for building -------> alts.csv')
    print('')
    #print('Done')

    return build_list
