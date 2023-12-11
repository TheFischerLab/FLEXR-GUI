# -*- coding: utf-8 -*-
#
#    FLEXR:
#    This script detects conformational changes from Ringer measurements and
#    automatically builds them into a multiconformer model for refinement.
#
#    Authors: Tim Stachowski & Marcus Fischer
#    Email: tim.stachowski@stjude.org
#    Copyright 2022 St. Jude Children's Research Hospital, Memphis, TN
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os

def main(ARGS):

    from src.flexrpkg import top_level
    from src.flexrpkg.top_level import intro_messages,args,create_log

    intro_messages()

    #ARGS = args()
    #ARGS = arguments

    create_log(ARGS)

    FS = False

    if ARGS.mode == "DELTA":
        from src.delta import ringer_delta
        from src.delta.ringer_delta import deltamain
        deltamain(ARGS.colors,ARGS.reload,ARGS.chi,ARGS.safety,ARGS.pearson,ARGS.render)

    if ARGS.mode == "FLEXRSCORE":
        from src.em import flexrscore
        from src.em.flexrscore import flexrscoremain
        bestscore, FS, peakheight, peakprominence = flexrscoremain(ARGS.filename,\
                        ARGS.sigmathreshold,ARGS.plot,ARGS.peakheight,\
                        ARGS.peakprominence,ARGS.peakwidth,ARGS.peakdistance,ARGS.step,\
                        ARGS.mode,ARGS.geotolerance,ARGS.build_limit,ARGS.ligand,\
                        ARGS.pdb_in,ARGS.distance,ARGS.singleconfs)

    if (FS):
        from src.flexrpkg import flexr_find
        from src.flexrpkg.flexr_find import flexrfindmain
        build_list = flexrfindmain(ARGS.filename,bestscore,ARGS.plot,peakheight,\
                    peakprominence,ARGS.peakwidth,ARGS.peakdistance,ARGS.step,\
                    ARGS.mode,ARGS.geotolerance,ARGS.build_limit,ARGS.ligand,\
                    ARGS.pdb_in,ARGS.distance,ARGS.singleconfs)
        # needs building bit here too

    if (ARGS.mode == "FLEXR"):
        from src.flexrpkg import flexr_find
        from src.flexrpkg.flexr_find import flexrfindmain
        build_list = flexrfindmain(ARGS.filename,ARGS.sigmathreshold,ARGS.plot,ARGS.peakheight,\
                    ARGS.peakprominence,ARGS.peakwidth,ARGS.peakdistance,ARGS.step,\
                    ARGS.mode,ARGS.geotolerance,ARGS.build_limit,ARGS.ligand,\
                    ARGS.pdb_in,ARGS.distance,ARGS.singleconfs)
        if (ARGS.build == 'True') & (ARGS.pdb_in is not None):
            from src.flexrpkg import top_level
            from src.flexrpkg.top_level import get_coot_loc
            libraryloc, cootloc,cootexe = get_coot_loc()
            try:
                buildpath = cootloc+'src/building/flexr_build.py'
                os.system('%s --script %s %s %s %s %s %s' % (cootexe,buildpath,build_list,ARGS.pdb_in,ARGS.branching,ARGS.cootmolnum,ARGS.exitcoot))
            #os.system('COOT1 --script flexr_build.py %s %s %s %s %s' % (build_list,ARGS.pdb_in,ARGS.branching,ARGS.cootmolnum,ARGS.exitcoot))
            except:
                print('Cannot find flexr building script.')
        else:
            print('Done.')

if __name__ == '__main__':
    from src.flexrpkg.top_level import args
    ARGS = args()
    main(ARGS)
    print('FLEXR is finished.')
    print('')








