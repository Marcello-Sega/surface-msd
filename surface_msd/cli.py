def main():
    import argparse
    from surface_msd.msd_core import LayerMSD
    
    
    dsc = """
           Compute mean square displacement of the center of
           mass of residues within the surface molecular layer.
           For the MSD calculation to make sense,
           configuration in the trajectory have to be stored
           every integration timestep.
    """
    
    parser = argparse.ArgumentParser(description=dsc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    """ Adding arguments """
    parser.add_argument ('config', metavar = 'configuration',
                                 type=str, help    = 'Input (.gro)' )
    parser.add_argument ('traj'  , metavar = 'trajectory'   ,
                                 type=str, help    = 'Input (.trr or .xtc)' )
    parser.add_argument ('residue' , metavar = 'residue'   ,
                                 type=str, help    = 'Residue type to compute the MSD of' )
    
    """ Adding optional arguments """
    
    parser.add_argument ('--output', type=str, default = 'msd.dat',
                         help    = 'Calculated MSD (.dat)' )
    
    parser.add_argument ('--algorithm', type=str, default = 'ITIM', 
                         help    = 'choose the algorithm to use' )
    
    
    parser.add_argument ('--pytim-select', type=str, default = 'name OW',
                         help    = "pass the string used to select the group for the surface analysis (as in, e.g., group=u.select_atoms('name OW')" )
    
    parser.add_argument ('--pytim-params', type=str, default = '',
                         help    = 'pass the comma-separater list of your parameters for the algorithm of choice (string value do not need to be enclosed with quoutation)' )
    
    parser.add_argument('--direction' ,type=str, default='xy', metavar='DIR',
                         help = "Calculate the MSD of these components only (e.g., if normal points along z, 'xy' computes the in-plane MSD, 'z' the out-of-plane MSD, and  'xyz' is the averaged MSD in 3D). Change according to your normal direction.")
    parser.add_argument('--length' ,type=int, default=100,metavar='L',
                         help = 'Size of the window over which the MSD is sampled, in timesteps')
    parser.add_argument ('--start', type=int, default=0,
                         help    = 'First frame to use')
    parser.add_argument ('--stop', type=int, default=-1,
                         help    = 'Last frame to use')
    parser.add_argument ('--stride', type=int, default=1,
                         help    = 'Stride between frames (typically as low as possible)')
    parser.add_argument ('--verbose', action='store_true',
                         help    = 'add some info to the output')
    # internal check
    parser.add_argument ('--test', action='store_true', help = argparse.SUPPRESS )
    
    args = parser.parse_args()
    
    params = {}
    if args.pytim_params.strip():
        for item in args.pytim_params.replace(" ", "").split(","):
            k, v = item.split("=")
            try: params[k] = int(v)
            except:
                try: params[k] = float(v)
                except: params[k] = v.replace("'", "").replace('"', "")
    
    
    if args.verbose: print('loading packages...')
    import numpy as np
    import MDAnalysis as mda
    import pytim
    from pytim.observables import Observable, Position, Mass
    
    def tqdm(x): return x
    try: 
        if args.verbose: 
            from tqdm import tqdm
    except: 
        pass
    
    if args.verbose: print('...done.')
    
    if args.test: 
        from pytim.datafiles import WATER_XTC, WATER_GRO
        args.residue = 'SOL'
        args.config  = WATER_GRO
        args.traj    = WATER_XTC
        
    
    u = mda.Universe(args.config,args.traj)
    refgroup = u.select_atoms('resname '+args.residue)
    MSD = LayerMSD(refgroup, direction = args.direction, npoints=args.length)
    
    if args.algorithm == 'ITIM': alg = pytim.ITIM
    elif args.algorithm == 'GITIM': alg = pytim.GITIM
    
    inter = alg(u, group=u.select_atoms(args.pytim_select), **params)
    
    
    if args.test:
        args.start, args.stride, args.stop =0,1,10
    
    time = np.round(u.trajectory.dt,6) * args.stride * np.arange(args.length+1)
    
    for i,ts in tqdm(enumerate(u.trajectory[args.start:args.stop:args.stride])):
        inp = inter.atoms
        if args.test:
            inp = u.residues[:1].atoms
            inp.atoms.positions =  1.*i
            if i == 3 or i == 4 : 
                inp = u.residues[1:2].atoms
                inp.atoms.positions =  0
        MSD.sample(inp)
    np.savetxt(args.output,np.vstack([time, MSD.data]).T)
